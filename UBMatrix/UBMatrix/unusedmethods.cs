using System;
using System.Collections.Generic;
using System.Text;

namespace DAS.UBMatrix
{
    class UnusedUBMatrixFunctions
    {
        /// <summary>
        /// [It is not used by any client/user function now]
        /// 
        /// Generate UB Matrix by input unit cell parameter and some reflections
        /// Refinement (including linear least square fitting and nonlinear (in future)
        /// will be applied if the number of input reflections are larger than 2. 
        /// </summary>
        ///<source>UnitCell.for: SUBROUTINE GenUB(GI,NCOUNT,TOL)</source> 
        /// <param name="IncidentAngles"></param>
        /// <returns></returns>
        public bool GenerateUBFromLatticeRefine(MotorIncidentAngles[] IncidentAngles) 
        {
#if ENABLED
            UBStorage mStorage = new UBStorage(IncidentAngles);

            // 1. Storage
            Matrix mUBBest = new Matrix(3);
            Matrix mInvUBBest = new Matrix(3);

            // 2. Find all possible combinations of 2 reflections                        
            Vector vH1 = new Vector(3);
            Vector vH2 = new Vector(3);
            Vector uH1 = new Vector(3);
            Vector uH2 = new Vector(3);
            int ncount = IncidentAngles.Length;

            MillerIndices[] calHKLList = new MillerIndices[ncount];
            double[] diff2ThetaList = new double[ncount];
            double[] diffUPhiList = new double[ncount];

            // 2.1 Set up storage
            int nfit = 0;

            // 2.2 Brute force looping starts
            bool continuesearch = true;
            for (int n1 = 0; n1 < ncount - 1; n1++)
            {
                // old: ICalcH(IncidentAngles[n1], vH1);
                MotorIncidentAngles.UnitVectorFromAngle(IncidentAngles[n1], vH1);
                vH1.UnitVector(uH1);

                if (false)
                {
                    Console.WriteLine("(_DB)   Incident Angle 1 and u_phi1:");
                    Console.WriteLine("{0} | {1}\n", IncidentAngles[n1].ToString(), uH1.ToString());

                    m_log.Record(String.Format("Incident Angle 1", IncidentAngles[n1].ToString()));
                    m_log.Record(String.Format("u_phi 1", uH1.ToString()));
                }


                for (int n2 = n1 + 1; n2 < ncount; n2++)
                {
                    // old: ICalcH(IncidentAngles[n2], vH2);
                    MotorIncidentAngles.UnitVectorFromAngle(IncidentAngles[n2], vH2);
                    vH2.UnitVector(uH2);

                    // a) measure angle between these two reflections
                    double dot = Vector.DotProduct(uH1, uH2);
                    double ang12 = DasMath.AcosD(dot);

                    if (ang12 <= 20.0)
                    {
                        // if angle is too small, continue to try next reflection (n2)                        
                        continue;
                    }

                    if (false)
                    {
                        Console.WriteLine("\t(_DB)   Incident Angle 2 and u_phi2:");
                        Console.WriteLine("\t{0} | {1}", IncidentAngles[n2].ToString(), uH2.ToString());
                        Console.WriteLine("\t(u_phi1, u_phi2) = {0} degree", ang12);

                        m_log.Record(String.Format("Incident Angle 2", IncidentAngles[n2].ToString()));
                        m_log.Record(String.Format("u_phi 2", uH2.ToString()));
                        m_log.Record(String.Format("Angle between u_phi 1 and u_phi 2", Convert.ToString(ang12)));
                    }

                    // b)  Find Max. Theta (of n1 and n2), if angle >= 20.0
                    // Note: anglst(1, x) is the array for 2theta of all angles
                    double tmax = 0.0;
                    if (Math.Abs(IncidentAngles[n1].twotheta) > Math.Abs(IncidentAngles[n2].twotheta))
                    {
                        tmax = Math.Abs(IncidentAngles[n1].twotheta);
                    }
                    else
                    {
                        tmax = Math.Abs(IncidentAngles[n2].twotheta);
                    }
                    // Error check

                    if (tmax < 1.0E-16)
                    {
                        string msg = "Error!  Maximum 2Theta is ZERO.  Input Reflections Error! Job aborted.";
                        Console.WriteLine(msg);
                        return false;
                    }
                    tmax = tmax / 2.0 + 1.0;

                    // c) Calculate HKL limits from tmax
                    int[] Lim = new int[3];
                    for (int i = 0; i < 3; i++)
                    {
                        Lim[i] = Convert.ToInt32((2.0 * DasMath.SinD(tmax)) / (mWaveLength * Math.Sqrt(mGI.M[i][i])) + 1.0);
                    }

                    Console.WriteLine("\t(_DB)  Limit of HKL = ({0}, {1}, {2})", Lim[0], Lim[1], Lim[2]);
                    Logging("Limit of (HKL) of n1 and n2",
                        String.Format("{0}, {1}, {2})", Lim[0], Lim[1], Lim[2]));

                    // d) FIND HKL FOR FIRST REFLECTION 
                    //    Algorithm: it is go over all the (h, k, l) and 
                    //    calculate 2theta from (hkl) and compare to experimental 2theta
                    List<MillerIndices> hkllist = GenerateReflectionListInverse(Lim);
                    int sizehklist = hkllist.Count;
                    for (int ihkl = 0; ihkl < sizehklist; ihkl++)
                    {
                        MillerIndices hkl1 = hkllist[ihkl];
                        double d1;
                        double twotheta_cal_1 = Calculate2ThetaFromHKL(hkl1, out d1); //Calculate2Theta(hkl1, mGI, wavelength, out d1);

                        if (Math.Abs(IncidentAngles[n1].twotheta - twotheta_cal_1) <
                            m_arrTolerance[0])
                        {

                            // ii. FIND HKL FOR SECOND REFLECTION 
                            for (int jhkl = 0; jhkl < sizehklist; jhkl++)
                            {
                                MillerIndices hkl2 = hkllist[jhkl];
                                double d2;
                                double twotheta_cal_2 = Calculate2ThetaFromHKL(hkl2, out d2);

                                if (Math.Abs(IncidentAngles[n2].twotheta - twotheta_cal_2) <
                                    m_arrTolerance[0])
                                {
                                    // Match in 2theta for 2nd reflection. **      

                                    // iii. Check: Is the angle between the relections correct?
                                    double da1 = Vector.DotProduct(hkl1, mGI.MultVector(hkl2));
                                    double da2 = DasMath.AcosD(da1 / Math.Sqrt(d1 * d2));
                                    if (Math.Abs(da2 - ang12) < m_arrTolerance[1])
                                    {   // ***
                                        // then ! yes
                                        // trial indices and calculated 2 theta for second reflection
                                        // means: ALL TESTS ARE PASSED!!!

#if _DBOUTPUT
                                        Console.WriteLine("\t(_DB)  (HKL)_1 = {0}, calculated 2theta = {1}",
                                            hkl1.ToString(), twotheta_cal_1);
                                        Console.WriteLine("\t(_DB)  (HKL)_2 = {0}, calculated 2theta = {1}",
                                            hkl2.ToString(), twotheta_cal_2);
                                        Logging("(HKL)_1", hkl1.ToString());
                                        Logging("2theta_1", Convert.ToString(twotheta_cal_1));
                                        Logging("(HKL)_2", hkl2.ToString());
                                        Logging("2theta_2", Convert.ToString(twotheta_cal_2));
#endif

                                        // rate 'apparent indices' for rest of reflections
                                        int ierr;
                                        // old: MakeUB(mB, this, hkl1, hkl2, IncidentAngles[n1], IncidentAngles[n2], out ierr);
#if true
                                        _MakeUB(hkl1, hkl2, IncidentAngles[n1], IncidentAngles[n2], out ierr);
                                        if (ierr == 1)
                                        {
                                            string msg = "Error In Generating UB Matrix!  Quit GenUB()";
                                            Console.WriteLine(msg);
                                            return false;
                                        }

#else
                                        Reflection r1, r2;
                                        r1 = new Reflection();
                                        r1.motor_angles = IncidentAngles[n1];
                                        r1.rlvector = hkl1;

                                        r2 = new Reflection();
                                        r2.motor_angles = IncidentAngles[n2];
                                        r2.rlvector = hkl2;

                                        CalUBMatrix(r1, r2);
                                        // Console.WriteLine("UB 2:");
                                        // Console.WriteLine(this.ToString());
#endif

                                        double det;
                                        this.InvertMatrix(m_invertUB, out det);

                                        int nfound;
                                        // old: TestRef(mTolerance, mUB, mInvUB, mGI, IncidentAngles, out nfound);
                                        TestUBByReflections(IncidentAngles, calHKLList, diff2ThetaList, diffUPhiList, out nfound);

                                        int[] anglepair = new int[2];
                                        anglepair[0] = n1;
                                        anglepair[1] = n2;

                                        int overlimit = 0;
                                        for (int ifc = 0; ifc < nfound; ifc++)
                                            if (diff2ThetaList[ifc] > m_arrTolerance[0])
                                                overlimit++;

                                        if (nfound >= 2 && overlimit < nfound - 1)
                                            mStorage.AddUBSuite(this, anglepair, calHKLList, diff2ThetaList, diffUPhiList);

                                        if (nfound > nfit)
                                        {
                                            // 1. store the result to some temporary solution
                                            this.CopyToMatrix(mUBBest);
                                            m_invertUB.CopyToMatrix(mInvUBBest);
                                            // mUB.CopyToMatrix(UBMatrix);
                                            // mInvUB.CopyToMatrix(InvUBMatrix);
                                            nfit = nfound;
                                        } // endif 
#if true
                                        Console.WriteLine("----- Summary -----  nfound = {0}", nfound);
#endif
                                        if (nfit == ncount)
                                        {
                                            continuesearch = false;
                                            break;
                                        }
                                    } // endif ***
                                } // endif **
                            } // for: jhkl
                        } // if math... ??? enddo *
                        if (continuesearch == false)
                            break;
                    }  // for: ihkl
                    if (continuesearch == false)
                        break;
                } // 90  enddo  ! Loop on n2
                if (continuesearch == false)
                    break;
            } // enddo  ! loop on n1

            // 3. Check Result and set to the best UB matrix 
            if (nfit == 0)
            {
                Console.WriteLine("No UB Matrix Found");
                return false;
            }
            else
            {
                mUBBest.CopyToMatrix(this);
                mInvUBBest.CopyToMatrix(m_invertUB);

                Console.WriteLine("Best UB Can Fit {0} Reflections", nfit);
                Logging("Number of Reflections Fit By Best UB", Convert.ToString(nfit));

            }

            // 4. Optional Output Result
            for (int icount = 0; icount < ncount; icount++)
            {
                // a) Incident angle --> (HKL)
#if false
                Vector vH = new Vector(3);
                ICalcH(IncidentAngles[icount], vH);
                double dexp = wavelength / (2.0 * DasMath.SinD(0.5 * IncidentAngles[icount].twotheta));
                Vector vHKL = UBMatrix.MultVector(vH);
                for (int i = 0; i < 3; i++)
                {
                    vHKL.V[i] = DasMath.NearestInteger(vHKL.V[i] / dexp);
                }
#else
                Vector vH = new Vector(3);
                MillerIndices vHKL = new MillerIndices();
                CalculateHKLFromIncidentAngle(IncidentAngles[icount], vHKL, vH);
#endif

                // b) (HKL) --> 2theta
#if false
                double dcal = Vector.DotProduct(vHKL, mGI.MultVector(vHKL));
                double twothetacal = 2.0 * DasMath.AsinD(wavelength * Math.Sqrt(dcal) / 2.0);
#else
                double dcal;
                double twothetacal = Calculate2ThetaFromHKL(vHKL, out dcal);
#endif
                // c) (HKL) --> IncidentAngle
                Vector vC = this.MultVector(vHKL);
                Vector uC = vC.UnitVector();
                Vector uH = vH.UnitVector();
                double dot = Vector.DotProduct(uC, uH);
                if (dot >= 1.0)
                    dot = 1.0;
                // e) cal diff
                double diff2theta = twothetacal - IncidentAngles[icount].twotheta;
                double diffeuler = DasMath.AcosD(dot);
                // f) print out and log
                string title1 = String.Format("Reflection {0}", icount);
                string value1 = String.Format("{0}  --> Caculated HKL = {1}", IncidentAngles.ToString(), vHKL.ToString());
                string title2 = "Difference in 2Theta and Euler";
                string value2 = String.Format("{0},  {1}", diff2theta, diffeuler);

                Console.WriteLine(String.Format("{0}:  {1}", title1, value1));
                Console.WriteLine(String.Format("{0}:  {1}", title2, value2));
            }

            mStorage.Stat();
            mStorage.Print("allub.txt");
#endif
            return true;
        } // GenUB 


        // FIXME - Need to justify physics
        ///<source>SUBROUTINE PARAMS</source>
        ///<VZ>Guess XTX = tild(UB) * (UB) = G^-1</VZ>
        public void Params(Matrix mX, UnitCell cA1, UnitCell cA2, double[] angles)
        {
            // EXTRACTS 
            // (1) DIRECT LATTICE PARAMETERS A1 AND VOL1 AND 
            // (2) RECIPROCAL LATTICE PARAMETERS A2 AND VOL2 
            // FROM ORIENTATION MATRIX X = UB
            // 	DIMENSION XTX(3,3)
            // 	Global: common /clsq/x(3,3),a1(9),a2(6),vol1,vol2,xa(15,16),xb(12)
            // 	DATA RAD/57.2958/
            // C  WORK RECIPROCAL LATTICE FROM METRIC

            throw new Exception("The codes below will be unified ... ");

            //Matrix mXtx = new Matrix(3);

            //// FIXME - ? XTX = X x X
            //for (int i = 0; i < 3; i++)
            //{
            //    for (int j = 0; j < 3; j++)
            //    {
            //        mXtx.M[i][j] = 0.0;
            //        for (int k = 0; k < 3; k++)
            //        {
            //            mXtx.M[i][j] += mX.M[k][i] * mX.M[k][j];
            //        }
            //    } // DO 40 J=1,3
            //}

            //// 2. Form (reciprocal?) unit cell
            //for (int i = 0; i < 3; i++)
            //{
            //    cA2.Axiles[i] = Math.Sqrt(mXtx.M[i][i]);
            //} // DO 50 I=1,3

            //cA2.alpha = mXtx.M[1][2] / (cA2.b * cA2.c);
            //cA2.beta = mXtx.M[2][0] / (cA2.c * cA2.a);
            //cA2.gamma = mXtx.M[0][1] / (cA2.a * cA2.b);

            //cA2.CalculateReciprocalUnitCell(cA1);

            //// recip(a2,vol2,a1,vol1);  FIXME:  what is this for???
            //for (int i = 0; i < 3; i++)
            //{
            //    double tx = Math.Sqrt(1.0 - cA1.Angles[i] * cA1.Angles[i]);
            //    double ty = cA1.Angles[i];
            //    angles[i] = DasMath.rad * DasMath.Atan2(tx, ty);
            //}
            //return;
        }


        // FIXME - Need to know physics
        // Relocate - Put to class Matrix
        ///<source>SUBROUTINE MXLNEQ(A,NRANK,IDA,DET,JRANK,EPS,NTEMP,NCV)</source>
        public void MxLNeq(Matrix mA, int nrank, int ida, int ncv, int[] ntemp, out double det,
            out int jrank, out int eps)
        {
            // 
            //  MATRIX INVERSION AND EQUATIONS SOLUTION ROUTINE FROM R. L. HOTCHKISS
            //     UNIVERSITY OF MINNESOTA.
            // 
            //  ROUTINE USES METHOD OF GAUSS-JORDAN DOUBLE PIVOTING.
            // 
            //  - A IS MATRIX OF LINEAR SYSTEM WITH COEFFICIENTS AUGMENTED BY VECTOR
            //  - NRANK IS INPUT ORDER OF SYSTEM
            //  - IDA IS FIRST DIMENSION OF A IN CALLING PROGRAM
            //  - DET IS OUTPUT VALUE OF DETERMINANT
            //  - JRANK IS OUTPUT RANK OF SYSTEM
            //  - EPS IS TEST VALUE USED TO DETERMINE SINGULARITY
            //  - NTEMP IS INTEGER ARRAY OF SIZE NRANK
            //  - NCV IS NUMBER OF COLUMNS AUGMENTED. IF NCV=0, INVERSION ONLY.
            // 
            // 	DIMENSION A(2),NTEMP(2)				! CHANGE FEB 2003 TO FIX
            // 	DIMENSION A(256),NTEMP(NRANK)		! OUT-OF-BOUNDS IN DEBUG EXE
            //
            // INITIALIZE
            throw new Exception("Will be implemented in Phase-2");
//            int nc = nrank + ncv;
//            det = 1.0;
//            double[] dA = new double[256];
//            int ij, kk;
//            double pivot;
//#if true
//            ij = 0;
//            kk = 0;
//            pivot = 0;
//            Console.WriteLine("evaluating to ij, kk, pivot is not true");
//#endif

//            //  START MASTER REDUCTION	LOOP
//            int m, l, jsub, lj, kj;
//            double temp;
//            for (int k = 0; k < nrank; k++)
//            {
//                // DO 240 K=1,NRANK
//                double piv = 0.0;
//                int ksub = ida * (k - 1);
//                // START SEARCH FOR PIVOT
//                for (int j = k; j < nrank; j++)
//                { // 110
//                    jsub = ida * (j - 1);
//                    for (int i = k; i < nrank; i++)
//                    { // do 110 
//                        ij = jsub + i;
//                        // check if abs(a(i,j)).gt.abs(pivot)
//                        double p = Math.Abs(mA.M[i][j]);
//                        if ((p - piv) > 0)
//                        {
//                            // found new value for trial pivot 
//                            piv = p;
//                            pivot = mA.M[i][j];
//                            m = j;
//                            l = i;
//                        }
//                    }
//                } // 110 
//                //  CHECK SIZE OF PIVOT
//#if true
//                Console.WriteLine("The evaluating to eps and pivot is wrong!");
//                eps = 0;
//                pivot = 0;
//                l = 0;
//                m = 0;
//                kj = 0;
//                kk = 0;
//#endif
//                if ((piv - eps) <= 0)
//                {
//                    // 120,130,130	//
//                    // SINGULARITY RETURN 
//                    jrank = k - 1;
//                    eps = 0;
//                    return;
//                }
//                else
//                {
//                    //  VALID PIVOT FOUND, UPDATE DETERMINANT
//                    det = det * pivot;
//                    //  SAVE ORIGINAL ROW AND COLUMN OF PIVOT 
//                    ntemp[k] = 256 * l + m;
//                }
//                if ((l - k) != 0)
//                { // 140,160,140
//                    //  NEED TO PUT PIVOT ELEMENT INTO	A(K,K),	INTERCHANGE ROWS
//                    det = -det;
//                    for (int j = 0; j < nc; j++)
//                    { // DO 150 J=1,NC
//                        jsub = ida * (j - 1);
//                        lj = l + jsub;
//                        kj = k + jsub;
//                        temp = mA.M[l][j];
//                        mA.M[l][j] = mA.M[k][j];
//                        mA.M[k][j] = temp;
//                    } //150 

//                }
//                // 160 
//                if ((m - k) != 0)
//                {  // 170,190,170
//                    //  interchange columns
//                    det = -det;
//                    int jm = ida * (m - 1);
//                    int jk = ksub;
//                    for (int j = 0; j < nrank; j++)
//                    { // do 180 j
//                        jm = jm + 1;
//                        jk = jk + 1;
//                        temp = mA.M[j][m];
//                        mA.M[j][m] = mA.M[j][k];
//                        mA.M[j][k] = temp;
//                    } //180 

//                    //  reduce	pivot row
//                } // 190 
//                pivot = 1.0 / pivot;
//                for (int j = 0; j < nc; j++)
//                { // do 200 j
//                    kj = ida * (j - 1) + k;
//                    mA.M[k][j] = pivot * mA.M[k][j];
//                } // 200 

//                kk = k + ksub;
//                mA.M[k][k] = 0.0;
//                //  REDUCE	NON-PIVOT ROWS
//                for (int i = 0; i < nrank; i++)
//                { // DO 230 I=1,NRANK
//                    int ik = i + ksub;
//                    temp = mA.M[i][k];
//                    if ((i - k) != 0)
//                    { // 210,230,210
//                        // 210 
//                        for (int j = 0; j < nc; j++)
//                        { // DO 220 J
//                            ij = ida * (j - 1) + i;
//                            kj = ida * (j - 1) + k;
//                        } // 220 
//                        dA[ij] = dA[ij] - temp * dA[kj];
//                        dA[ik] = -temp * pivot;
//                    }
//                } // 230
//            } // 	240 
//            dA[kk] = pivot;

//            //  INVERSION COMPLETE, RESTORE ORIGINAL ROW AND COLUMN ORDER
//            jrank = nrank;
//            int kx = nrank;
//            for (int j = 0; j < nrank; j++)
//            { // do 300 j
//                int ksub = ida * (kx - 1);
//                m = ntemp[kx] / 256;
//                // l=mod(ntemp(kx),256);
//                l = ntemp[kx] % 256;
//                if ((m - kx) != 0)
//                { // 250,270,250 //
//                    // INTERCHANGE ROWS 
//                    int msub = ida * (m - 1);
//                    for (int i = 0; i < nrank; i++)
//                    { // DO 260 I
//                        int ik = ksub + i;
//                        int im = msub + i;
//                        temp = dA[ik];
//                        dA[ik] = dA[im];
//                        dA[im] = temp;
//                    } // 260 
//                } // 270 
//                if ((l - kx) != 0)
//                { // 280,300,280
//                    // INTERCHANGE COLUMNS
//                    for (int i = 0; i < nc; i++)
//                    { //280 DO 290
//                        int ki = ida * (i - 1) + kx;
//                        int li = ida * (i - 1) + l;
//                        temp = dA[ki];
//                        dA[ki] = dA[li];
//                        dA[li] = temp;
//                    } // 290 
//                }
//                kx = kx - 1;
//            } // 300 

//            eps = 0;
//#if true
//            Console.WriteLine("eps = 0 is faked");
//#endif
//            return;

        }

#if false

	///<source>SUBROUTINE CONSTR(CRYSTP,IERR,NPARM)</source> 
        public void Constrain(int crystp, out int ierr, out int nparm, Matrix mXA, Matrix mXAI, Matrix mB, double[] a2,
            int[] ncnst, int[][] cnstyp, int[] ncoef)
        {

            // mXA:  No need for input... will be intialized to zero at the start of subroutine

            //C
            //C SUBROUTINE TO APPLY CONSTRAINTS FOR CRYSTAL OF 'CRYSTP' TO NORMAL
            //C  EQUATIONS MATRIX XA, SOLVE SYSTEM AND RETURN XA INVERSE IN
            //C  XA, SOLUTION IN XB. METHOD OF LAGRANGIAN MULTIPLIERS USED
            //C  FOR CONSTRAINTS. IERR = 0 IF NO ERROR, NON-ZERO MEANS ERROR.
            //C
            //	DIMENSION NTEMP(14)
            //	common /clsq/x(3,3),a1(9),a2(6),vol1,vol2,xa(15,16),xb(12)
            //	INTEGER CRYSTP,NCNST(6),CNSTYP(5,6),NCOEF(10)
            //	1,ROWTP(3,10),COEFTP(3,10)
            //	REAL COEF(15)
            //C
            //C CRYSTP = CRYSTAL TYPE: 0-TRICLINIC, 1-MONOCLINIC(B-UNIQUE)
            //C    2-MONOCLINIC (C-UNIQUE), 3-ORTHORHOMBIC, 4-TETRAGONAL
            //C    5-HEXAGONAL, 6-CUBIC
            //C
            //C NCNST(J) = NUMBER OF CONSTRAINTS FOR CRYSTAL TYPE J   
            //C
            //C CNSTYP(J,I) = CONSTRAINT TYPE FOR J'TH CONSTRAINT OF   
            //C     I'TH CRYSTAL TYPE
            //C
            //C NCOEF(J) = NUMBER OF COEFFICIENTS IN CONSTRAINT EQUATION I   
            //C
            //C ROWTP(J,I) = ROW (COLUMN) NUMBER IN NORMAL EQUATIONS MATRIX
            //C     OF J'TH PARAMETER IN I'TH CONSTRAINT EQUATION
            //C
            //C COEFTP(J,I) = COEFFICIENT TYPE OF J'TH TERM IN I'TH
            //C     CONSTRAINT EQUATION
            //C
            //	DATA NCNST/2,2,3,4,4,5/         
            //	DATA CNSTYP/9,10,3*0,7,8,3*0,1,2,3,2*0,4,1,2,3,0	
            //	1,6,1,7,8,0,4,5,1,2,3/
            //	DATA NCOEF/2,2,2,2,2,3,3,3,3,3/
            //	DATA ROWTP/2,4,0,3,7,0,6,8,0,1,5,0,1,9,0
            //	1,1,4,5,7,3,6,7,8,6,6,8,2,6,4,2/
            //	DATA COEFTP/1,2,15,3,4,15,5,6,15,7,8,15,7,8,15
            //	1,7,7,8,4,3,9,10,6,11,5,6,12,13,2,14/
            //C
            //C EVALUATE MULTIPLIER COEFFICIENTS
            //C
            //
            double[] coef = new double[15];
            coef[0] = a2[0] / a2[1];
            coef[1] = a2[1] / a2[0];
            coef[2] = a2[0] / a2[2];
            coef[3] = a2[2] / a2[0];
            coef[4] = a2[1] / a2[2];
            coef[5] = a2[2] / a2[1];
            coef[6] = 1.0;
            coef[7] = -1.0;
            coef[8] = a2[1] * a2[7] / a2[2];
            coef[9] = -a2[2] * a2[5] / a2[0];
            coef[10] = a2[1] * (1.0 - a2[5] * a2[5]) / a2[2];
            coef[11] = a2[0] * a2[4] / a2[1];
            coef[12] = -a2[1] * a2[4] / a2[2];
            coef[13] = a2[0] * (1.0 - a2[4] * a2[4]) / a2[1];
            coef[14] = 0.0;

            // CLEAR CONSTRAINT PORTION OF MATRIX XA
            //
            for (int j = 9; j < 16; j++)
            {
                for (int i = 0; i < 9; i++)
                {
                    mXA.M[i][j] = 0.0;
                }
            }
            for (int j = 0; j < 16; j++)
            {
                for (int i = 9; i < 15; i++)
                {
                    mXA.M[i][j] = 0.0;
                }
            }

            // Constraint 
            if (crystp >= 1 && crystp <= 6)
            {

                // GET NUMBER OF CONSTRAINT EQUATIONS
                int neq = ncnst[crystp];
                nparm = 9 - neq;

                // C LOOP THROUGH CONSTRAINTS
                for (int i = 0; i < neq; i++)
                {

                    // GET CONSTRAINT TYPE AND NUMBER OF COEFFICIENTS 
                    int ityp = cnstyp[i][crystp];
                    int numb = ncoef[ityp];

                    // LOOP THROUGH TERMS IN CONSTRAINT EQUATION
                    int irow;

                    for (int j = 0; j < numb; j++)
                    {

                        // GET ROW AND COLUMN OF TERM TO MODIFY AND COEFFICIENT OF MODIFICATION
                        irow = rowtp(j, ityp); // FIXME - should minus 1?
                        ncf = coeftp(j, ityp); // FIXME - should minus 1?
                        mXA[irow][8 + i] = coef(ncf);
                        mXA[8 + i][irow] = coef(ncf);
                    } // 120
                } // 130

                nrank = 9 + neq; // FIXME - should minus 1?

                // COPY VECTOR TO AUGMENTED POSITION, INVERT MATRIX AND SOLVE EQUATIONS
                for (int i = 0; i < 9; i++)
                    mXA(i, nrank + 1) = xb(i);	// 140

                Mxlneq(mXA, nrank, 15, out det, jrank, 1.0e-6, ntemp, 1);
                if (jrank == nrank)
                {
                    for (int i = 0; i < 9; i++)
                    {
                        vXB.V[i] = mXA.M[i][nrank + 1];
                    }
                    ierr = 0;
                    return;
                }
                else
                {
                    //go to 200
                    //
                    // MATRIX SINGULAR
                    ierr = 1;
                    return;
                }

            }
            else
            {
                // CRYSTAL TYPE NOT PROPER
                // goto 210
                ierr = 2;
                return;
            }

        }

#endif

        ///<source>subroutine lsvec(ICRSTP)</source>
        public bool RefineUnitCellLsq(Reflection[] reflectslist, string celltype)
        {
#if false
            //
            // reads reflection data from reflection list and performs
            // a constrained least-squares refinement of the
            // unit cell parameters, writing full results to the screen.
            // coded 6-May-81 by Russell Ralph and Larry Finger
            // Geophysical Laboratory, Washington DC
            //
            // method from Shoemaker and Bassi (1970), Acta Cryst. A26, 97.
            //
            // Modified for Diffractometer Control System by LWF,  6-Feb-92
            // MODIFIED FOR GUI CALLBACK STRUCTURE, RJA JAN 2004

            // include 'single.inc'
            // dimension a(3,3),yb(3,3),ss(3),vect(3,100),aa1(6),as1(6)
            // 1 ,bb(3,3),u(3,3),v(3,3),v0(3,3),v0inv(3,3),d(3,3),angl(4)
            // 2 ,anglc(4),hh(3,100),atemp(7,9),
            // 3 dcov(9,9),derivs(7,9),sigs(7),tmp1(3),tmp2(3)
            // integer h(3,100),wk3(12),ISEQ(listlen)
            // common /clsq/x(3,3),a1(9),a2(6),vol1,vol2,xa(15,16),xb(12)
            // CHARACTER*80 TEXT
            // DATA RAD/57.2958/,FFACT/5.E-4/
            //
            // variable name table:
            // vect --> uphilist

            // 1. INITIALIZE COUNTERS AND BEGIN UNCONSTRAINED PROBLEM
            Matrix mA = new Matrix(3);
            Matrix mYB = new Matrix(3);
            Matrix mX = new Matrix(3);
            mA.SetToZero();
            mYB.SetToZero();
            mX.SetToZero();

            Matrix mXA = new Matrix(12);
            mXA.SetToZero();
            Vector vXB = new Vector(12);
            vXB.SetToZero();

            // 2.   Transfer (1) Incident angles and (2) Miller indices 
            //      from input list to matrices for 
            //      Lease Square Fitting (linear)
            int nlst = reflectslist.Length;
            Vector[] hphis = new Vector[nlst];
            MillerIndices[] hkls = new MillerIndices[nlst];

            for (int i = 0; i < nlst; i++)
            {
                hkls[i] = new MillerIndices();
                for (int j = 0; j < 3; j++)
                {
                    hkls[i].V[j] = DasMath.NearestInteger(hkls[i].V[j]);
                }
                hphis[i] = new Vector(3);
                MotorIncidentAngles.VectorFromAngles(reflectslist[i].MotorAngles, mWaveLength, hphis[i]);
            }

            // 3. Do initial unconstrained lsq form normal equations matrices a, yb, and xa
            for (int k = 0; k < nlst; k++)
            {
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        mA.M[j][i] += hphis[k].V[j] * hphis[k].V[j];  // H.M[i][k] * H.M[j][k];
                        mYB.M[j][i] += hkls[k].V[i] * hkls[k].V[j];
                    }
                }
            }
            // 
            for (int j = 0; j < 3; j++)
            {
                for (int i = 0; i < 3; i++)
                {
                    mXA.M[i][j] = mA.M[i][j];
                    mXA.M[i + 3][j + 3] = mA.M[i][j];
                    mXA.M[i + 6][j + 6] = mA.M[i][j];
                }
            }
            // 4. INVERT MATRIX A
            Matrix mAI = new Matrix(3);
            double det;
            mA.InvertMatrix(mAI, out det);
            if (det < 1.0E-6)
            {
                string errmsg = "Matrix is Singular";
                Console.WriteLine(errmsg);
                return false;
            }

            // 5. SOLVE FOR ORIENTATION MATRIX X
            double[] vSS = new double[3];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        mX.M[i][j] += mA.M[j][k] * mYB.M[k][i];
                    }
                }
                vSS[i] = 0.0;
            }

            // Get Cell Parameters 
            // <source>call params</source>
            m_directUnitCell.CalculateFromG(mG);

            // 6.  C C  BEGIN CONSTRAINED LEAST-SQUARES  construct the true orientation matrix bb corresponding to the unconstrained cell parameters

            // c C  FORM AND INVERT BB AND OBTAIN THE TRUE ORTHOGONAL ROTATION U FROM X c
            // <source>CALL BCALC(BB,A1,A2)</source>
            Matrix mBI = new Matrix(3);
            _CalculateUnitCellFromUBMatrix();
            mB.InvertMatrix(mBI, out det);
            Matrix.MulMatrix(mX, mBI, mU);

#endif

#if false

	// constrain reciprocal lattice to have proper symmetry alpha never refined in cases treated
	a2[3]=0.0;
	
	switch (icrstp)
	{
		case "MONOCLINIC - B UNIQUE": 
		{
			a2[5] = 0.0;
			break;
		}
		case "MONOCLINIC - C UNIQUE":
		{
			a2[4]=0.0;
			break;
		}
		case "ORTHORHOMBIC":
		{
			a2(4)=0.0; 
			a2(5)=0.0;
			break;
		}
		case "TETRAGONAL":
		{
			a2[0]=Math.Sqrt(a2[0]*a2[1]);
			a2[1]=a2[0];
			a2[4]=0.0;
			a2[4]=0.0;
			break;
		}
		case "HEXAGONAL":
		{
			a2[0] = Math.Sqrt(a2[0]*a2[1]);
			a2[1] = a2[0];
			a2[4] = 0.0;
			a2[5] = 0.5;
			break;
		}
		case "CUBIC":
		{
			a2[0]=(a2[0]*a2[1]*a2[2])**0.33333333;
			a2[1]=a2[0];
			a2[2]=a2[0];
			a2[4]=0.0;
			a2[5]=0.0;
			break;
		}
		case "RHOMBOHEDRAL":
		{
			a2[0] = (a2[0]*a2[1]*a2[2])**0.33333333;
			a2[1] = a2[0];
			a2[2] = a2[0];
			a2[3] = (a2[3]+a2[4]+a2[5])/3.0;
			a2[4] = a2[3];
			a2[5] = a2[3];
			break;
		}
	}
	
	// CONSTRUCT CONSTRAINED B
	recip(a2,vol2,a1,vol1);
	Bcalc(mBI,a1,a2);

	// CONSTRUCT CONSTRAINED GUESS V0 FOR ORIENTATION MATRIX
	Matrix.MulMatrix(mU, mBI, mV0);

	// C  INVERT V0
	mV0.InvertMatrix(mV0Inv, out det);
	if (Math.Abs(det) < 0.00001){
		string errmsg = "Normal Equations Matrix is Singular.";
		Console.WriteLine(errmsg);
		return false;
	}

	// CALCULATE "OBSERVED" INDICES HH
	//
	for (int i = 0; i < nlst; i ++){
		for (int j = 0; j < 3; j ++){
			mHH[j][i] = 0.0;
			for (int k = 0; k < 3; k ++){
				mHH[j][i] += mV0Inv[j][k] * uphilist;
			}
		}
	}

	// SET UP EXPANDED VECTOR XB
	for (int k = 0; k < nlst; k ++){
		for (int i = 0; i < 3; i ++){
	    		mXB[i]   = mXB[i]   + mH[i][k]*mHH[0][k];
	    		mXB[i+3] = mXB[i+3] + mH[i][k]*mHH[1][k];
	    		mXB[i+6] = mXB[i+6] + mH[i][k]*mHH[2][k];
		}
	}

	// TREAT CONSTRAINED PROBLEM
	// CALL CONSTR(ICRSTP,IERR,NPARM)
	Constrain(icrstp, out ierr, out nparam);

	if(ierr != 0){
		string errmsg = "Normal Equations Matrix is Singular.";
		Console.WriteLine(errmsg);
		return false;
	}

	for (int i = 0; i < 3; i ++){
	  mD[0][i]=vXB[i];
	  mD[1][i]=vXB[i+3];
	  mD[2][i]=vXB[i+6];
	}

	// EXTRACT ORIENTATION MATRIX V AND CELL PARAMETERS FROM D
	//
	Matrix.MulMatrix(mV0, mD, mV);
	mV.CopyTo(mX);
	Params();

	//  BACK CALCULATE "OBSERVED" INDICES AND EVALUATE VARIANCE- 
	//   -COVARIANCE MATRIX DCOV FOR D
	double fact=0.0;
	for (int k = 0; k < nlst; k ++){
		for (int i = 0; i < 3; i ++){
			vTmp1[i] = 0.0;
			vTmp2[i] = -mHH[i][k];
			for (int jj = 0; jj < 3; jj ++){
		    		vTmp1[i] += mX[i][jj]*mH[jj][k];
		    		vTmp2[i] += mD[i][jj]*mH[jj][k];
			}
		}
		Calca(uphilist[k], angl, out ierr);
	    	Calca(vTmp1, anglc, out ierr);
		double  temp=0.0;
		for (int jj = 0; jj < 3; jj ++){
			vTmp1[jj] = uphilist[k].V[jj] - mTmp1[jj];
			Temp=Temp+vTmp1[jj]**2;
		}
	    	fact=fact+vTmp2.V[0]**2+ vTmp2[1]**2+vTmp2[2]**2;
	    	temp=Math.Sqrt(temp);
	}

	temp=fact;
	fact=fact/(3*num-nparm);
	mXA.Multiply(fact, mDcov);

	for (int i = 0; i < 3; i++){
        cell[i] = a1[i];
	  cell[i+3] = a1[i+6];
	}
	cell[6] = vol1;

	// PERTURB D ONE ELEMENT AT A TIME, PROPAGATE RESULTS TO
	// ORIENTATION MATRIX V, AND ACCUMULATE PARTIAL DERIVATIVES
	// FOR EACH D.L. PARAMETER IN MATRIX DERIVS
	//
	vA1.CopyTo(vAA1);
	voll1=vol1;

	for (int j = 0; j < 3; j ++){
		for (int i = 0; i < 3; i ++){
			// !only original line of code
			mD.M[i][j] += ffact;
			for (int jj = 0; jj < 3; jj ++){
				for (int ii = 0; ii < 3; ii ++){
					mX[ii][jj] =0.0;
					for (int k = 0; k < 3; k ++){
						mX.M[ii][jj] += mV0.M[ii][k]*mD.M[k][jj];
					}
					Param();
					for (int k = 0; k < 6; k ++){
						vAs1.V[k] = 0.0;
					}
					l=j+3*(i-1);
					for (int k = 0; k < 6; k ++){
						derivs.M[k][l]=(a1[k]-aa1[k])/ffact;
					} // 1530
					derivs.M[6,l] = (vol1-voll1)/ffact;
					//!only original
					mD.M[i][j] = mD.M[i][j]-ffact;
				}
			}
		}
	} // 1690

			// CALCULATE AND PRINT STANDARD DEVIATIONS
			for (int k = 0; k < 7; k ++){
				for (int i = 0; i < 9; i ++){
					mAtemp.M[k][i] = 0.0;
					for (int j = 0; j < 9; j ++){
						mAtemp.M[k][i] += mDerivs.M[k][j] * mDcov[j][i];
					}
				}
			}

			for (int k = 0; k < 7; k ++){
				for (int i = 0; i < 7; i ++){
					mDcov.M[k][i] = 0.0;
					for (int j = 0; j < 9; j ++){
						mDcov.M[k][i] += mAtemp.M[k][j]*mDerivs.M[i][j];
					}
				}
			} // 1710

			// CLEAR NON-REFINED ROWS AND COLUMNS OF DCOV
			for (int i = 0; i < 6; i ++){
				mDcov.M[i,3] = 0.0;
				mDcov.M[3,i] = 0.0;
				if (icrstp != 2){
					mDcov.M[5][i]=0.0;
					mDcov.M[i][5]=0.0;
				}
				// 1712
				if(icrstp != 1){
					mDcov.M[4][i]=0.0;
					mDcov.M[i][4]=0.0;
				}
			} // 1715

			for (int k = 0; k < 7; k ++){
				fact=mDcov.M[k][k];
				if(fact <= 0.0)
					fact=1.0e-24;
				sigs[k] = Math.Sqrt(fact);
			} // 1720

			for (int j = 0; j < 7; j ++)
				for (int k = 0; k < 7; k ++)
					mAtempM[j][k]=mDcov.M[j][k]/(sigs[j]*sigs[k]); // 1725

			for (int k = 3; k < 6; k ++){
				sigs[k]=DasMath.Rad*sigs[k]*Math.Sqrt(1.0-aa1[k]**2);
			} // 1730

			for (int i = 0; i < 7; i ++){
				esdcell[i]=sigs[i];
            }
#endif
            return true;
        }

    } // End of class


#if false
        /// <summary>
        /// Calcualte (Miller) indices corresponding to a (instrument) angle
        /// <source>UBCALS.FOR -> HCALC()</source>
        /// </summary>
        /// <param name="angl">Instrument angle</param>
        /// <returns>Miller indices</returns>
        public MillerIndices CalHc(MotorIncidentAngles angl)
        {
            // 1. If the inversed UB is not up-to-date, Invert it!
            if (mInvUBCounter != mCalUBCounter)
            {
                double det;
                this.InvertMatrix(mInvUB, out det);
                mInvUBCounter = mCalUBCounter;
            }

            // 2. Calculate h_phi from input angle
            Vector h_phi = new Vector(3);
            MotorIncidentAngles.VectorFromAngles(angl, mWaveLength, h_phi);

            // 3. Calculate h
            MillerIndices h = new MillerIndices();
            Matrix.MultVector(mInvUB, h_phi, h);

            return h;
        }

        
        /// <summary>
        /// Calculate UB matrix with 2 reflections and given crystal unit cell parameters 
        /// Matrix B should be calculated in the previous steps
        /// It will calcualte U and UB (this)
        /// ROUTINE TO CALCULATE AN ORIENTATION MATRIX FROM THE UNIT CELL
        /// AND TWO REFLECTIONS The orientation matrix is U in b
        /// need checks to see if the h,k,l vectors are input.
        /// </summary>
        /// <param name="r1"></param>
        /// <param name="r2"></param>
        /// <returns></returns>
        public int CalUBMatrix(Reflection r1, Reflection r2)
        {
            Vector u_phi1, u_phi2, h_c1, h_c2;
            
            // 1. Determine wavelength (Not sure how to put this to code)
            // mWaveLength = 1.0;
            //	wavelength=GetWaveLength();  need a GetWaveLength() function call into the class that creates this class.
            
            // 2. Generate h_phi from incident angles
            u_phi1 = new Vector(3);
            MotorIncidentAngles.VectorFromAngles(r1.motor_angles, mWaveLength, u_phi1);
            u_phi2 = new Vector(3);
            MotorIncidentAngles.VectorFromAngles(r2.motor_angles, mWaveLength, u_phi2);

#if _DBOUTPUT
            Console.WriteLine("\nCalculate UB Matrix");            
            Console.Write("u_phi1: with length {0}", u_phi1.Length());
            u_phi1.Print();
            Console.Write("u_phi2: with length {0}", u_phi2.Length());
            u_phi2.Print();
#endif

            // 3. Generete h_c
            //    i.e., calculate equation (23) in busing and levy.
            h_c1 = B.MultVector(r1.rlvector);
            h_c2 = B.MultVector(r2.rlvector);
#if _DBOUTPUT
            Console.WriteLine("\nCalculate UB Matrix");
            Console.Write("hc1: ");
            h_c1.Print();
            Console.Write("hc2: ");
            h_c2.Print();
#endif

#if false
            // 4. Not understand why this! 
            //is o in single...
            //calculate the angle between the vectors, and the dot pr.
            double dotpr = Vector.DotProduct(h1, h2);
            double temp1 = Vector.DotProduct(h1, h1);
            double temp2 = Vector.DotProduct(h2, h2);
            double v11 = DasMath.rad * Math.Acos(dotpr / Math.Sqrt(temp1 * temp2));

            //now calculate vl2 which is same as vl1 expect with the c vectors.
            dotpr = Vector.DotProduct(hc1, hc2);
            temp1 = Vector.DotProduct(hc1, hc1);
            temp2 = Vector.DotProduct(hc2, hc2);
            double vl2 = DasMath.rad * Math.Acos(dotpr / Math.Sqrt(temp1 * temp2));
#endif

            //5. calculate a matrix T_phi (used to name O)
            //O is matrix related to the motor angles. This is Tphi in Busing and Levy
            Matrix T_phi = new Matrix(3);
            Orthog(u_phi1, u_phi2, T_phi);
   
            // 6. Caculate matrix T_c (used to name C)
            //    C is matrix related to the crystal axis. This is Tc in Busing and Levy
            Matrix T_c = new Matrix(3);
            Orthog(h_c1, h_c2, T_c); 
            //equation 27 in busing and levy
#if _DBOUTPUT
            Console.WriteLine("Matrix T_phi (O):");
            T_phi.Print();
            Console.WriteLine("Matrix T_c (C)");
            T_c.Print();
#endif

            //7. now calcualte UB matrix and save....
            Matrix T_ctilt = new Matrix(3);
            double temp1;
            int errorcode = T_c.InvertMatrix(T_ctilt, out temp1);
            if (errorcode == 1)
            {
                //could not invert the matrix
                return -1;
            }
#if _DBOUTPUT
            Console.WriteLine("Matrix T_c^-1 (H):");
            T_ctilt.Print();            
#endif

            // 8. Calculate U & UB
            MulMatrix(T_phi, T_ctilt, U);
            MulMatrix(U, B, this);

#if _DBOUTPUT
            Console.WriteLine("Matrix U: ");
            U.Print();
            Console.WriteLine("Matrix UB: ");
            this.Print();
#endif

            // 9. Update calculation counter
            mCalUBCounter++;

            return 0;
        }
#endif




}
