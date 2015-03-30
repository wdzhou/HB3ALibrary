// #define _DBOUTPUT

using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Xml;

namespace DAS.UBMatrix
{
    /// <summary>
    /// UB matrix including all operations with UB matrix
    /// </summary>
    public class UBMatrix:Matrix
    {
        #region Static variable
        static bool _DBOUTPUT = false;

        #endregion



        #region private members
        /// <summary>
        /// wave length
        /// </summary>
        private double mWaveLength;
        private Matrix mU;
        private Matrix mB;
        private Matrix mG;
        private Matrix mGI;
        /// <summary>
        /// inversed UB 
        /// </summary>
        private Matrix m_invertUB;
        private Matrix matCovariance;
        private UnitCell mDirectCell;
        private UnitCell mRecipCell;
        private UnitCell mStoredDirectCell;
        // private int mInvUBCounter;
        private int mCalUBCounter;
        private double[] mTolerance;
        private double mChi2; 
        private Log mLog;

        #endregion

        #region Constructors
        /// <summary>
        /// Constructors
        /// </summary>
        public UBMatrix(double uiwavelength):this()
        {
            mWaveLength = uiwavelength;
            return;
        }

        /// <summary>
        /// Construct UB matrix without wavelength
        /// </summary>
        public UBMatrix()
            : base(3)
        {
            // 1. Matrix initialization
            mU = new Matrix(3);
            mB = new Matrix(3);
            m_invertUB = new Matrix(3);
            mG = new Matrix(3);
            mGI = new Matrix(3);
            matCovariance = new Matrix(3);
            mTolerance = new double[2];

            // 2. Wavelength and Unit cell
            mWaveLength = -1; // unphysical value
            for (int i = 0; i < 2; i++)
                mTolerance[i] = 1.0;
            mDirectCell = new UnitCell();
            mRecipCell = new UnitCell();
            mStoredDirectCell = new UnitCell();
            mChi2 = -1.0;

            PrevUB = new Matrix(3);
            mLog = new Log("c:\\ublog.txt");

            // 5. Set counter
            // mInvUBCounter = 0;
            mCalUBCounter = 0;

            return;
        }
        #endregion

        #region Public properties
        /// <summary>
        /// Wavelength
        /// </summary>
        public double WaveLength
        {
            get { return mWaveLength; }
            set { mWaveLength = value; }
        }
        /// <summary>
        /// UB matrix before refinement or etc. 
        /// </summary>
        public Matrix PrevUB;
        /// <summary>
        /// U matrix of this UB matrix
        /// </summary>
        public Matrix U
        {
            get { return mU; }
        }
        /// <summary>
        /// B matrix of this UB matrix
        /// </summary>
        public Matrix B
        {
            get { return mB; }
        }
        /// <summary>
        /// Covariance matrix
        /// </summary>
        public Matrix CovarianceMatrix
        {
            get { return matCovariance; }
        }
        /// <summary>
        /// Unit cell (direct space)
        /// </summary>
        public UnitCell CrystallCell
        {
            get { return mDirectCell; }
            set { mDirectCell = value; }
        }
        /// <summary>
        /// Tolerance for 2theta
        /// </summary>
        public double Tolerance_2Theta
        {
            set { SetTolerance(0, value); }
        }
        /// <summary>
        /// Tolerance for Eulerian angle
        /// </summary>
        public double Tolerance_Euler
        {
            set { SetTolerance(1, value); }
        }
        /// <summary>
        /// Chi^2 of least square fitting
        /// </summary>
        public double Chi2
        {
            get { return mChi2; }
        }
        #endregion

        #region private methods

        /// <summary>
        /// FIXME - should compare with Matrix.CalUB()... They are similar. 
        /// FIXME - Reorginaize the code so that it will reflect the correct solution
        /// </summary>
        /// <source>UnitCell.FOR -> MAKEUB </source>
        /// <param name="vH1"></param>
        /// <param name="vH2"></param>
        /// <param name="ang1"></param>
        /// <param name="ang2"></param>
        /// <param name="ierr"></param>
        private void _MakeUB(MillerIndices vH1, MillerIndices vH2, INST_ANGLES ang1, INST_ANGLES ang2, out int ierr)
        {
            string logmessage = String.Format("Generating UB Matrix (_MakeUB) From\n{0}: {1}\n{2}: {3}\n",
                vH1, ang1, vH2, ang2);

            ierr = 0;

#if true
            // This algorithm comes from Feng's (or say original 4-circle paper)
            // 1. Calcualte Tp
            Vector h1p = new Vector(3);
            Vector h2p = new Vector(3);
            Vector h3p = new Vector(3);
            INST_ANGLES.UnitVectorFromAngle(ang1, h1p);
            INST_ANGLES.UnitVectorFromAngle(ang2, h2p);
            Vector.CrossProduct(h1p, h2p, h3p);

            Vector t1p = new Vector(3);
            Vector t2p = new Vector(3);
            Vector t3p = new Vector(3);
            h1p.UnitVector(t1p);
            h3p.UnitVector(t3p);
            Vector.CrossProduct(t3p, t1p, t2p);

            Matrix mTp = new Matrix(3);
            for (int i = 0; i < 3; i++)
                mTp.M[i][0] = t1p.V[i];
            for (int i = 0; i < 3; i++)
                mTp.M[i][1] = t2p.V[i];
            for (int i = 0; i < 3; i++)
                mTp.M[i][2] = t3p.V[i];

            // 2. Calculate Tc'
            Vector vQc1 = new Vector(3);
            Vector vQc2 = new Vector(3);
            Vector vQc3 = new Vector(3);
            mB.MultVector(vH1, vQc1);
            mB.MultVector(vH2, vQc2);
            Vector.CrossProduct(vQc1, vQc2, vQc3);
            // Unify 
            Vector t1c = new Vector(3);
            Vector t2c = new Vector(3);
            Vector t3c = new Vector(3);
            vQc1.UnitVector(t1c);
            vQc3.UnitVector(t3c);
            Vector.CrossProduct(t3c, t1c, t2c);
            Matrix mTc_t = new Matrix(3);
            for (int i = 0; i < 3; i++)
                mTc_t.M[0][i] = t1c.V[i];
            for (int i = 0; i < 3; i++)
                mTc_t.M[1][i] = t2c.V[i];
            for (int i = 0; i < 3; i++)
                mTc_t.M[2][i] = t3c.V[i];

            Matrix.MulMatrix(mTp, mTc_t, mU);

#else
            This is previous way to calculate UB matrix from 2 measurements.  
            Algorithm was inherited from SINGLE code
            It is proved not that correct
            // 1. Get matrix for crystal Cartesian system, DC          
            Vector h1c = mB.MultVector(vH1);
            Vector h2c = mB.MultVector(vH2);
            logmessage += String.Format("B Matrix:\n{0}", mB);

            // 2. Make basis
            Matrix mDC = new Matrix(3);
            Matrix.MakeBasis(h1c, h2c, mDC);
            logmessage += String.Format("Basis Matrix (generated from Miller Indices):\n{0}", mDC);
#if _DBOUTPUT
            Console.WriteLine("h1c = {0}\nh2c = {1}", h1c.ToString(), h2c.ToString());
            Console.WriteLine("\t(_DB)  Matrix From h1c and h2c:");
            Console.WriteLine("{0}", mDC.ToString());
#endif

            // 3. Work on matrix in phi-coordinates
            Vector h1p = new Vector(3);
            Vector h2p = new Vector(3);
            INST_ANGLES.UnitVectorFromAngle(ang1, h1p);
            INST_ANGLES.UnitVectorFromAngle(ang2, h2p);

            Matrix mDP = new Matrix(3);
            Matrix.MakeBasis(h1p, h2p, mDP);
            logmessage += String.Format("Basis Matrix (generated from motor angles):\n{0}", mDP);

#if _DBOUTPUT
            Console.WriteLine("\t(_DB)  Matrix From h1phi and h2phi:");
            Console.WriteLine("{0}", mDP.ToString());
#endif

            // 3. Construct the U matrix, U*DC=DP
            Matrix invDC = new Matrix(3);
            double vol2;
            mDC.InvertMatrix(invDC, out vol2);

            // U = DP*DC^{-1} and Construct UB
            Matrix mU = new Matrix(3);
            Matrix.MulMatrix(mDP, invDC, mU);

            Console.WriteLine("Dp:");
            mDP.Print();
            Console.WriteLine("Dc-1:");
            invDC.Print();
            Console.WriteLine("U:");
            mU.Print();
            throw new Exception("Debug Stop");
#endif


            Matrix.MulMatrix(mU, mB, this);

            logmessage += String.Format("UB Matrix (generated)\n{0}", this);

            // 4. Write Log
            mLog.Record(logmessage);

            return;
        }

        /// <summary>
        /// Base step to calculate error of a UB matrix to an observed reflection
        /// </summary>
        /// <param name="iangle"></param>
        /// <param name="mindex"></param>
        /// <param name="diffangle"></param>
        /// <param name="cal2theta"></param>
        private void _CalculateError(INST_ANGLES iangle, MillerIndices mindex, out double diffangle, out double cal2theta)
        {
            MillerIndices calHKL = new MillerIndices();
            Vector caluPhi = new Vector(3);

            // 1. Caclulate (H, K, L) from incident angle
            CalculateHKLFromIncidentAngle(iangle, calHKL, caluPhi);

            // Step 2: Use calculated (HKL) to calculate 2theta
            double dtc;
            cal2theta = Calculate2ThetaFromHKL(calHKL, out dtc);

            // Step 3: Calculate difference and record
            Vector vC = this.MultVector(calHKL);
            double lvc = vC.Length();
            
            double dot = Vector.DotProduct(caluPhi, vC) / lvc;
            if (dot > 1.0)
                dot = 1.0;
            diffangle = DasMath.AcosD(dot);

            // Step 4: Outpout
            Console.Write("(H, K, L)_cal = ({0}, {1}, {2})\t", calHKL.H, calHKL.K, calHKL.L);
            Console.Write("2theta_obs = {0:0.000};  2theta_cal = {1:0.000} From UB matrix\t",
                iangle.twotheta, cal2theta);
            Console.WriteLine("Angle b/w u_phi_cal and u_phi_obs = {0:0.000}  (from {1:0.00})", 
                diffangle, dot);

            return;
        }

        /// <summary>
        /// Set tolerance
        /// </summary>
        /// <param name="index">int, 0 for 2theta, 1 for Eulerian angles</param>
        /// <param name="value">double, the value to set</param>
        private void SetTolerance(int index, double value)
        {
            if (value > 1.0E-10)
            {
                mTolerance[index] = value;
            }
            else
            {
                mTolerance[index] = 1.0;
            }
            return;
        }

        /// <summary>
        /// Calculate Cell parameters after the UB matrix is calculated or refined
        /// (Won't be necessary if UB matrix is generated from known cell parameters 
        /// and 2 reflections)
        /// </summary>
        private void _CalculateUnitCellFromUBMatrix()
        {
            Matrix mUBT = new Matrix(3);
            this.Transpose(mUBT);

            Matrix.MulMatrix(mUBT, this, mGI);
            mRecipCell.CalculateFromG(mGI);
            mRecipCell.CalculateReciprocalUnitCell(mDirectCell);
            mDirectCell.CalculateG(mG);

#if false
            Console.WriteLine("UB1146-1: Unit Cell {0}", mDirectCell);
            Console.WriteLine("UB1146-2: G^-1\n{0}", mGI);
#endif 

            return;

        }

        /// <summary>
        /// Calculate chi^2 of observed motor angles and miller indices 
        /// </summary>
        /// <param name="tReflections"></param>
        /// <returns></returns>
        private double _CalculateChi2(Reflection[] tReflections)
        {
            double diff2sum = 0;
            Vector uphiobs = new Vector(3);
            Vector uphihkl = new Vector(3);
            //Console.WriteLine("Length = {0}", tReflections.Length);

            for (int i = 0; i < tReflections.Length; i++)
            {
                INST_ANGLES.VectorFromAngles(tReflections[i].MotorAngles, mWaveLength, uphiobs);
                this.MultVector(tReflections[i].MillerIndex, uphihkl);
                double subsum = 0;
                for (int j = 0; j < 3; j++)
                {
                    subsum += (uphihkl.V[j] - uphiobs.V[j]) * (uphihkl.V[j] - uphiobs.V[j]);
                }
                diff2sum += subsum;
                //Console.WriteLine("{0}:  {1}  -->  {2}", i, subsum, diff2sum);
                //Console.WriteLine("\tObs:  {0}\n\tHKL:  {1}", uphiobs, uphihkl);
            }

            mChi2 = diff2sum;
            return diff2sum;
        }

        #endregion

        #region public methods

        /// <summary>
        /// Save UB matrix to file, including
        /// 1. direct cell
        /// 2. ub matrix
        /// </summary>
        /// <param name="pathfilename"></param>
        public void Save(string pathfilename)
        {
            // 1. Get time set
            DateTime currtime = DateTime.Now;

            XmlTextWriter myXMLTextWriter = new XmlTextWriter(pathfilename, System.Text.Encoding.UTF8);
            myXMLTextWriter.Formatting = Formatting.Indented;
            myXMLTextWriter.WriteStartDocument(false);

            myXMLTextWriter.WriteComment("UB Matrix Information");

            myXMLTextWriter.WriteStartElement("ubmatrix");

            // a. time                
            myXMLTextWriter.WriteStartElement("time");
            myXMLTextWriter.WriteAttributeString("time", currtime.ToString());
            myXMLTextWriter.WriteEndElement();
            // b. unit cell              
            myXMLTextWriter.WriteStartElement("unitcell");
            myXMLTextWriter.WriteAttributeString("a", mDirectCell.a.ToString());
            myXMLTextWriter.WriteAttributeString("b", mDirectCell.b.ToString());
            myXMLTextWriter.WriteAttributeString("c", mDirectCell.c.ToString());
            myXMLTextWriter.WriteAttributeString("alpha", mDirectCell.alpha.ToString());
            myXMLTextWriter.WriteAttributeString("beta", mDirectCell.beta.ToString());
            myXMLTextWriter.WriteAttributeString("gamma", mDirectCell.gamma.ToString());
            myXMLTextWriter.WriteEndElement();
            // c. wavelength
            myXMLTextWriter.WriteStartElement("wavelength");
            myXMLTextWriter.WriteAttributeString("lambda", mWaveLength.ToString());
            myXMLTextWriter.WriteEndElement();
            // d. ub matrix
            myXMLTextWriter.WriteStartElement("matrix");
            myXMLTextWriter.WriteAttributeString("matrix", XmlUtility.ConvertToString(this.M));
            myXMLTextWriter.WriteEndElement();

            myXMLTextWriter.WriteEndElement();

            myXMLTextWriter.Flush();
            myXMLTextWriter.Close();

            return;
        }


        /// <summary>
        /// Load from 
        /// </summary>
        /// <param name="pathfilename"></param>
        public void Load(string pathfilename)
        {
            List<string> keys = new List<string>();
            List<string> values = new List<string>();

            // 1. Parse from XML file
            XmlTextReader reader = new XmlTextReader(pathfilename);

            // a. To 'dictionary'
            while (reader.Read())
            {
                switch (reader.NodeType)
                {
                    case XmlNodeType.Element:
                        // Console.WriteLine(reader.Name);
                        while (reader.MoveToNextAttribute())
                        {
                            string attribute = reader.Name;
                            string value = reader.Value;
                            keys.Add(attribute);
                            values.Add(value);
                            // Console.WriteLine("  > {0} = {1}", attribute, value);
                        }
                        break;
                    default:
                        break;
                }
            }

            // b. Construct values
            UnitCell uc = new UnitCell();
            uc.a = Convert.ToDouble(values[keys.IndexOf("a")]);
            uc.b = Convert.ToDouble(values[keys.IndexOf("b")]);
            uc.c = Convert.ToDouble(values[keys.IndexOf("c")]);
            uc.alpha = Convert.ToDouble(values[keys.IndexOf("alpha")]);
            uc.beta = Convert.ToDouble(values[keys.IndexOf("beta")]);
            uc.gamma = Convert.ToDouble(values[keys.IndexOf("gamma")]);
            double lambda = Convert.ToDouble(values[keys.IndexOf("lambda")]);

            double[][] ubmatrix = new double[3][];
            for (int i = 0; i < 3; i++)
                ubmatrix[i] = new double[3];
            XmlUtility.ConvertTo2DDoubleArray(values[keys.IndexOf("matrix")], ubmatrix);

#if false
            Console.WriteLine("Loading UnitCell:\n\t> {0}", uc);
            Console.WriteLine("\t> Lambda = {0}\n > Matrix: \n", lambda);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    Console.Write("{0:0.000000}  ", ubmatrix[i][j]);
                Console.WriteLine();
            }
#endif

            // 2. Load unitcell, wavelength, and UB matrix, Inverted UB matrix to it
            mWaveLength = lambda;
            SetDirectUnitCell(uc);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    this.M[i][j] = ubmatrix[i][j];
            double det;
            InvertMatrix(m_invertUB, out det);

            string logmessage = String.Format("Loading UB Matrix From File {0}", pathfilename);
            logmessage += String.Format("UB Matrix Loaded:\n{0}", this);
            logmessage += String.Format("Inverted UB Matrix:\n{0}", m_invertUB);
            mLog.Record(logmessage);

            return;
        }

        /// <summary>
        /// Set the direct unit cell for this UB matrix
        /// This uses the method of CBF to determine the reciprocal lattice parameters 
        /// without a matrix inversion.
        /// </summary>
        /// <param name="inpcell"></param>
        public void SetDirectUnitCell(UnitCell inpcell)
        {
            // 0. Check prerequisit
            if (mWaveLength <= 0)
            {
                string expmsg = String.Format("Wavelength = {0} Is Not Physical", mWaveLength);
                throw new Exception(expmsg);
            }

            // 1. Set up unit cell
            inpcell.CopyTo(mDirectCell);

            // 2. Calculate G matrix
            mDirectCell.CalculateG(mG);

            // 3. Set up reciprocal lattice
            mDirectCell.CalculateReciprocalUnitCell(mRecipCell);
            mRecipCell.CalculateG(mGI);

#if _DEBUG
            Console.WriteLine("Debug Output 1017 (UBMatrix)");
            Console.WriteLine("Direct Cell:{0}", mDirectCell);
            Console.WriteLine("G Matrix:\n{0}", mG);
            Console.WriteLine("Reciprocal Lattice: {0}", mRecipCell);
            Console.WriteLine("G^-1 matrix:\n{0}", mGI);

            Matrix testG = new Matrix(3);
            UnitCell testRecipCell = new UnitCell();
            UnitCell testDirecCell = new UnitCell();

            testRecipCell.CalculateFromG(mGI);
            testRecipCell.CalculateReciprocalUnitCell(testDirecCell);
            testDirecCell.CalculateG(testG);

            Console.WriteLine("Reciprocal Lattice <- G^-1: {0}", testRecipCell);
            Console.WriteLine("Realspace Lattice  <- Reciprocal lattice: {0}", testDirecCell);
            Console.WriteLine("G <- Real Space Lattice:\n{0}", testG);

            throw new Exception("Debug Stop"); 
#endif

            // 4. Generate B matrix
            CalBMatrix();

            // 5. Log
            mLog.Record(String.Format("Set Wave-length = {0}", mWaveLength.ToString()));
            mLog.Record(String.Format("Set Unit Cell:\n{0}", mDirectCell));
            mLog.Record(String.Format("Set Reciprocal Cell:\n{0}", mRecipCell));
            mLog.Record(String.Format("Set G Matrix: \n{0}", mG));
            mLog.Record(String.Format("set G^-1 Matrix:\n{0}", mGI));
            mLog.Record(String.Format("Set B matrix:\n{0}", mB));

            return;
        }

        /// <summary>
        /// Get direct cell information
        /// </summary>
        /// <param name="edges"></param>
        /// <param name="angles"></param>
        public void GetDirectUnitCell(double[] edges, double[] angles)
        {
            string logmessage = String.Format("UB Matrix:  Unit Cell = {0}", mDirectCell);
            mLog.Record(logmessage);

            for (int i = 0; i < 3; i++)
                edges[i] = mDirectCell.Axiles[i];
            for (int i = 0; i < 3; i++)
                angles[i] = mDirectCell.Angles[i];
            
            return;
        }

        /// <summary>
        /// Sotre the current Unit Cell
        /// </summary>
        public void StoreUnitCell()
        {
            mDirectCell.CopyTo(mStoredDirectCell);
            return;
        }

        /// <summary>
        /// Restore the Unit Cell and related objects
        /// </summary>
        public void RestoreUnitCell()
        {
            SetDirectUnitCell(mStoredDirectCell);
            return;
        }

        /// <summary>
        /// Calcualte B matrix from the metrical matrix of reciprocal unit cell
        /// <source>UnitCell.For -> MAKEB(...)</source>
        /// B is caluclated directly has describe in Equation-3 of bushing and levy
        /// </summary>
        public void CalBMatrix()
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    mB.M[i][j] = 0.0;

            mB.M[0][0] =  mRecipCell.a;
            mB.M[0][1] =  mRecipCell.b * DasMath.CosD(mRecipCell.gamma);
            mB.M[1][1] =  mRecipCell.b * DasMath.SinD(mRecipCell.gamma);
            mB.M[0][2] =  mRecipCell.c * DasMath.CosD(mRecipCell.beta);
            mB.M[1][2] = -mRecipCell.c * DasMath.SinD(mRecipCell.beta) * DasMath.CosD(mDirectCell.alpha);
            mB.M[2][2] = 1 / mDirectCell.c;

            return;
        }

        /// <summary>
        /// Make UB matrix by 2 reflections (motor angles and miller index)
        /// This is a public wrapping
        /// </summary>
        /// <param name="vH1"></param>
        /// <param name="vH2"></param>
        /// <param name="ang1"></param>
        /// <param name="ang2"></param>
        /// <param name="ierr"></param>
        /// <returns></returns>
        public void MakeUB(MillerIndices vH1, MillerIndices vH2, INST_ANGLES ang1, INST_ANGLES ang2, out int ierr)
        {
            // 1. Calculate UB matrix from 2 reflections
            _MakeUB(vH1, vH2, ang1, ang2, out ierr);
            double det;
            InvertMatrix(m_invertUB, out det);

            // 2. Logging
            mLog.Record(String.Format("UB Matrix Calculated:\n{0}\nFrom 1: {1}  {2}\n     2: {3}  {4}\n", 
                this, vH1, ang1, vH2, ang2));

            return;
        }

        /// <summary>
        /// CALCULATE ANGLES FROM VECTOR - IERR NON-ZERO IF NOT ACCESSIBLE            
        /// This is the algorithm to calculate the incident angle from given (HKL) and UB matrix
        /// </summary>
        /// <param name="inphkl">Input Miller index</param>
        /// <param name="angl">Output motor angles calculated from HKL and UB matrix</param>
        /// <param name="fixed_phi"></param>
        /// <param name="ierror"></param>
        public void CalMotorAnglesFromMillerIndex(MillerIndices inphkl, INST_ANGLES angl, bool fixed_phi, out int ierror)
        {
            // 1. Calculate h_phi
            Vector thphi = new Vector(3);
            this.MultVector(inphkl, thphi);

            // 2. calculate sin theta
            double dinv = thphi.Length();
            double dinvhkl = inphkl.Length();
            double stheta = dinv * mWaveLength * 0.5;
            // Console.WriteLine("Debug ubmatrix-107: d^-1 = {0}, sin(theta) = {1}", dinv, stheta);

            if (dinv <= 1.0e-6 || stheta > 1.0)
            {
                // Input vector is not good
                angl.chi = 0.0;
                angl.omega = 0.0;
                angl.phi = 0.0;
                angl.twotheta = 0.0;

                ierror = -1;
            }
            else
            {
                // 1. calculate two theta
                angl.twotheta = 2.0 * DasMath.AsinD(stheta);

                // 2. calculate peak position
                if (fixed_phi == false)
                {
                    // 2.1 in bisecting mode as (38)
                    angl.omega = 0.0;
                    angl.chi = Math.Atan2(thphi.k, thphi.Project(0, 1)) * DasMath.radin;
                    angl.phi = Math.Atan2(thphi.j, thphi.i) * DasMath.radin;
                }
                else
                {
                    // 2.2 NOT bisecting mode - calculate in fixed phi mode
                    angl.omega = Math.Asin(thphi.j / dinv) * DasMath.rad;
                    angl.chi = Math.Atan2(thphi.k, thphi.i) * DasMath.rad;
                    angl.phi = 0.0;
                }

                ierror = 0;
            }

            mLog.Record(String.Format("Calculate Motor Angle:   {0} ==> {1}\n",
                inphkl, angl));

            return;
        }

        /// <summary>
        /// ROUTINE TO CALCULATE AN ORIENTATION MATRIX FROM THREE reflections
        /// need checks to see if the h,k,l vectors are input.
        /// </summary>
        /// <param name="r1">Reflection</param>
        /// <param name="r2">Reflection</param>
        /// <param name="r3">Reflection</param>
        /// <returns></returns>
        public int CalUBMatrix(Reflection r1, Reflection r2, Reflection r3)
        {
            Matrix Hphi, Hi, HiInv;  //h is a temp matrix to hold 
            Vector[] vHs = new Vector[3]; // h1, h2, h3
            Reflection[] refls = new Reflection[3];

            refls[0] = r1;
            refls[1] = r2;
            refls[2] = r3;
  
            // FIXME - It is not quite right with the equation... Equation (28) is missed!!!!!!!!!!!!!!!!!
            // need a GetWaveLength() function call into the class that creates this class.
            for (int i = 0; i < 3; i ++)
                INST_ANGLES.VectorFromAngles(refls[i].MotorAngles, mWaveLength, vHs[i]);

            //form matrices...
            Hphi = new Matrix(3);
            for (int i = 0; i < 3; i++)
                Hphi.SetColumn(i, vHs[i]);

            Hi = new Matrix(3);
            for (int i = 0; i < 3; i++)
                Hi.SetColumn(i, refls[i].MillerIndex);

            double temp1;
            HiInv = new Matrix(3);
            int errorcode = Hi.InvertMatrix(HiInv, out temp1);
            if (errorcode == 1)
            {
                //could not invert the matrix
                return -1;
            }

            MulMatrix(Hphi, HiInv, this);

            // -1 Update counter
            mCalUBCounter++;

            // -1 Calculate unit cell and B matrix
            _CalculateUnitCellFromUBMatrix();

            return 0;
        }

        /// <summary>
        /// h_{hkl} =  u x (UB)^-1 /( lambda / (2 sin theta ) ) = u x (UB)^-1 x q
        /// Class variable mInvUB will be used, other than UB matrix
        /// </summary>
        /// <param name="iangle">(In) incident angle</param>
        /// <param name="hkl">(Out) miller indices</param>
        /// <param name="uPhi">(Out) u_phi</param>
        public void CalculateHKLFromIncidentAngle(INST_ANGLES iangle, MillerIndices hkl, Vector uPhi)
        {
            INST_ANGLES.UnitVectorFromAngle(iangle, uPhi);
            double dt = mWaveLength / (2.0 * DasMath.SinD(0.5 * iangle.twotheta));
            Vector vHH = m_invertUB.MultVector(uPhi);

            // Console.WriteLine("[CalHKL] UB Matrix: \n{0}", this);
            // Console.WriteLine("[CalHKL] UB^-1    : \n{0}", m_invertUB);

            // Console.WriteLine("Input motor angle: {0}. ", iangle);

            if (false)
            {
                Console.WriteLine("U_phi:");
                uPhi.Print();
                Console.WriteLine("DB 425A  vHH:");
                vHH.Print();
                Matrix mTest = new Matrix(3);
                Matrix.MulMatrix(m_invertUB, this, mTest);
                Console.WriteLine("This is supposed to be unit");
                mTest.Print();
            }
            
            // Change to nearest integer
            for (int i = 0; i < 3; i++)
            {
                vHH.V[i] = vHH.V[i] / dt;
            }

            // hkl.V[i] = DasMath.NearestInteger(vHH.V[i] / dt);
            hkl.SetFromCalculatedValue(vHH);

            return;
        }

        /// <summary>
        /// Calculate 2theta from given (HKL) and UB matrix
        /// </summary>
        /// <param name="vhkl"></param>
        /// <param name="d">1/q (in fact) for inverse of d-spacing</param>
        /// <returns></returns>
        public double Calculate2ThetaFromHKL(MillerIndices vhkl, out double d)
        {
            d = Vector.DotProduct(vhkl, mGI.MultVector(vhkl));
            double d2 = mWaveLength * Math.Sqrt(d) / 2.0;
            double ttc = 2.0 * DasMath.AsinD(d2);
            return ttc;
        }

        /// <summary>
        /// Calculate the angle between 2 miller indices
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        public double CalculateAngle2MillerIndicies(MillerIndices m1, MillerIndices m2)
        {
            // 1. Normalize
            MillerIndices n1 = new MillerIndices();
            MillerIndices n2 = new MillerIndices();
            m1.Normalize(n1, mGI);
            m2.Normalize(n2, mGI);

            // 2. Dot product m1 x G^-1  x m2
            Vector temp1 = new Vector(3);
            mGI.MultVector(n2, temp1);
            double dotprod = Vector.DotProduct(n1, temp1);

            // 3. Acos
            double ang12 = DasMath.AcosD(dotprod);

            return ang12;
        }

        /// <summary>
        /// Generate a list of reflections in the order from largest (H, K, L)
        /// </summary>
        /// <param name="limits"></param>
        /// <returns></returns>
        public List<MillerIndices> GenerateReflectionListInverse(int[] limits)
        {
            List<MillerIndices> hkllist = new List<MillerIndices>();
            for (int h = limits[0]; h >= -limits[0]; h--)
            {
                for (int k = limits[1]; k >= -limits[1]; k--)
                {
                    for (int l = limits[2]; l >= -limits[2]; l--)
                    {
                        MillerIndices ml = new MillerIndices();
                        ml.H = h;
                        ml.K = k;
                        ml.L = l;
                        hkllist.Add(ml);
                    }
                }
            }
            return hkllist;
        }

        /// <summary>
        /// Check Error of a UB matrix for
        /// a. calculated 2theta vs. observed 2theta
        /// b. |uphi_obs dot uphi_cal| == 1???
        /// </summary>
        /// <param name="iangle"></param>
        /// <param name="mindex"></param>
        /// <param name="toleranceiangle"></param>
        /// <param name="tolerancemindex"></param>
        /// <param name="diff2theta"></param>
        /// <param name="diffangle"></param>
        /// <returns></returns>
        public bool CheckError(INST_ANGLES iangle, MillerIndices mindex, double toleranceiangle, double tolerancemindex,
            out double diff2theta, out double diffangle)
        {
            bool returnbool = true;

            // 1. Calculate 2theta and angle of u_phi_cal dot u_phi_obs
            double ttc, angdot;
            _CalculateError(iangle, mindex, out angdot, out ttc);

            // 2.  Compare 2theta
            diff2theta = Math.Abs(ttc - iangle.twotheta);
            diffangle = angdot;

            if (diff2theta > tolerancemindex)
                returnbool = false;
            if (angdot > toleranceiangle)
                returnbool = false;

            return returnbool;
        }

        /// <summary>
        /// Routine to test how many reflections in list can be fitted by matrix UB
        /// </summary>
        /// <source>UnitCell.For: TestRef()</source>
        /// <param name="IncidentAngles"></param>
        /// <param name="calHKLs"></param>
        /// <param name="diff2theta"></param>
        /// <param name="diffuphi"></param>
        /// <param name="nfound"></param>
        /// <returns></returns>
        public bool TestUBByReflections(INST_ANGLES[] IncidentAngles, MillerIndices[] calHKLs, 
            double[] diff2theta, double[] diffuphi, out int nfound)
        {
            nfound = 0;
            int ncount = calHKLs.Length;

            for (int icount = 0; icount < ncount; icount++)
            {
                // Step 1: Use UB^-1 to computer (HKL) from Incident angle
                MillerIndices vhkl = new MillerIndices();
                Vector uPhi = new Vector(3); // u_phi
                CalculateHKLFromIncidentAngle(IncidentAngles[icount], vhkl, uPhi);

                // Step 2: Use calculated (HKL) to calculate 2theta
                double dtc;
                double ttc = Calculate2ThetaFromHKL(vhkl, out dtc);

                // Step 3: Calculate difference and record
                Vector vC = this.MultVector(vhkl);
                double lvc = vC.Length();
                double dot = Vector.DotProduct(uPhi, vC) / lvc;
                if (dot > 1.0)
                    dot = 1.0;
                double angdot = DasMath.AcosD(dot);

                // 4. Record
                calHKLs[icount] = vhkl;
                diff2theta[icount] = Math.Abs(ttc - IncidentAngles[icount].twotheta);
                diffuphi[icount] = angdot;

                // Test if both two-theta and Eulerian angle are within tolerance
                if ( diff2theta[icount] <= mTolerance[0] )
                {
                    if (angdot < mTolerance[1])
                    {
                        // The Eulerian angle is within tolerance 
                        nfound = nfound + 1;
                    }
                }
            } // icount

            bool rfound = false;
            if (nfound > 0)
            {
                rfound = true;
            }

            return rfound;
        } // TestRef()

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
        public bool GenerateUBFromLatticeRefine(INST_ANGLES[] IncidentAngles)
        {
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
                INST_ANGLES.UnitVectorFromAngle(IncidentAngles[n1], vH1);
                vH1.UnitVector(uH1);

                if (false)
                {
                    Console.WriteLine("(_DB)   Incident Angle 1 and u_phi1:");
                    Console.WriteLine("{0} | {1}\n", IncidentAngles[n1].ToString(), uH1.ToString());

                    mLog.Record(String.Format("Incident Angle 1", IncidentAngles[n1].ToString()));
                    mLog.Record(String.Format("u_phi 1", uH1.ToString()));
                }


                for (int n2 = n1 + 1; n2 < ncount; n2++)
                {
                    // old: ICalcH(IncidentAngles[n2], vH2);
                    INST_ANGLES.UnitVectorFromAngle(IncidentAngles[n2], vH2);
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

                        mLog.Record(String.Format("Incident Angle 2", IncidentAngles[n2].ToString()));
                        mLog.Record(String.Format("u_phi 2", uH2.ToString()));
                        mLog.Record(String.Format("Angle between u_phi 1 and u_phi 2", Convert.ToString(ang12)));
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
#if _DBOUTPUT
                    Console.WriteLine("\t(_DB)  Limit of HKL = ({0}, {1}, {2})", Lim[0], Lim[1], Lim[2]);
                    Logging("Limit of (HKL) of n1 and n2",
                        String.Format("{0}, {1}, {2})", Lim[0], Lim[1], Lim[2]));
#endif

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
                            mTolerance[0])
                        {

                            // ii. FIND HKL FOR SECOND REFLECTION 
                            for (int jhkl = 0; jhkl < sizehklist; jhkl++)
                            {
                                MillerIndices hkl2 = hkllist[jhkl];
                                double d2;
                                double twotheta_cal_2 = Calculate2ThetaFromHKL(hkl2, out d2);

                                if (Math.Abs(IncidentAngles[n2].twotheta - twotheta_cal_2) <
                                    mTolerance[0])
                                {
                                    // Match in 2theta for 2nd reflection. **      

                                    // iii. Check: Is the angle between the relections correct?
                                    double da1 = Vector.DotProduct(hkl1, mGI.MultVector(hkl2));
                                    double da2 = DasMath.AcosD(da1 / Math.Sqrt(d1 * d2));
                                    if (Math.Abs(da2 - ang12) < mTolerance[1])
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
                                            if (diff2ThetaList[ifc] > mTolerance[0])
                                                overlimit++;

                                        if (nfound >= 2 && overlimit < nfound-1)                                        
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
#if _DBOUTPUT
                Console.WriteLine("Best UB Can Fit {0} Reflections", nfit);
                Logging("Number of Reflections Fit By Best UB", Convert.ToString(nfit));
#endif
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

            return true;
        } // GenUB 

        /// <summary>
        /// VZ's version to generate UB matrix from a series of reflections and 
        /// known lattice parameters
        /// The best fitting UB matrix is stored in mM[][], while
        /// the next 10 best fitting UB matrices are stored in a list too. 
        /// </summary>
        /// <param name="IncidentAngles"></param>
        /// <param name="numreflections">int, number of reflections (in case input list is larger)</param>
        /// <returns></returns>
        public bool GenerateUBFromUnIndexedReflectiuons(INST_ANGLES[] IncidentAngles, int numreflections)
        {

            // 1. Initialization
            int ncount;
            if (numreflections >= 0)
                ncount = numreflections;
            else
                ncount = IncidentAngles.Length;
            Matrix mUBBest = new Matrix(3);
            Matrix mInvUBBest = new Matrix(3);
            int[] foundcount = new int[ncount + 1];
            List<int>[] mFoundList = new List<int>[ncount + 1];
            List<double>[] mErrorList = new List<double>[ncount + 1];
            for (int i = 0; i < foundcount.Length; i++)
            {
                mFoundList[i] = new List<int>();
                mErrorList[i] = new List<double>();
            }

            // 2. Find all possible combinations of 2 reflections                        
            Vector vH1 = new Vector(3);
            Vector vH2 = new Vector(3);


            MillerIndices[] calHKLList = new MillerIndices[ncount];
            double[] diff2ThetaList = new double[ncount];
            double[] diffUPhiList = new double[ncount];

            // 2.2 Find the angle pairs
            int c12 = ncount * (ncount - 1) / 2;
            double[] ang12list = new double[c12];
            int[] pair1list = new int[c12];
            int[] pair2list = new int[c12];
            ListPairAngles(IncidentAngles, ncount, ang12list, pair1list, pair2list);

            // 3. Calculate HKL limit list
            int[][] Limits = new int[ncount][];
            for (int n = 0; n < ncount; n++)
                Limits[n] = new int[3];
            for (int n1 = 0; n1 < ncount; n1++)
            {
                double tmax = Math.Abs(IncidentAngles[n1].twotheta) / 2.0 + 1.0;
                for (int i = 0; i < 3; i++)
                    Limits[n1][i] = Convert.ToInt32((2.0 * DasMath.SinD(tmax)) /
                        (mWaveLength * Math.Sqrt(mGI.M[i][i])) + 1.0);
            }

            // 4. Brute force loop around all possible combination of HKL
            // Note that sorting is from smallest value to the largest
            bool continueserach = true;
            int pairindex = ang12list.Length - 1;
            while (continueserach == true)
            {
                // a) Generate reflection list
                int n1 = pair1list[pairindex];
                int n2 = pair2list[pairindex];
                int[] Lim;
                double max2theta;
                if (IncidentAngles[n1].twotheta > IncidentAngles[n2].twotheta)
                {
                    Lim = Limits[n1];
                    max2theta = n1;
                }
                else
                {
                    Lim = Limits[n2];
                    max2theta = n2;
                }
                if (max2theta < 1.0E-16)
                {
                    string msg = "Error!  Maximum 2Theta Of Two Measured Reflections is ZERO.  Input Reflections Error! Job aborted.";
                    Console.WriteLine(msg);
                    throw new Exception(msg);
                }

                List<MillerIndices> hkllist = GenerateReflectionListInverse(Lim);
                int sizehklist = hkllist.Count;

                Console.WriteLine("UB223-1 Reflection Pair \t{0}: {1}\n\t\t\t\t{2}: {3}", n1, IncidentAngles[n1],
                    n2, IncidentAngles[n2]);
                Console.WriteLine("UB223-2 \tAngle = {8}  (HKL) Limit for {0} = {1}, {2}, {3}\t\t(HKL) Limit for {4} = {5}, {6}, {7}",
                    n1, Limits[n1][0], Limits[n1][1], Limits[n1][2], n2, Limits[n2][0], Limits[n2][1], Limits[n2][2], ang12list[pairindex]);

                Console.WriteLine("UB225:  Press Any Key To Continue");
                Console.ReadLine();

                // b)loop around possible HKL for n1 and n2
                int maxfound = 0;
                for (int i = 0; i < foundcount.Length; i++)
                {
                    foundcount[i] = 0;
                    mFoundList[i].Clear();
                    mErrorList[i].Clear();
                }

                for (int ihkl = 0; ihkl < sizehklist; ihkl++)
                {
                    // n1: Guess HKL for First reflection
                    MillerIndices hkl1 = hkllist[ihkl];
                    double d1;
                    double twotheta_cal_1 = Calculate2ThetaFromHKL(hkl1, out d1); //Calculate2Theta(hkl1, mGI, wavelength, out d1);

                    if (Math.Abs(IncidentAngles[n1].twotheta - twotheta_cal_1) <
                        mTolerance[0])
                    {

                        // n2:  Guess HKL for Second reflection
                        for (int jhkl = 0; jhkl < sizehklist; jhkl++)
                        {
                            MillerIndices hkl2 = hkllist[jhkl];
                            double d2;
                            double twotheta_cal_2 = Calculate2ThetaFromHKL(hkl2, out d2);

                            if (Math.Abs(IncidentAngles[n2].twotheta - twotheta_cal_2) <
                                mTolerance[0])
                            {
                                // d) 3rd Check: Is the angle between the relections correct?
                                double da1 = Vector.DotProduct(hkl1, mGI.MultVector(hkl2));
                                double da2 = DasMath.AcosD(da1 / Math.Sqrt(d1 * d2));

                                if (Math.Abs(da2 - ang12list[pairindex]) < mTolerance[1])
                                {
                                    // e) All loose requirements meet, then caluclate UB matrix, UB^-1                                
                                    int ierr;
                                    _MakeUB(hkl1, hkl2, IncidentAngles[n1], IncidentAngles[n2], out ierr);
                                    if (ierr == 1)
                                        throw new Exception("Error In Generating UB Matrix!  Quit Program!");
                                    double det;
                                    this.InvertMatrix(m_invertUB, out det);

                                    Console.WriteLine("D555-1:  n1 = Reflection {0} Try: {1} \t\t n2 = Reflection {2}  Try: {3}   ({0}, {2}) = {4}",
                                        n1, hkl1, n2, hkl2, ang12list[pairindex]);
                                    this.Print();

                                    // f) Test with other reflections                                    
                                    int nfound;
                                    TestUBByReflections(IncidentAngles, calHKLList, diff2ThetaList, diffUPhiList, out nfound);
                                    double sum2theta = 0.0;
                                    double sumangles = 0.0;
                                    for (int r = 0; r < ncount; r++)
                                    {
                                        sum2theta += diff2ThetaList[r] * diff2ThetaList[r];
                                        sumangles += diffUPhiList[r] * diffUPhiList[r];
                                    }
                                    double sigma2theta = Math.Sqrt(sum2theta/ncount);
                                    double sigmaangles = Math.Sqrt(sumangles/ncount);                                   

                                    if (nfound > maxfound)
                                        maxfound = nfound;
                                    foundcount[nfound] += 1;
                                    mFoundList[nfound].Add(ihkl);
                                    mFoundList[nfound].Add(jhkl);
                                    mErrorList[nfound].Add(sigma2theta);
                                    mErrorList[nfound].Add(sigmaangles);

                                    Console.WriteLine("D555-2:  NFound = {0}\n", nfound);
                                    // Console.ReadLine();
                                }

                            }

                        } // j-hkl
                    } // if hkl-1 meets tolerance
                } // i-hkl

                // Compare the result
                for (int i = mFoundList.Length - 1; i >= 0; i--)
                {
                    if (mFoundList[i].Count > 0)
                    {
                        Console.WriteLine("UB348:  {0} ... {1}/{2}", i, mFoundList[i].Count/2, mErrorList[i].Count/2);
                        for (int j = 0; j < mFoundList[i].Count / 2; j++)
                        {
                            int t1 = mFoundList[i][2 * j];
                            int t2 = mFoundList[i][2 * j + 1];
                            Console.WriteLine("{0}:  {1} = {3}     {2} = {4}", 
                                j, hkllist[t1], hkllist[t2], mErrorList[i][2*j], mErrorList[i][2*j+1]);
                        }
                    }

                }

                Console.WriteLine("Debug 457   Max Found = {0}  For ({1}, {2})= {3} of {4} Loop",
                    maxfound, n1, n2, ang12list[pairindex], pairindex);
                for (int i = 0; i < foundcount.Length; i++)
                {
                    Console.WriteLine("Number of Matrix To Fit {0} Input Reflections = {1}",
                        i, foundcount[i]);
                }

                Console.WriteLine("Pre Exit For Testing 152");
                Environment.Exit(0);

                Console.ReadLine();
                if (pair1list.Length - pairindex > 5)
                    break;

                // -1: Loop control
                if (pairindex == 0)
                    continueserach = false;
                else
                {
                    pairindex--;
                    if (ang12list[pairindex] < 20.0)
                        continueserach = false;
                }

            }

            return true;

        } // End of Function   


        /// <summary>
        /// Calcualte UB matrix error as
        /// \sum_{reflections} ([HKL] - [hkl])^2
        /// as (hkl) is the nearest integer version of (HKL)
        /// </summary>
        /// <param name="reflecslist"></param>
        /// <returns></returns>
        public double CalculateUBMatrixError(Reflection[] reflecslist)
        {
            // Check
            if (reflecslist.Length == 0)
                throw new Exception("User input empty reflection list for error.  Incorrect.");

            double error = 0;
            Vector cv = new Vector(3);                            
            MillerIndices temphkl = new MillerIndices();   
            MillerIndices temphklint = new MillerIndices();
            for (int i = 0; i < reflecslist.Length; ++i)
            {
                INST_ANGLES tempmotorangle = reflecslist[i].MotorAngles;

                // Calculate (hkl) from UB matrix
                CalculateHKLFromIncidentAngle(tempmotorangle, temphkl, cv);
                // Console.WriteLine("Calculated reuslt: {0}", temphkl);

                // Calcualte difference
                double temperror = 0;
                for (int j = 0; j < 3; ++j)
                {
                    double dev = temphkl.Error;
                    temperror += dev * dev;
                }

                error += temperror;
            }

            error = Math.Sqrt(error)/reflecslist.Length;

            return error;
        }

        /// <summary>
        /// Refine UB matrix withh Linear least square fitting
        /// Code is cleaned from SINGLE (used be versioned by VZ)
        /// a) linear least square fitting
        /// b) reserve covariant matrix
        /// c) calculate chi^2
        /// d) calculate refined lattice parameter
        /// </summary>
        /// <param name="reflectsList"></param>
        /// <param name="nlsq"></param>
        /// <returns></returns>
        public bool RefineUBLinearLSF(Reflection[] reflectsList, int nlsq)
        {
            string message = "";

            //************************************************************************************
            // Check Input's Validity & Logging
            //************************************************************************************
            // Check # of input reflections
            if (nlsq < 3)
            {
                string errmsg = "At least 3 Reflections Are Needed For Least Square Fitting";
                Console.WriteLine(errmsg);
                mLog.Record(String.Format("Refining UB Matrix Error Message: {0}", errmsg));
                return false;
            }
            else if (nlsq < reflectsList.Length)
            {
                string errmsg = String.Format("Input number of reflections ({0}) is smaller than 'nlsq' ({1}), {2}.",
                    reflectsList.Length, nlsq, "which will cause index error");
                Console.WriteLine(errmsg);
                mLog.Record(String.Format("Refining UB Matrix Error Message: {0}", errmsg));
                return false;
            }
            else if (nlsq > reflectsList.Length)
            {
                Console.WriteLine("Refining UB matrix:  using first {0} reflections out of {1} input reflections.",
                    nlsq, reflectsList.Length);
            }
              

            // Check consistent of input (H, K, L), motor angles and current UB matrix
            message += "Check Consistency Among Input (H, K, L), Incident Angles and UB Matrix\n";
            MillerIndices cmi = new MillerIndices();
            Vector cv = new Vector(3);
            for (int i = 0; i < nlsq; i++)
            {
                // Verify whether input HKL values are same as calculated HKL values from present UB
                MillerIndices inputhkl = reflectsList[i].MillerIndex;
                // Console.Write("Reflection {0}: user input: {1}; ", i, inputhkl);
                
                CalculateHKLFromIncidentAngle(reflectsList[i].MotorAngles, cmi, cv);
                // Console.Write("calculated value: {0}, ", cmi);

                cmi.ConvertToIntegers();
                // Console.WriteLine("converted to integer: {0}", cmi); 

                // Use the calculated version instead
                if (!inputhkl.EqualIntegerValue(cmi))
                {
                    message += String.Format("Input Reflection {0} With Invalid (HKL) = {1}",
                        i, reflectsList[i].MillerIndex);
                    cmi.CopyTo(reflectsList[i].MillerIndex);
                    message += String.Format(" Replaced By {0}\n", reflectsList[i].MillerIndex);
                }
            }

            if (_DBOUTPUT)
            {
                message += "Reflection For Least Square Fitting UB Matrix:\t\t(DB335)\n";
                for (int n = 0; n < nlsq; n++)
                    message += String.Format("{0}\t{1}\n", n, reflectsList[n]);
                message += String.Format("UB Matrix Before Refinement:\n{0}", this);
                Console.WriteLine(message);
            }

            //************************************************************************************
            // Refine UB matrix linearly
            //************************************************************************************
            // Zero Arrays
            Matrix mYT, mH, mY;
            mYT = new Matrix(3);
            mH = new Matrix(3);
            mY = new Matrix(3);

            // Transfer angles and indices from list to temp arrays for Least square fitting
            // Original code:  tempa (angle)  tempi (reflection hkl)
            Vector[] hphiList = new Vector[nlsq];
            for (int i = 0; i < nlsq; i++)
            {
                // 2.1 Incident angle and h_phi
                hphiList[i] = new Vector(3);
                INST_ANGLES.VectorFromAngles(reflectsList[i].MotorAngles, mWaveLength, hphiList[i]);
            }

            // 3. form matrices H(phi) and H of rank 3 x N and do matrix products
            for (int n = 0; n < nlsq; n++)
            {
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        mYT.M[i][j] += hphiList[n].V[j] * reflectsList[n].MillerIndex.V[i];
                        mH.M[i][j] += reflectsList[n].MillerIndex.V[i] * reflectsList[n].MillerIndex.V[j];
                    }
                }
            }

            // 4. Inverse matrix H7, Transpose Matrix H6 and use euqation
            //    UB = Y X H^(-1) 
            //    UB = H6^T/H7 = (hkl)(UB)(hkl)^T/(hkl)(hkl)^T
            //    H^{-1} = covariance matrix
            Matrix mHI = new Matrix(3);
            double det;
            mH.InvertMatrix(mHI, out det);
            mHI.CopyToMatrix(matCovariance);

            // 5. Least square fitting
            bool returnbool = true;
            if (Math.Abs(det) <= 1.0E-4)
            {
                // Singular matrix
                string errormsg = "Singular Matrix:  LSQ terminated";
                message += String.Format("{0}\n", errormsg);
                Console.WriteLine(errormsg);
                returnbool = false;

            }
            else
            {
                // Least square fitting
                // a. least squae fitting for UB matrix
                mYT.Transpose(mY);
                this.CopyToMatrix(PrevUB);

                //Console.WriteLine("DB312-2:  Matrix H:\n{0}", mH);
                //Console.WriteLine("DB312-3:  Matrix Y:\n{0}", mY);

                Matrix.MulMatrix(mY, mHI, this);
                //Console.WriteLine("UB1130:  Refined  UB Matrix:\n{0}", this);
                //Console.WriteLine("UB1130:  Previous UB Matrix:\n{0}", PrevUB);
                
                // b. Statistic of fitted UB matrix
                UnitCell preUnitCell = new UnitCell();
                mDirectCell.CopyTo(preUnitCell);
                _CalculateUnitCellFromUBMatrix();

                // b1. Calculate chi^2
                _CalculateChi2(reflectsList);

                // b2. Logging
                message += String.Format("Refined UB Matrix with chi^2 = {1}:\n{0}", this, mChi2);
                message += String.Format("Refined Unit Cell: {0}  vs \nPrevious Unit Cell:  {1}",
                    mDirectCell, preUnitCell);
                message += String.Format("Covariance Matrix:\n{0}", matCovariance);

                // c. If # reflection >= 4, then calculate standard deviation and etc. 
                if (nlsq >= 4)
                {
                    // Calculate delta_h_phi between observed and calculated by UB
                    double[] delta1 = new double[3];
                    double[] delta2 = new double[3];

                    for (int i = 0; i < 3; i++)
                    {
                        delta1[i] = 0.0;
                        delta2[i] = 0.0;
                    }
                    Vector hphical = new Vector(3);
                    Vector hphi = new Vector(3);
                    for (int n = 0; n < nlsq; n++)
                    {
                        this.MultVector(reflectsList[n].MillerIndex, hphical);                        
                        INST_ANGLES.VectorFromAngles(reflectsList[n].MotorAngles, mWaveLength, hphi);
                        for (int i = 0; i < 3; i++)
                        {
                            delta1[i] = hphi.V[i] - hphical.V[i];
                            delta2[i] += delta1[i] * delta1[i];
                        }
                    }

                    // double tthc = Math.Asin(Math.Sqrt(sum) * mWaveLength * 0.5) * 114.5916;
                    // b. Calculate uncertainty for UB
                    //    mP:  perturbation matrix
                    Matrix mP = new Matrix(3);
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                        {
                            mP.M[j][i] = Math.Sqrt(Math.Abs(mY.M[i][i]) * delta2[j] / (nlsq - 3)); 
                        }
                    // Console.WriteLine("UB143-2  Matrix Z2 (Pertubation matrix):\n{0}", mP);

                    // c. Store UB matrix
                    Matrix storedUB = new Matrix(3);
                    this.CopyToMatrix(storedUB);

                    Matrix matUBT = new Matrix(3);
                    Matrix matGI = new Matrix(3);
                    Matrix matG = new Matrix(3);
                    UnitCell cDA = new UnitCell();
                    UnitCell cDB = new UnitCell();  // sigma^2
                    UnitCell sigmaCell = new UnitCell();
                    UnitCell storedUnitCell = new UnitCell();
                    mDirectCell.CopyTo(storedUnitCell);

                    // d. Perturb UB matrix for uncertainty
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                        {
                            // i) Perturb UB and calculate G
                            mM[i][j] += mP.M[i][j];
                            this.Transpose(matUBT);
                            Matrix.MulMatrix(matUBT, this, matGI);
                            matGI.InvertMatrix(matG, out det);

                            // ii) calculate direct space volume and lattice constants
                            //     Calculate unit cell from G^-1
                            double vol = 1.0 / Math.Sqrt(det);                            
                            for (int k = 0; k < 3; k++)
                                cDA.Axiles[k] = Math.Sqrt(matG.M[k][k]);
                            cDA.Angles[0] = DasMath.AcosD(matG.M[2][1] / (cDA.Axiles[2] * cDA.Axiles[1]));
                            cDA.Angles[1] = DasMath.AcosD(matG.M[2][0] / (cDA.Axiles[2] * cDA.Axiles[0]));
                            cDA.Angles[2] = DasMath.AcosD(matG.M[0][1] / (cDA.Axiles[0] * cDA.Axiles[1]));

                            for (int k = 0; k < 3; k++)
                            {
                                double tempe = (cDA.Axiles[k] - mDirectCell.Axiles[k]);
                                double tempa = (cDA.Angles[k] - mDirectCell.Angles[k]);
                                cDB.Angles[k] += tempa * tempa;
                                cDB.Axiles[k] += tempe * tempe;
                            }
                            cDB.Volume += (vol - mDirectCell.Volume) * (vol - mDirectCell.Volume);

                            // c restore UB element
                            mM[j][i] = mM[j][i] - mP.M[j][i];

                        } // 300	continue

                    // e. Calculate Sigma
                    for (int k = 0; k < 3; k++)
                    {
                        sigmaCell.Axiles[k] = Math.Sqrt(cDB.Axiles[k]);
                        sigmaCell.Angles[k] = Math.Sqrt(cDB.Angles[k]);
                        // double sigangle = Math.Sqrt(cDB.Angles[k]);
                        // sigmaCell.Angles[k] = DasMath.SinD(mDirectCell.Angles[k]) * sigangle * DasMath.rad;
                    }
                    sigmaCell.Volume = Math.Sqrt(cDB.Volume);

                    // f. Back
                    storedUB.CopyToMatrix(this);
                    storedUnitCell.CopyTo(mDirectCell);

                    // -1: Output
                    if (_DBOUTPUT)
                    {
                        Console.WriteLine("Sigma Of Cell Parameters: \n\t{0}", sigmaCell);
                    }
                }
                else
                {
                    Console.WriteLine("Reflection < 4:  No Statistic otuput");
                }
            }

            // Set up inverted UB matrix
            this.InvertMatrix(m_invertUB, out det);

            return returnbool;
        }


        /// <summary>
        /// Grid refine UB Matrix
        /// </summary>
        /// <param name="reflectionlist">Reflections (angle+motor) for refinement on UB matrix</param>
        /// <param name="motoranglerange"></param>
        /// <param name="motoranglestep"></param>
        /// <param name="bestmotorrangle">[Output] best motor anges</param>
        /// <returns></returns>
        public bool GridRefineUBMatrix(Reflection[] reflectionlist, INST_ANGLES motoranglerange, INST_ANGLES motoranglestep,
            INST_ANGLES bestmotorrangle)
        {
            // Determine the grid
            int num2theta = (int)((2 * motoranglerange.twotheta) / motoranglestep.twotheta) + 1;
            int numomega = (int)(2 * motoranglerange.omega / motoranglestep.omega) + 1;
            int numchi = (int)(2 * motoranglerange.chi / motoranglestep.chi) + 1;

            // TODO - Add a method to store the original UB matrix
            // StoreUBMatrix();
            // InitBestUBMatrix();

            // Refine UB matrix (linearly)
            for (int i = 0; i < num2theta; ++i)
            {
                double delta2theta = -motoranglerange.twotheta + (double)i * motoranglestep.twotheta;
                for (int j = 0; j < numomega; ++j)
                {
                    double deltaomega = -motoranglerange.omega + (double)j * motoranglestep.omega;
                    for (int k = 0; k < numchi; ++k)
                    {
                        double deltachi = -motoranglerange.chi + (double)k * motoranglestep.chi;
                        Console.WriteLine("d(2theta) = {0},  d(omega) = {1},  d(chi) = {2}",
                            delta2theta, deltaomega, deltachi);

                        // TODO - Refine UB matrix

                        // TODO - Check whether need to save this as best UB matrix (UB, motor offset)

                        // TODO - Restore UB matrix and UB^-1
                        // RestoreUBMatrix();

                    }
                }
            }

            // Set up output
            bestmotorrangle.twotheta = 0.1;
            bestmotorrangle.omega = 0.1;
            bestmotorrangle.chi = 0.1;

            // TODO - Implement ASAP Set best UB matrix and UB^-1
            // SetUBMatrix(..., ...);

            return true;
        }



        /// <summary>
        /// List all pairs of angles of the given incident angle
        /// </summary>
        /// <param name="IncidentAngles"></param>
        /// <param name="ncount"></param>
        /// <param name="ang12list"></param>
        /// <param name="pair1list"></param>
        /// <param name="pair2list"></param>
        public void ListPairAngles(INST_ANGLES[] IncidentAngles, int ncount, double[] ang12list,
            int[] pair1list, int[] pair2list)
        {
            Vector uH1 = new Vector(3);
            Vector uH2 = new Vector(3);

            int index = 0;
            for (int n1 = 0; n1 < ncount - 1; n1++)
            {
                INST_ANGLES.UnitVectorFromAngle(IncidentAngles[n1], uH1);
                for (int n2 = n1 + 1; n2 < ncount; n2++)
                {
                    INST_ANGLES.UnitVectorFromAngle(IncidentAngles[n2], uH2);
                    double dotp = Vector.DotProduct(uH1, uH2);
                    pair1list[index] = n1;
                    pair2list[index] = n2;
                    ang12list[index] = DasMath.AcosD(dotp);
                    index++;
                }
            }

            // 2. 3 Sort in increasing order
            int[][] sortservarrays = new int[2][];
            sortservarrays[0] = pair1list;
            sortservarrays[1] = pair2list;
            double[][] dsortservarrays = new double[0][];
            DasMath.Sort(ang12list, sortservarrays, dsortservarrays);

            return;
        }


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

            Matrix mXtx = new Matrix(3);

            // FIXME - ? XTX = X x X
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    mXtx.M[i][j] = 0.0;
                    for (int k = 0; k < 3; k++)
                    {
                        mXtx.M[i][j] += mX.M[k][i] * mX.M[k][j];
                    }
                } // DO 40 J=1,3
            }

            // 2. Form (reciprocal?) unit cell
            for (int i = 0; i < 3; i++)
            {
                cA2.Axiles[i] = Math.Sqrt(mXtx.M[i][i]);
            } // DO 50 I=1,3

            cA2.alpha = mXtx.M[1][2] / (cA2.b * cA2.c);
            cA2.beta = mXtx.M[2][0] / (cA2.c * cA2.a);
            cA2.gamma = mXtx.M[0][1] / (cA2.a * cA2.b);

            cA2.CalculateReciprocalUnitCell(cA1);

            // recip(a2,vol2,a1,vol1);  FIXME:  what is this for???
            for (int i = 0; i < 3; i++)
            {
                double tx = Math.Sqrt(1.0 - cA1.Angles[i] * cA1.Angles[i]);
                double ty = cA1.Angles[i];
                angles[i] = DasMath.rad * DasMath.Atan2(tx, ty);
            }
            return;
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
            int nc = nrank + ncv;
            det = 1.0;
            double[] dA = new double[256];
            int ij, kk;
            double pivot;
#if true
            ij = 0;
            kk = 0;
            pivot = 0;
            Console.WriteLine("evaluating to ij, kk, pivot is not true");
#endif

            //  START MASTER REDUCTION	LOOP
            int m, l, jsub, lj, kj;
            double temp;
            for (int k = 0; k < nrank; k++)
            {
                // DO 240 K=1,NRANK
                double piv = 0.0;
                int ksub = ida * (k - 1);
                // START SEARCH FOR PIVOT
                for (int j = k; j < nrank; j++)
                { // 110
                    jsub = ida * (j - 1);
                    for (int i = k; i < nrank; i++)
                    { // do 110 
                        ij = jsub + i;
                        // check if abs(a(i,j)).gt.abs(pivot)
                        double p = Math.Abs(mA.M[i][j]);
                        if ((p - piv) > 0)
                        {
                            // found new value for trial pivot 
                            piv = p;
                            pivot = mA.M[i][j];
                            m = j;
                            l = i;
                        }
                    }
                } // 110 
                //  CHECK SIZE OF PIVOT
#if true
                Console.WriteLine("The evaluating to eps and pivot is wrong!");
                eps = 0;
                pivot = 0;
                l = 0;
                m = 0;
                kj = 0;
                kk = 0;
#endif
                if ((piv - eps) <= 0)
                {
                    // 120,130,130	//
                    // SINGULARITY RETURN 
                    jrank = k - 1;
                    eps = 0;
                    return;
                }
                else
                {
                    //  VALID PIVOT FOUND, UPDATE DETERMINANT
                    det = det * pivot;
                    //  SAVE ORIGINAL ROW AND COLUMN OF PIVOT 
                    ntemp[k] = 256 * l + m;
                }
                if ((l - k) != 0)
                { // 140,160,140
                    //  NEED TO PUT PIVOT ELEMENT INTO	A(K,K),	INTERCHANGE ROWS
                    det = -det;
                    for (int j = 0; j < nc; j++)
                    { // DO 150 J=1,NC
                        jsub = ida * (j - 1);
                        lj = l + jsub;
                        kj = k + jsub;
                        temp = mA.M[l][j];
                        mA.M[l][j] = mA.M[k][j];
                        mA.M[k][j] = temp;
                    } //150 

                }
                // 160 
                if ((m - k) != 0)
                {  // 170,190,170
                    //  interchange columns
                    det = -det;
                    int jm = ida * (m - 1);
                    int jk = ksub;
                    for (int j = 0; j < nrank; j++)
                    { // do 180 j
                        jm = jm + 1;
                        jk = jk + 1;
                        temp = mA.M[j][m];
                        mA.M[j][m] = mA.M[j][k];
                        mA.M[j][k] = temp;
                    } //180 

                    //  reduce	pivot row
                } // 190 
                pivot = 1.0 / pivot;
                for (int j = 0; j < nc; j++)
                { // do 200 j
                    kj = ida * (j - 1) + k;
                    mA.M[k][j] = pivot * mA.M[k][j];
                } // 200 

                kk = k + ksub;
                mA.M[k][k] = 0.0;
                //  REDUCE	NON-PIVOT ROWS
                for (int i = 0; i < nrank; i++)
                { // DO 230 I=1,NRANK
                    int ik = i + ksub;
                    temp = mA.M[i][k];
                    if ((i - k) != 0)
                    { // 210,230,210
                        // 210 
                        for (int j = 0; j < nc; j++)
                        { // DO 220 J
                            ij = ida * (j - 1) + i;
                            kj = ida * (j - 1) + k;
                        } // 220 
                        dA[ij] = dA[ij] - temp * dA[kj];
                        dA[ik] = -temp * pivot;
                    }
                } // 230
            } // 	240 
            dA[kk] = pivot;

            //  INVERSION COMPLETE, RESTORE ORIGINAL ROW AND COLUMN ORDER
            jrank = nrank;
            int kx = nrank;
            for (int j = 0; j < nrank; j++)
            { // do 300 j
                int ksub = ida * (kx - 1);
                m = ntemp[kx] / 256;
                // l=mod(ntemp(kx),256);
                l = ntemp[kx] % 256;
                if ((m - kx) != 0)
                { // 250,270,250 //
                    // INTERCHANGE ROWS 
                    int msub = ida * (m - 1);
                    for (int i = 0; i < nrank; i++)
                    { // DO 260 I
                        int ik = ksub + i;
                        int im = msub + i;
                        temp = dA[ik];
                        dA[ik] = dA[im];
                        dA[im] = temp;
                    } // 260 
                } // 270 
                if ((l - kx) != 0)
                { // 280,300,280
                    // INTERCHANGE COLUMNS
                    for (int i = 0; i < nc; i++)
                    { //280 DO 290
                        int ki = ida * (i - 1) + kx;
                        int li = ida * (i - 1) + l;
                        temp = dA[ki];
                        dA[ki] = dA[li];
                        dA[li] = temp;
                    } // 290 
                }
                kx = kx - 1;
            } // 300 

            eps = 0;
#if true
            Console.WriteLine("eps = 0 is faked");
#endif
            return;

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

            throw new Exception("Will be implemented in Phase-2");

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
                INST_ANGLES.VectorFromAngles(reflectslist[i].MotorAngles, mWaveLength, hphis[i]);
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
            mDirectCell.CalculateFromG(mG);

            // 6.  C C  BEGIN CONSTRAINED LEAST-SQUARES  construct the true orientation matrix bb corresponding to the unconstrained cell parameters

            // c C  FORM AND INVERT BB AND OBTAIN THE TRUE ORTHOGONAL ROTATION U FROM X c
            // <source>CALL BCALC(BB,A1,A2)</source>
            Matrix mBI = new Matrix(3);
            _CalculateUnitCellFromUBMatrix();
            mB.InvertMatrix(mBI, out det);
            Matrix.MulMatrix(mX, mBI, mU);



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

        #endregion

        #region static public method
        /// <summary>
        /// Routine to select the two trial reflections for generation of the orientation matrix
        /// Usually the first two are selected
        /// </summary>
        /// <source>subroutine RefPick(ncount,n1,n2,ang12)</source>
        /// <param name="IncidentAngles">Array of incident angles</param>
        /// <param name="n1">first pickup</param>
        /// <param name="n2">second pickup</param>
        /// <param name="ang12">angles b/w two pickups in coordinator (phi)</param>
        /// <returns></returns>
        static public bool RefPick(INST_ANGLES[] IncidentAngles, out int n1, out int n2, out double ang12)
        {
            // integer*4 ipoint
            // common /useref/ipoint(listlen)
            // real*4 h1(3),h2(3)

            int ncount = IncidentAngles.Length;

            Vector vH1 = new Vector(3);
            Vector vH2 = new Vector(3);
            Vector uH1 = new Vector(3);
            Vector uH2 = new Vector(3);

            n1 = -1;
            n2 = -1;

            bool found = false;
            ang12 = 0.0;
            for (int t1 = 0; t1 < ncount - 1; t1++)
            {
                // 1. Get observed angle between the two reflections. 
                INST_ANGLES.UnitVectorFromAngle(IncidentAngles[t1], vH1);
                vH1.UnitVector(uH1);

                // Loop over n2
                for (int t2 = t1 + 1; t2 < ncount; t2++)
                {
                    // 2. Get observed angle between the two reflections
                    INST_ANGLES.UnitVectorFromAngle(IncidentAngles[t2], vH2);
                    vH2.UnitVector(uH2);

                    double dot = Vector.DotProduct(uH1, uH2);
                    ang12 = DasMath.AcosD(dot);

                    if (ang12 > 20)
                    {
                        found = true;
                        n1 = t1;
                        n2 = t2;
                        break;
                    }
                }

                // break if found
                if (found == true)
                {
                    break;
                }
            } // for i

            if (found == false)
            {
                ang12 = -1.0E8;
            }

            return found;

        }
        #endregion

        #region Unused
#if false
        /// <summary>
        /// Calcualte (Miller) indices corresponding to a (instrument) angle
        /// <source>UBCALS.FOR -> HCALC()</source>
        /// </summary>
        /// <param name="angl">Instrument angle</param>
        /// <returns>Miller indices</returns>
        public MillerIndices CalHc(INST_ANGLES angl)
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
            INST_ANGLES.VectorFromAngles(angl, mWaveLength, h_phi);

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
            INST_ANGLES.VectorFromAngles(r1.motor_angles, mWaveLength, u_phi1);
            u_phi2 = new Vector(3);
            INST_ANGLES.VectorFromAngles(r2.motor_angles, mWaveLength, u_phi2);

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
        #endregion
    }
}
