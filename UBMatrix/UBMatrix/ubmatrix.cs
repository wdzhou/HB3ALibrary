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
        /// <summary>
        /// U matrix
        /// </summary>
        private Matrix mU;
        /// <summary>
        /// B matrix
        /// </summary>
        private Matrix mB;
        /// <summary>
        /// G matrix
        /// </summary>
        private Matrix mG;
        /// <summary>
        /// Inverted G matrix
        /// </summary>
        private Matrix mGI;
        /// <summary>
        /// inversed UB 
        /// </summary>
        private Matrix m_invertUB;
        /// <summary>
        /// Covariance matrix
        /// </summary>
        private Matrix m_covarMatrix;
        /// <summary>
        /// Direct unit cell
        /// </summary>
        private UnitCell m_directUnitCell;
        /// <summary>
        /// Reciprocal unit cell
        /// </summary>
        private UnitCell m_reciprocalUnitCell;
        /// <summary>
        /// Stored direct unit cell
        /// </summary>
        private UnitCell m_storedDirectUnitCell;
        /// <summary>
        /// UB matrix calculation counter
        /// </summary>
        // private int mCalUBCounter;
        private double[] m_arrTolerance;
        /// <summary>
        /// Chi square
        /// </summary>
        private double m_chi2; 

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
            // Matrix initialization
            mU = new Matrix(3);
            mB = new Matrix(3);
            m_invertUB = new Matrix(3);
            mG = new Matrix(3);
            mGI = new Matrix(3);
            m_covarMatrix = new Matrix(3);

            PrevUB = new Matrix(3);

            // Wavelength and Unit cell
            mWaveLength = -1; // unphysical value
            m_arrTolerance = new double[2];
            for (int i = 0; i < 2; i++)
                m_arrTolerance[i] = 1.0;
            m_directUnitCell = new UnitCell();
            m_reciprocalUnitCell = new UnitCell();
            m_storedDirectUnitCell = new UnitCell();

            m_chi2 = -1.0;

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
            get { return m_covarMatrix; }
        }
        /// <summary>
        /// Unit cell (direct space)
        /// </summary>
        public UnitCell CrystallCell
        {
            get { return m_directUnitCell; }
            set { m_directUnitCell = value; }
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
            get { return m_chi2; }
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
        private void _MakeUB(MillerIndices vH1, MillerIndices vH2, MotorIncidentAngles ang1, MotorIncidentAngles ang2, out int ierr)
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
            MotorIncidentAngles.UnitVectorFromAngle(ang1, h1p);
            MotorIncidentAngles.UnitVectorFromAngle(ang2, h2p);
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
            MotorIncidentAngles.UnitVectorFromAngle(ang1, h1p);
            MotorIncidentAngles.UnitVectorFromAngle(ang2, h2p);

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
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix._MakeUB(). " + logmessage, "Info");
            return;
        }

        /// <summary>
        /// Base step to calculate error of a UB matrix to an observed reflection
        /// </summary>
        /// <param name="iangle"></param>
        /// <param name="mindex"></param>
        /// <param name="diffangle"></param>
        /// <param name="cal2theta"></param>
        private void _CalculateError(MotorIncidentAngles iangle, MillerIndices mindex, out double diffangle, out double cal2theta)
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
                m_arrTolerance[index] = value;
            }
            else
            {
                m_arrTolerance[index] = 1.0;
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
            m_reciprocalUnitCell.CalculateFromG(mGI);
            m_reciprocalUnitCell.CalculateReciprocalUnitCell(m_directUnitCell);
            m_directUnitCell.CalculateG(mG);

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
                MotorIncidentAngles.VectorFromAngles(tReflections[i].MotorAngles, mWaveLength, uphiobs);
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

            m_chi2 = diff2sum;
            return diff2sum;
        }

        #endregion

        #region public methods

        /// <summary>
        /// Save UB matrix to an XML file, including
        /// 1. direct cell
        /// 2. ub matrix
        /// </summary>
        /// <param name="pathfilename"></param>
        public void SaveToXMLFile(string pathfilename)
        {

            // Initialize XML writer
            XmlTextWriter xmlTestWriter = new XmlTextWriter(pathfilename, System.Text.Encoding.UTF8);
            xmlTestWriter.Formatting = Formatting.Indented;
            xmlTestWriter.WriteStartDocument(false);

            // Comment
            xmlTestWriter.WriteComment("UB Matrix Information");
            // Start element
            xmlTestWriter.WriteStartElement("ubmatrix");
            // time                
            xmlTestWriter.WriteStartElement("time");
            DateTime currtime = DateTime.Now;
            xmlTestWriter.WriteAttributeString("time", currtime.ToString());
            xmlTestWriter.WriteEndElement();
            // unit cell              
            xmlTestWriter.WriteStartElement("unitcell");
            xmlTestWriter.WriteAttributeString("a", m_directUnitCell.a.ToString());
            xmlTestWriter.WriteAttributeString("b", m_directUnitCell.b.ToString());
            xmlTestWriter.WriteAttributeString("c", m_directUnitCell.c.ToString());
            xmlTestWriter.WriteAttributeString("alpha", m_directUnitCell.alpha.ToString());
            xmlTestWriter.WriteAttributeString("beta", m_directUnitCell.beta.ToString());
            xmlTestWriter.WriteAttributeString("gamma", m_directUnitCell.gamma.ToString());
            xmlTestWriter.WriteEndElement();
            // wavelength
            xmlTestWriter.WriteStartElement("wavelength");
            xmlTestWriter.WriteAttributeString("lambda", mWaveLength.ToString());
            xmlTestWriter.WriteEndElement();
            // ub matrix
            xmlTestWriter.WriteStartElement("matrix");
            xmlTestWriter.WriteAttributeString("matrix", XmlUtility.ConvertToString(this.M));
            xmlTestWriter.WriteEndElement();
            // End element
            xmlTestWriter.WriteEndElement();

            // Close
            xmlTestWriter.Flush();
            xmlTestWriter.Close();

            return;
        }


        /// <summary>
        /// Load from an XML File
        /// </summary>
        /// <param name="pathfilename"></param>
        public void LoadFromXMLFile(string pathfilename)
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

            if (_DBOUTPUT)
            {
                Console.WriteLine("Loading UnitCell:\n\t> {0}", uc);
                Console.WriteLine("\t> Lambda = {0}\n > Matrix: \n", lambda);
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                        Console.Write("{0:0.000000}  ", ubmatrix[i][j]);
                    Console.WriteLine();
                }
            }
            
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
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.LoadFromXMLFile(). " + logmessage, "Info");

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
            inpcell.CopyTo(m_directUnitCell);

            // 2. Calculate G matrix
            m_directUnitCell.CalculateG(mG);

            // 3. Set up reciprocal lattice
            m_directUnitCell.CalculateReciprocalUnitCell(m_reciprocalUnitCell);
            m_reciprocalUnitCell.CalculateG(mGI);

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
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.SetDirectUnitCell(). " + String.Format("Set Wave-length = {0}", mWaveLength.ToString()), "Info");
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.SetDirectUnitCell(). " + String.Format("Set Unit Cell:\n{0}", m_directUnitCell), "Info");
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.SetDirectUnitCell(). " + String.Format("Set Reciprocal Cell:\n{0}", m_reciprocalUnitCell), "Info");
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.SetDirectUnitCell(). " + String.Format("Set G Matrix: \n{0}", mG), "Info");
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.SetDirectUnitCell(). " + String.Format("set G^-1 Matrix:\n{0}", mGI), "Info");
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.SetDirectUnitCell(). " + String.Format("Set B matrix:\n{0}", mB), "Info");

            return;
        }

        /// <summary>
        /// Get direct cell information
        /// </summary>
        /// <param name="edges"></param>
        /// <param name="angles"></param>
        public void GetDirectUnitCell(double[] edges, double[] angles)
        {
            string logmessage = String.Format("UB Matrix:  Unit Cell = {0}", m_directUnitCell);
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.GetDirectUnitCell(). " + logmessage, "Info");

            for (int i = 0; i < 3; i++)
                edges[i] = m_directUnitCell.Axiles[i];
            for (int i = 0; i < 3; i++)
                angles[i] = m_directUnitCell.Angles[i];
            
            return;
        }

        /// <summary>
        /// Sotre the current Unit Cell
        /// </summary>
        public void StoreUnitCell()
        {
            m_directUnitCell.CopyTo(m_storedDirectUnitCell);
            return;
        }

        /// <summary>
        /// Restore the Unit Cell and related objects
        /// </summary>
        public void RestoreUnitCell()
        {
            SetDirectUnitCell(m_storedDirectUnitCell);
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

            mB.M[0][0] =  m_reciprocalUnitCell.a;
            mB.M[0][1] =  m_reciprocalUnitCell.b * DasMath.CosD(m_reciprocalUnitCell.gamma);
            mB.M[1][1] =  m_reciprocalUnitCell.b * DasMath.SinD(m_reciprocalUnitCell.gamma);
            mB.M[0][2] =  m_reciprocalUnitCell.c * DasMath.CosD(m_reciprocalUnitCell.beta);
            mB.M[1][2] = -m_reciprocalUnitCell.c * DasMath.SinD(m_reciprocalUnitCell.beta) * DasMath.CosD(m_directUnitCell.alpha);
            mB.M[2][2] = 1 / m_directUnitCell.c;

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
        public void MakeUB(MillerIndices vH1, MillerIndices vH2, MotorIncidentAngles ang1, MotorIncidentAngles ang2, out int ierr)
        {
            // 1. Calculate UB matrix from 2 reflections
            _MakeUB(vH1, vH2, ang1, ang2, out ierr);
            double det;
            InvertMatrix(m_invertUB, out det);

            // 2. Logging
            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.MakeUB(). " + String.Format("UB Matrix Calculated:\n{0}\nFrom 1: {1}  {2}\n     2: {3}  {4}\n",
                this, vH1, ang1, vH2, ang2), "Info");

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
        public void CalMotorAnglesFromMillerIndex(MillerIndices inphkl, MotorIncidentAngles angl, bool fixed_phi, out int ierror)
        {
            // Calculate h_phi
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

            Log.Write(Log.Info, "In DAS.UBMatrix.UBMatrix.CalMotorAnglesFromMillerIndex(). " + String.Format("Calculate Motor Angle:   {0} ==> {1}\n", inphkl, angl), "Info");

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
                MotorIncidentAngles.VectorFromAngles(refls[i].MotorAngles, mWaveLength, vHs[i]);

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
            // mCalUBCounter++;

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
        public void CalculateHKLFromIncidentAngle(MotorIncidentAngles iangle, MillerIndices hkl, Vector uPhi)
        {
            MotorIncidentAngles.UnitVectorFromAngle(iangle, uPhi);
            double dt = mWaveLength / (2.0 * DasMath.SinD(0.5 * iangle.twotheta));
            Vector vHH = m_invertUB.MultVector(uPhi);

            // Console.WriteLine("[CalHKL] UB Matrix: \n{0}", this);
            // Console.WriteLine("[CalHKL] UB^-1    : \n{0}", m_invertUB);

            // Console.WriteLine("Input motor angle: {0}. ", iangle);

            if (_DBOUTPUT)
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
        public bool CheckError(MotorIncidentAngles iangle, MillerIndices mindex, double toleranceiangle, double tolerancemindex,
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
        public bool TestUBByReflections(MotorIncidentAngles[] IncidentAngles, MillerIndices[] calHKLs, 
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
                if ( diff2theta[icount] <= m_arrTolerance[0] )
                {
                    if (angdot < m_arrTolerance[1])
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
        /// VZ's version to generate UB matrix from a series of reflections and 
        /// known lattice parameters
        /// The best fitting UB matrix is stored in mM[][], while
        /// the next 10 best fitting UB matrices are stored in a list too. 
        /// </summary>
        /// <param name="IncidentAngles"></param>
        /// <param name="numreflections">int, number of reflections (in case input list is larger)</param>
        /// <returns></returns>
        public bool GenerateUBFromUnIndexedReflectiuons(MotorIncidentAngles[] IncidentAngles, int numreflections)
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
                        m_arrTolerance[0])
                    {

                        // n2:  Guess HKL for Second reflection
                        for (int jhkl = 0; jhkl < sizehklist; jhkl++)
                        {
                            MillerIndices hkl2 = hkllist[jhkl];
                            double d2;
                            double twotheta_cal_2 = Calculate2ThetaFromHKL(hkl2, out d2);

                            if (Math.Abs(IncidentAngles[n2].twotheta - twotheta_cal_2) <
                                m_arrTolerance[0])
                            {
                                // d) 3rd Check: Is the angle between the relections correct?
                                double da1 = Vector.DotProduct(hkl1, mGI.MultVector(hkl2));
                                double da2 = DasMath.AcosD(da1 / Math.Sqrt(d1 * d2));

                                if (Math.Abs(da2 - ang12list[pairindex]) < m_arrTolerance[1])
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
                MotorIncidentAngles tempmotorangle = reflecslist[i].MotorAngles;

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
                Log.Write(Log.Error, String.Format("In DAS.UBMatrix.UBMatrix.RefineUBLinearLSF(). Refining UB Matrix Error Message: {0}", errmsg), "Error");
                return false;
            }
            else if (nlsq < reflectsList.Length)
            {
                string errmsg = String.Format("Input number of reflections ({0}) is smaller than 'nlsq' ({1}), {2}.",
                    reflectsList.Length, nlsq, "which will cause index error");
                Console.WriteLine(errmsg);
                Log.Write(Log.Error, String.Format("In DAS.UBMatrix.UBMatrix.RefineUBLinearLSF(). Refining UB Matrix Error Message: {0}", errmsg), "Error");
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
                MotorIncidentAngles.VectorFromAngles(reflectsList[i].MotorAngles, mWaveLength, hphiList[i]);
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
            mHI.CopyToMatrix(m_covarMatrix);

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
                m_directUnitCell.CopyTo(preUnitCell);
                _CalculateUnitCellFromUBMatrix();

                // b1. Calculate chi^2
                _CalculateChi2(reflectsList);

                // b2. Logging
                message += String.Format("Refined UB Matrix with chi^2 = {1}:\n{0}", this, m_chi2);
                message += String.Format("Refined Unit Cell: {0}  vs \nPrevious Unit Cell:  {1}",
                    m_directUnitCell, preUnitCell);
                message += String.Format("Covariance Matrix:\n{0}", m_covarMatrix);

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
                        MotorIncidentAngles.VectorFromAngles(reflectsList[n].MotorAngles, mWaveLength, hphi);
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
                    m_directUnitCell.CopyTo(storedUnitCell);

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
                                double tempe = (cDA.Axiles[k] - m_directUnitCell.Axiles[k]);
                                double tempa = (cDA.Angles[k] - m_directUnitCell.Angles[k]);
                                cDB.Angles[k] += tempa * tempa;
                                cDB.Axiles[k] += tempe * tempe;
                            }
                            cDB.Volume += (vol - m_directUnitCell.Volume) * (vol - m_directUnitCell.Volume);

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
                    storedUnitCell.CopyTo(m_directUnitCell);

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
        /// Store present UB matrix and inverted UB matrix 
        /// </summary>
        /// <param name="storedUB"></param>
        /// <param name="storeInvUB"></param>
        private void StoreUBMatrix(Matrix storedUB, Matrix storeInvUB)
        {
            this.CopyToMatrix(storedUB);
            m_invertUB.CopyToMatrix(storeInvUB);

            return;
        }

        /// <summary>
        /// Restore to present UB matrix and inverted UB matrix
        /// </summary>
        /// <param name="storedUB"></param>
        /// <param name="storedInvUB"></param>
        private void RestoreUBMatrix(Matrix storedUB, Matrix storedInvUB)
        {
            storedUB.CopyToMatrix(this);
            storedInvUB.CopyToMatrix(m_invertUB);

            return;
        }


        /// <summary>
        /// Refine UB matrix on a 3D grid of 2theta, omega and chi in order to find the 
        /// motors offsets with smallest UB matrix error
        /// </summary>
        /// <param name="reflectionlist"></param>
        /// <param name="motoranglerange"></param>
        /// <param name="motoranglestep"></param>
        /// <param name="bestMotorShift">[Output] best motor anges</param>
        /// <returns></returns>
        public bool GridRefineUBMatrix(Reflection[] reflectionlist, MotorIncidentAngles motoranglerange, MotorIncidentAngles motoranglestep,
            MotorIncidentAngles bestMotorShift)
        {
            // Determine the grid
            int num2theta = (int)((2 * motoranglerange.twotheta) / motoranglestep.twotheta) + 1;
            int numomega = (int)(2 * motoranglerange.omega / motoranglestep.omega) + 1;
            int numchi = (int)(2 * motoranglerange.chi / motoranglestep.chi) + 1;

            // Store current UB matrix
            Matrix storedUB = new Matrix(3, 3);
            Matrix storedInvUB = new Matrix(3, 3);
            Matrix bestUB = new Matrix(3, 3);
            Matrix bestInvUB = new Matrix(3, 3);
            StoreUBMatrix(storedUB, storedInvUB);

            bool restore = false;

            // Initialize best 
            double bestError = 1.0E100;
            double deltaphi = 0;

            DateTime starttime = DateTime.Now;

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


                        // Use the original UB matrix
                        if (restore)
                        {
                            RestoreUBMatrix(storedUB, storedInvUB);
                        }
                        else
                        {
                            restore = true;
                        }
                        
                        // Refine UB matrix
                        Reflection[] vecShiftedReflects = GenerateReflectionOffsetMotorPositions(reflectionlist, delta2theta,
                            deltaomega, deltachi, deltaphi);
                        bool refinegood = RefineUBLinearLSF(vecShiftedReflects, vecShiftedReflects.Length);

                        // Check whether need to save this as best UB matrix (UB, motor offset)
                        double error = CalculateUBMatrixError(vecShiftedReflects);
                        // Console.WriteLine("d(2theta) = {0},  d(omega) = {1},  d(chi) = {2}:  UB Error = {3}", 
                        //     delta2theta, deltaomega, deltachi, error);
                        if (error < bestError)
                        {
                            bestError = error;
                            bestMotorShift.twotheta = delta2theta;
                            bestMotorShift.omega = deltaomega;
                            bestMotorShift.chi = deltachi;

                            StoreUBMatrix(bestUB, bestInvUB);
                        }
                    }
                }
            }

            DateTime endtime = DateTime.Now;
            TimeSpan duration = endtime.Subtract(starttime);
            Console.WriteLine("Grid refinement spent {0}. Best error = {1}", duration, bestError);

            // Use best UB matrix and inverted UB matrix
            RestoreUBMatrix(bestUB, bestInvUB);

            return true;
        }

        /// <summary>
        /// Generate a new list of reflections with motor angle shifted
        /// </summary>
        /// <param name="reflectionlist"></param>
        /// <param name="delta2thata"></param>
        /// <param name="deltaomega"></param>
        /// <param name="deltachi"></param>
        /// <param name="deltaphi"></param>
        /// <returns></returns>
        public Reflection[] GenerateReflectionOffsetMotorPositions(Reflection[] reflectionlist,
            double delta2thata, double deltaomega, double deltachi, double deltaphi)
        {
            // Initialize data
            int numreflects = reflectionlist.Length;
            Reflection[] offsetreflist = new Reflection[numreflects];

            for (int i = 0; i < numreflects; ++i)
            {
                // Calculate new incident angles
                MotorIncidentAngles tempangles = new MotorIncidentAngles();
                tempangles.twotheta = reflectionlist[i].MotorAngles.twotheta + delta2thata;
                tempangles.omega = reflectionlist[i].MotorAngles.omega + deltaomega;
                tempangles.chi = reflectionlist[i].MotorAngles.chi + deltachi;
                tempangles.phi = reflectionlist[i].MotorAngles.phi + deltaphi;

                // Set up new reflections
                offsetreflist[i] = new Reflection();
                offsetreflist[i].MotorAngles = tempangles;
                offsetreflist[i].MillerIndex = reflectionlist[i].MillerIndex;
            }

            return offsetreflist;
        }


        /// <summary>
        /// List all pairs of angles of the given incident angle
        /// </summary>
        /// <param name="IncidentAngles"></param>
        /// <param name="ncount"></param>
        /// <param name="ang12list"></param>
        /// <param name="pair1list"></param>
        /// <param name="pair2list"></param>
        public void ListPairAngles(MotorIncidentAngles[] IncidentAngles, int ncount, double[] ang12list,
            int[] pair1list, int[] pair2list)
        {
            Vector uH1 = new Vector(3);
            Vector uH2 = new Vector(3);

            int index = 0;
            for (int n1 = 0; n1 < ncount - 1; n1++)
            {
                MotorIncidentAngles.UnitVectorFromAngle(IncidentAngles[n1], uH1);
                for (int n2 = n1 + 1; n2 < ncount; n2++)
                {
                    MotorIncidentAngles.UnitVectorFromAngle(IncidentAngles[n2], uH2);
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
        static public bool RefPick(MotorIncidentAngles[] IncidentAngles, out int n1, out int n2, out double ang12)
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
                MotorIncidentAngles.UnitVectorFromAngle(IncidentAngles[t1], vH1);
                vH1.UnitVector(uH1);

                // Loop over n2
                for (int t2 = t1 + 1; t2 < ncount; t2++)
                {
                    // 2. Get observed angle between the two reflections
                    MotorIncidentAngles.UnitVectorFromAngle(IncidentAngles[t2], vH2);
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

    }
}
