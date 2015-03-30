using System;
using System.IO;
using System.Collections.Generic;
using System.Text;

namespace DAS.UBMatrix
{
    /// <summary>
    /// Miller incides for a reflection
    /// </summary>
    public class MillerIndices:Vector
    {
        #region private properties
        double[] mOriginalInput;
        bool mInputIsCalculated;
        double m_Error;
        #endregion

        #region Constructors
        /// <summary>
        /// Constructor of an (HKL)
        /// </summary>
        public MillerIndices():base(3)
        {
            mOriginalInput = new double[3];
            mInputIsCalculated = false;
            m_Error = 0.0;
            return;
        }
        #endregion

        #region public properties
        /// <summary>
        /// Error of calculated (HKL)
        /// </summary>
        public double Error
        {
            get { return m_Error; }
        }
        /// <summary>
        /// H
        /// </summary>
        public double H
        {
            get { return mVector[0]; }
            set { mVector[0] = value; }
        }
        /// <summary>
        /// K
        /// </summary>
        public double K
        {
            get { return mVector[1]; }
            set { mVector[1] = value; }
        }
        /// <summary>
        /// L
        /// </summary>
        public double L
        {
            get { return mVector[2]; }
            set { mVector[2] = value; }
        }
        /// <summary>
        /// Original input HKL
        /// </summary>
        public double[] OriginalHKL
        {
            get { return mOriginalInput; }
        }
        #endregion

        #region public methods
        /// <summary>
        /// Set the (HKL) from a caculated value, i.e., may not be correct
        /// </summary>
        /// <param name="inputhkl"></param>
        public void SetFromCalculatedValue(Vector inputhkl)
        {
            mInputIsCalculated = true;

            double diff = 0.0;
            for (int i = 0; i < 3; i++)
            {
                mOriginalInput[i] = inputhkl.V[i];
                mVector[i] = DasMath.NearestInteger(inputhkl.V[i]);
                diff += (mOriginalInput[i] - mVector[i]) * (mOriginalInput[i] - mVector[i]);                         
            }
            m_Error = Math.Sqrt(diff);

            return;
        }

        /// <summary>
        /// Convert the values to nearest integers (and cast back to double)
        /// </summary>
        public void ConvertToIntegers()
        {
            for (int i = 0; i < 3; ++i)
            {
                mVector[i] = (double)(Convert.ToInt32(mVector[i]));
            }
        }

        /// <summary>
        /// Customized format output
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string rst = "(";
            for (int i = 0; i < 3; i++)
                rst += String.Format("{0,2:0.}  ", mVector[i]);
            if (mInputIsCalculated == true)
            {
                rst += ") <---- (";
                for (int i = 0; i < 3; i++)
                    rst += String.Format("{0,7:0.000}  ", mOriginalInput[i]);
                rst += String.Format(")  dev = {0,7:0.00}", m_Error);
            } else 
                rst += ")";
            return rst;
        }

        /// <summary>
        /// Evaluate whether two Miller Indices are of same values
        /// </summary>
        /// <param name="m2"></param>
        /// <returns></returns>
        public bool EqualValue(MillerIndices m2)
        {
            double tolerance = 1.0E-6;

            bool same = true;
            for (int i = 0; i < 3; i++)
                if (Math.Abs(V[i] - m2.V[i]) > tolerance)
                    same = false;

            return same;
        }


        /// <summary>
        /// Evaluate whether two Miller Indices are of same values
        /// </summary>
        /// <param name="m2"></param>
        /// <returns></returns>
        public bool EqualIntegerValue(MillerIndices m2)
        {
            bool same = true;
            for (int i = 0; i < 3; i++)
                if (Convert.ToInt32(V[i]) != Convert.ToInt32(m2.V[i]))
                    same = false;

            return same;
        }


        /// <summary>
        /// Test whether it is a (0, 0, 0) HKL
        /// </summary>
        /// <returns></returns>
        public bool ZeroHKL()
        {
            bool nullhkl = true;
            for (int i = 0; i < 3; i++)
                if (Math.Abs(V[i]) > 1.0E-6)
                    nullhkl = false;

            return nullhkl;
        }

        /// <summary>
        /// Normalize a Miller index with G^-1 matrix
        /// </summary>
        /// <param name="normindex"></param>
        /// <param name="GInv"></param>
        public void Normalize(MillerIndices normindex, Matrix GInv)
        {
            // 1. get length
            Vector temp1 = new Vector(3);
            GInv.MultVector(this, temp1);
            double length = Math.Sqrt(Vector.DotProduct(this, temp1));

            // 2. Normalize
            for (int i = 0; i < 3; i++)
                normindex.V[i] = this.V[i] / length;

            return;
        }

        #endregion

    }

    /// <summary>
    /// Motor incident angle
    /// </summary>
    public class MotorIncidentAngles
    {
        #region private properties
        double[] mAngles;
        #endregion

        #region Constructors
        /// <summary>
        /// constructor
        /// </summary>
        public MotorIncidentAngles()
        {
            mAngles = new double[4];
        }
        #endregion

        #region public properties
        /// <summary>
        /// 2-theta
        /// </summary>
        public double twotheta
        {
            get { return mAngles[0]; }
            set { mAngles[0] = value; }
        }
        /// <summary>
        /// angle omega... (2) in SINGLE
        /// </summary>
        public double omega
        {
            get { return mAngles[1]; }
            set { mAngles[1] = value; }
        }
        /// <summary>
        /// angle chi
        /// </summary>
        public double chi
        {
            get { return mAngles[2]; }
            set { mAngles[2] = value; }
        }
        /// <summary>
        /// angle phi
        /// </summary>
        public double phi
        {
            get { return mAngles[3]; }
            set { mAngles[3] = value; }
        }        
        #endregion

        #region public method
        /// <summary>
        /// override customized format output
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string msg = String.Format("2theta: {0,8:0.0000};  omega: {1,8:0.0000};  chi: {2,8:0.0000};  phi: {3,8:0.0000}",
                twotheta, omega, chi, phi);
            return msg;
        }

        /// <summary>
        /// Copy the current value to antoher incident angle
        /// </summary>
        /// <param name="angle2"></param>
        public void CopyTo(MotorIncidentAngles angle2)
        {
            for (int i = 0; i < mAngles.Length; i++)
                angle2.mAngles[i] = mAngles[i];
            return;
        }

        /// <summary>
        /// See if all angles are zero
        /// </summary>
        /// <returns></returns>
        public bool ZeroAngle()
        {
            bool zeroangle = true;

            if (Math.Abs(twotheta) > 1.0E-6)
                zeroangle = false;
            else if (Math.Abs(omega) > 1.0E-6)
                zeroangle = false;
            else if (Math.Abs(chi) > 1.0E-6)
                zeroangle = false;
            else if (Math.Abs(phi) > 1.0E-6)
                zeroangle = false;

            return zeroangle;
        }
        #endregion

        #region static public method

        /// <summary>
        /// Calculate h_phi from incident angle omega, gamma, psi
        /// Same as Equation (22) in the Busing paper: 
        /// ICalcH --> Incident to calculate u_pHi
        /// Mapping: 2-> omega, 3 --> chi, 4 --> phi
        ///</summary>
        ///<source>SUBROUTINE ICalcH(Angl,hphi)</source>
        /// <param name="angl">Incident angle</param>
        /// <param name="uphi">Output u_phi</param>
        static public void UnitVectorFromAngle(MotorIncidentAngles angl, Vector uphi)
        {
            uphi.i = DasMath.CosD(angl.omega) * DasMath.CosD(angl.chi) * DasMath.CosD(angl.phi) -
                DasMath.SinD(angl.omega) * DasMath.SinD(angl.phi);
            uphi.j = DasMath.CosD(angl.omega) * DasMath.CosD(angl.chi) * DasMath.SinD(angl.phi) +
                DasMath.SinD(angl.omega) * DasMath.CosD(angl.phi);
            uphi.k = DasMath.CosD(angl.omega) * DasMath.SinD(angl.chi);

            return;
        }

        /// <summary>
        /// CALCULATE VECTOR COMPONENTS FROM ANGLES - 
        /// NOTE: OMEGA MUST BE DEVIATION FROM BISECTING POSITION DIMENSION ANGL(4),THPHI(3)
        /// Math: h = u x q as u is Equation(22) in Busing 
        /// </summary>
        /// <source>UBCALS.FOR -> CALCH</source>
        /// <param name="angl"></param>
        /// <param name="mWaveLength"></param>
        /// <param name="h_iphi"></param>
        /// <returns></returns>
        static public void VectorFromAngles(MotorIncidentAngles angl, double mWaveLength, Vector h_iphi)
        {
            Vector u_phi = new Vector(3);
            double sint = Math.Sin(angl.twotheta * 0.5 * DasMath.rad); 
#if false
            // Old Method... Not Used Anymore
            // 1. CALCULATE SIN THETA
            double sint = Math.Sin(angl.twotheta * 0.5 * DasMath.rad);

            double sino, sinx, sinp, coso, cosx, cosp;

            sino = Math.Sin(angl.omega * DasMath.rad);
            sinx = Math.Sin(angl.chi * DasMath.rad);
            sinp = Math.Sin(angl.phi * DasMath.rad);
            coso = Math.Cos(angl.omega * DasMath.rad);
            cosx = Math.Cos(angl.chi * DasMath.rad);
            cosp = Math.Cos(angl.phi * DasMath.rad);          
            

            // 2. FORM TERMS OF H-PHI
            /*
            h_iphi.i = 2.0 * sint * (coso * cosx * cosp - sino * sinp) / mWaveLength;
            h_iphi.j = 2.0 * sint * (coso * cosx * sinp + sino * cosp) / mWaveLength;
            h_iphi.k = 2.0 * sint * (coso * sinx) / mWaveLength;
            */
                        
            u_phi.i = coso * cosx * cosp - sino * sinp;
            u_phi.j = coso * cosx * sinp + sino * cosp;
            u_phi.k = coso * sinx;

            for (int i = 0; i < 3; i++)
                h_iphi.V[i] = 2.0 * sint * u_phi.V[i] / mWaveLength;
#else
            // New: Simpler      
            UnitVectorFromAngle(angl, u_phi);
            u_phi.Times(2.0*sint/mWaveLength, h_iphi);   
#endif

            return;
        }



        #endregion
    }

    /// <summary>
    /// A reflection
    /// </summary>
    public class Reflection
    {

        #region private properties
        private MillerIndices m_HKL;
        private MotorIncidentAngles m_motorAngles;
        private int trustdegreehkl; 
        #endregion

        #region Constructors
        /// <summary>
        /// constructor
        /// </summary>
        public Reflection()
        {
            m_HKL = new MillerIndices();
            m_motorAngles = new MotorIncidentAngles();
            trustdegreehkl = 0;
        }
        #endregion

        #region public properties
        /// <summary>
        /// How trustful the (HKL) can be
        /// -1: does  not make any sense
        /// 0 : the hkl type is correct
        /// 1 : it is correct
        /// </summary>
        public int TrustDegreeMillerIndices
        {
            set { trustdegreehkl = value; }
            get
            {
                return trustdegreehkl;
            }
        }
        /// <summary>
        /// miller indices (hkl) of the reflection 
        /// </summary>
        public MillerIndices MillerIndex
        {
            get { return m_HKL; }
            set { m_HKL = value; }
        }
        /// <summary>
        /// motor incident angle
        /// </summary>
        public MotorIncidentAngles MotorAngles
        {
            get { return m_motorAngles; }
            set { m_motorAngles = value; }
        }
        #endregion

        #region public methods
        /// <summary>
        /// Neat Print
        /// </summary>
        public void Print()
        {
            string ostr = "";
            ostr += "2theta = " + MotorAngles.twotheta.ToString("00.00") + "  ";
            ostr += "omega  = " + MotorAngles.omega.ToString("00.00") + "  ";
            ostr += "chi    = " + MotorAngles.chi.ToString("00.00") + "  ";
            ostr += "phi    = " + MotorAngles.phi.ToString("00.00") + "  ";
            ostr += String.Format("(H K L) = {0}, {1}, {2}", MillerIndex.H, MillerIndex.K, MillerIndex.L);
            Console.WriteLine(ostr);
            return;
            
        }

        /// <summary>
        /// Override ToSring()
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            String rbuf = String.Format("{0}   HKL: {1}", MotorAngles, m_HKL);
            return rbuf;
        }
        #endregion

        #region static public method

        /// <summary>
        /// Import a data file and convert to a 2-D array of strings
        /// </summary>
        /// <param name="filename">File name</param>
        /// <returns></returns>
        static public List<string[]> ImportDataFile(string filename)
        {
            // FIXME - This is a prototype
            StreamReader r = new StreamReader(filename);
            List<string[]> rlist = new List<string[]>();
            string line;

            while ((line = r.ReadLine()) != null)
            {
                string[] terms = line.Split();
                rlist.Add(terms);
#if false
                for (int i = 0; i < terms.Length; i++)
                {
                    Console.Write("{0},  ", terms[i]);
                }
                Console.WriteLine();
#endif
            }
            r.Close();

            return rlist;

        }

        /// <summary>
        /// Import reflections from a data file conformed to DAS format,i.e.,
        /// in the order of
        /// 2theta, omega, chi, phi, h, k, l for HB3A
        /// </summary>
        /// <param name="filename">Input file name</param>
        /// <returns></returns>
        static public List<Reflection> ImportFromFileDAS(string filename)
        {
            List<string[]> refstringlist = ImportDataFile(filename);
            List<Reflection> reflist = new List<Reflection>();

            Console.WriteLine("The data format is for HB3A Only.  Omega is NOT same as Busing's");

            for (int i = 0; i < refstringlist.Count; i++)
            {
                if (refstringlist[i].Length >= 5)
                {
                    // 1. Read from list
                    double twotheta = Convert.ToDouble(refstringlist[i][0]);
                    double omega = Convert.ToDouble(refstringlist[i][1]);
                    double chi = Convert.ToDouble(refstringlist[i][2]);
                    double phi = Convert.ToDouble(refstringlist[i][3]);
                    int hkl = Convert.ToInt16(refstringlist[i][4]);
                    int h = hkl / 100;
                    int k = (hkl - h * 100) / 10;
                    int l = hkl - h * 100 - k * 10;

                    // 2. Convert omega
                    omega = omega - twotheta / 2.0;

                    // 3. Put to list
                    Reflection newrefl = new Reflection();
                    newrefl.MotorAngles.twotheta = twotheta;
                    newrefl.MotorAngles.omega = omega;
                    newrefl.MotorAngles.chi = chi;
                    newrefl.MotorAngles.phi = phi;
                    newrefl.MillerIndex.H = h;
                    newrefl.MillerIndex.K = k;
                    newrefl.MillerIndex.L = l;

                    reflist.Add(newrefl);
                } // Ignore incorrect line
            }

            return reflist;

        }

        #endregion

    }

}