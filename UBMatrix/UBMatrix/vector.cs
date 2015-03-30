using System;
using System.Collections.Generic;
using System.Text;

namespace DAS.UBMatrix
{
    /// <summary>
    /// Vector
    /// </summary>
    public class Vector
    {
        #region protected properites
        /// <summary>
        /// 1D array for storing the vector
        /// </summary>
        protected double[] mVector;
        /// <summary>
        /// Vector size
        /// </summary>
        protected int mDimension; // in which space the vector located
        #endregion

        #region Constructors
        /// <summary>
        /// Constructor from (0, 0, 0)
        /// </summary>
        public Vector(int indimension)
        {
            mDimension = indimension;
            mVector = new double[mDimension];
            for (int i = 0; i < mDimension; i++)
            {
                mVector[i] = 0;
            }
            return;
        }

        /// <summary>
        /// Constructor with given values
        /// </summary>
        /// <param name="ii"></param>
        /// <param name="ij"></param>
        /// <param name="ik"></param>
        public Vector(double ii, double ij, double ik)
        {
            mVector = new double[3];
            mVector[0] = ii;
            mVector[1] = ij;
            mVector[2] = ik;
        }
        #endregion

        #region public properties
        /// <summary>
        /// scale at i (x) direction
        /// </summary>
        public double i
        {
            get{ return mVector[0]; }
            set { mVector[0] = value; }
        }

        /// <summary>
        /// scale at j (y) direction
        /// </summary>
        public double j
        {
            get { return mVector[1]; }
            set { mVector[1] = value; }
        }

        /// <summary>
        /// value at k (z) direction
        /// </summary>
        /// 
        public double k
        {
            get { return mVector[2]; }
            set { mVector[2] = value; }
        }
        /// <summary>
        /// size of the vector
        /// </summary>
        public int dimension
        {
            get { return mDimension; }
            set { mDimension = value; }
        }
        /// <summary>
        /// 1D array as stroage of the vector 
        /// </summary>
        public double[] V
        {
            get { return mVector; }
        }
        #endregion

        #region public methods

        /// <summary>
        /// Customerized output
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string rst = "";
            for (int i = 0; i < mVector.Length; i ++){
                rst += String.Format("{0,7:0.0000}, ", mVector[i]);
            }
            rst += String.Format("  |V| = {0}", Length());
            return rst;
        }

        /// <summary>
        /// Calculate the length/scale of the vector
        /// </summary>
        /// <returns></returns>
        public double Length()
        {
            double l2 = 0;
            for (int i = 0; i < 3; i ++){
                l2 += mVector[i] * mVector[i];
            }
            double l = Math.Sqrt(l2);
            return l;
        }

        /// <summary>
        /// Length of vector projected to a plane 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public double Project(int a, int b)
        {
            double l2 = mVector[a] * mVector[a] + mVector[b] * mVector[b];
            double l = Math.Sqrt(l2);
            return l; 
        }

        /// <summary>
        /// B = c X I
        /// </summary>
        /// <param name="c"></param>
        /// <param name="vB"></param>
        public void Times(double c, Vector vB)
        {
            for (int i = 0; i < vB.dimension; i++)
                vB.mVector[i] = c * mVector[i];
            return;
        }

        /// <summary>
        /// Normalize this vector
        /// </summary>
        /// <param name="vN"></param>
        public void Normalize(Vector vN)
        {
            for (int i = 0; i < 3; i++)
                vN.V[i] = this.V[i] / this.Length();
           
            return;
        }


        /// <summary>
        /// C = a x b
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <param name="w"></param>
        static public void CrossProduct(Vector u, Vector v, Vector w)
        {
            w.i = u.j * v.k - v.j * u.k;
            w.j = u.k * v.i - v.k * u.i;
            w.k = u.i * v.j - v.i * u.j;

            return;
        }

        /// <summary>
        /// prod = u .dot v
        /// </summary>
        /// <param name="u">Vector</param>
        /// <param name="v">Vector</param>
        static public double DotProduct(Vector u, Vector v)
        {
            double prod = 0;
            for (int i = 0; i < 3; i++)
            {
                prod += u.mVector[i] * v.mVector[i];
            }

            return prod;
        }

        /// <summary>
        /// V_c = V_a + V_b
        /// </summary>
        /// <param name="vA"></param>
        /// <param name="vB"></param>
        /// <param name="vC"></param>
        static public void Plus(Vector vA, Vector vB, Vector vC)
        {
            // Console.WriteLine("Length = {0}", vC.Length());
            for (int i = 0; i < vC.mDimension; i++)
            {
                vC.V[i] = vA.V[i] + vB.V[i];
                // Console.WriteLine("vC[{0}] = {1}", i, vC.V[i]);
            }
            // Console.Write("vC = {0}", vC.ToString());
            return;
        }

        /// <summary>
        /// B = c X A
        /// </summary>
        /// <param name="vA"></param>
        /// <param name="c"></param>
        /// <param name="vB"></param>
        static public void Times(Vector vA, double c, Vector vB)
        {
            for (int i = 0; i < vA.Length(); i++)
                vB.V[i] = c * vA.V[i];
            return;
        }

        /// <summary>
        /// Set a value to a specific position in vector
        /// </summary>
        /// <param name="pos"></param>
        /// <param name="v"></param>
        public void Set(int pos, double v)
        {
            mVector[pos] = v;
            return;
        }

        /// <summary>
        /// Set To Zero
        /// </summary>
        public void SetToZero()
        {
            for (int i = 0; i < mVector.Length; i++)
                mVector[i] = 0.0;
            return;
        }

        /// <summary>
        /// Get a value from a specific position
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public double Get(int index)
        {
            return mVector[index];
        }

        /// <summary>
        /// Copy the content to another Vector
        /// </summary>
        /// <param name="target"></param>
        public void CopyTo(Vector target)
        {
            for (int i = 0; i < mVector.Length; i++)
                target.mVector[i] = mVector[i];
            return;
        }

        /// <summary>
        /// Generate a unit vector
        /// </summary>
        /// <returns></returns>
        public Vector UnitVector()
        {
            Vector outV = new Vector(mDimension);

            UnitVector(outV);

            return outV;

        }

        /// <summary>
        /// Calculate unit vector
        /// </summary>
        /// <param name="vUnit">Vector</param>
        public void UnitVector(Vector vUnit)
        {
            double sum = 0;
            for (int i = 0; i < mDimension; i++)
            {
                sum += mVector[i] * mVector[i];
            }
            double length = Math.Sqrt(sum);

            if (length < 1.0E-6)
            {
                length = 1.0E-6;
            }

            for (int i = 0; i < mDimension; i++)
            {
                vUnit.mVector[i] = mVector[i] / length;
            }

            return;
        }

        /// <summary>
        /// Print the vector
        /// </summary>
        public void Print()
        {
            string outputstring = "";
            outputstring += String.Format("Vector in {0}-D:  ", mDimension);
            // outputstring += "Vector in " + mDimensiion.ToString() + "-D";
            for (int i = 0; i < mDimension; i++)
            {
                outputstring += String.Format("{0}  ,", mVector[i]);
            }
            Console.WriteLine(outputstring);
        }
    	
        #endregion

    }
}
