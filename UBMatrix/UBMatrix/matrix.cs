using System;
using System.Collections.Generic;
using System.Text;

namespace DAS.UBMatrix
{
    /// <summary>
    /// General matrix, including all operation with it
    /// </summary>
    public class Matrix
    {
        const double TINY = 1.0e-20;

        #region Protected properties
        /// <summary>
        /// matrix size for square matrix
        /// </summary>
        protected int mSize;
        /// <summary>
        /// row size (a for matrix a x b)
        /// </summary>
        protected int rowSize;
        /// <summary>
        /// column size (b for matrix a x b)
        /// </summary>
        protected int colSize;
        /// <summary>
        /// matrix storage 2D array
        /// </summary>
        protected double[][] mM;
        #endregion

        #region Constructors
        /// <summary>
        /// Constructor for square matrix
        /// </summary>
        /// <param name="size"></param>
        public Matrix(int size)
        {
            mSize = size;
            mM = new double[mSize][];
            for (int i = 0; i < mSize; i++)
            {
                mM[i] = new double[mSize];
            }
            rowSize = mSize;
            colSize = mSize;
        }

        /// <summary>
        /// General matrix constructor 
        /// </summary>
        /// <param name="inprowsize"></param>
        /// <param name="inpcolsize"></param>
        public Matrix(int inprowsize, int inpcolsize)
        {
            // 1. determine size
            rowSize = inprowsize;
            colSize = inpcolsize;
            if (rowSize != colSize)
                mSize = -1;
            else
                mSize = rowSize;

            // 2. 2D array
            mM = new double[rowSize][];
            for (int i = 0; i < rowSize; i++)
            {
                mM[i] = new double[colSize];
            }

            return;
        }
        #endregion

        #region Public properites
        /// <summary>
        /// row size
        /// </summary>
        public int RowSize
        {
            get { return rowSize; }
        }
        /// <summary>
        /// column size
        /// </summary>
        public int ColSize
        {
            get { return colSize; }
        }
        /// <summary>
        /// storage 2d array
        /// </summary>
        public double[][] M
        {
            get { return mM; }
        }
        #endregion

        #region Public static methods
        /// <summary>
        /// C = A X B
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="C"></param>
        static public void MulMatrix(Matrix A, Matrix B, Matrix C)
        {
            // 1. Check for the matrix size
            if (A.rowSize != B.colSize || A.colSize != B.rowSize || C.mSize != A.rowSize)
            {
                Console.WriteLine("Cannot multipy matrix with difference size");
                throw new System.Exception();
            }

            // 2. Do multiply
            for (int i = 0; i < A.rowSize; i++)
            {
                for (int j = 0; j < B.colSize; j++)
                {
                    C.mM[j][i] = 0.0;
                    for (int k = 0; k < A.colSize; k++)
                    {
                        C.mM[j][i] = C.mM[j][i] + A.mM[j][k] * B.mM[k][i];
                    }
                }
            }

            return;
        }

        /// <summary>
        /// LU Decomposition
        /// </summary>
        /// <note>Modified according to numerical recipe</note>
        /// <LastMajorModification>2009.09.30</LastMajorModification>
        /// <param name="A">Matrix In = A, Out = LU</param>
        /// <param name="indx">Out = row permutation</param>
        /// <param name="d">Out = flag of even or odd number of row interchanges</param>
        /// <returns>Error code</returns>
        static public int LUDcomp(Matrix A, int[] indx,out double d)
        {         	        
            //nvar in the single equivilent is = to np here.
	        double sum, aamax, dum;
            int np = A.mSize;
            double[] vv = new double[np];

	        d = 1.0;

            // Loop over rows to get the implicit scaling information
            for (int i = 0; i < np; i++)
            {
                aamax = 0.0;
                for (int j = 0; j < np; j++)
                {
                    if (Math.Abs(A.mM[i][j]) > aamax)
                        aamax = Math.Abs(A.mM[i][j]);
                }

                // all elements in Matrix A is zero
                if (aamax == 0.0)
                {
                    Console.WriteLine("Singular Matrix.  Error Code -2");
                    return -2; //singlar matrix....
                }

                vv[i] = 1.0 / aamax;
            } // i 

            // Loop over columns of Crout's method 
	        for (int j = 0; j < np; j++)
	        {
                // Equation 2.3.12
                for (int i = 0; i < j; i++)
                {
                    sum = A.mM[i][j];
                    for (int k = 0; k < i; k++)
                    {
                        sum = sum - A.mM[i][k] * A.mM[k][j];
                    }
                    A.mM[i][j] = sum;

                }  // i 

                // Init to Search for largest pivoting
		        aamax=0.0;
                int imax = -1;
                // i = j case of (2.3.12) and i > j case of (2.3.13)
                for (int i = j; i < np; i++)
                {
                    sum = A.mM[i][j];

                    for (int k = 0; k < j; k++)
                    {
                        sum = sum - A.mM[i][k] * A.mM[k][j];
                    }
                    A.mM[i][j] = sum;

                    dum = vv[i] * Math.Abs(sum);

                    // If the figure of merit for the pivot better than the best so far?
                    if (dum >= aamax)
                    {
                        imax = i;
                        aamax = dum;
                    }
                } // i 

                // Do we need to interchange rows?
		        if (j != imax)
                {
			        for (int k = 0; k < np; k++)
                    {
				        dum=A.mM[imax][k];
				        A.mM[imax][k] = A.mM[j][k];
				        A.mM[j][k]=dum;
			        }
			        d=-d;
			        vv[imax]=vv[j];
                } // if (j)

                // 
		        indx[j] = imax;
                if (A.mM[j][j] == 0.0)
                {
                    A.mM[j][j] = TINY;
                }

                // Divide by pivot element
		        if (j != np-1)
                {
			        dum=1.0/A.mM[j][j];
			        for (int i= j+1; i < np; i++)
			        {
				        A.mM[i][j]=A.mM[i][j]*dum;
                    }
                } // if j
            } // for j

	        return 0;
        }

        /// <summary>
        /// 2nd step of LU decomposition
        /// </summary>
        /// <param name="A">LU Matrix</param>
        /// <param name="indx">pivoting index</param>
        /// <param name="b">in = right-hand side vector.  out = solution</param>
        /// <returns></returns>
        static public int LUBkSb(Matrix A, int[] indx, double[] b)
        {  //n is always np in this program.....

	        int ii = -1,ll;	       
	        double sum;
            int np = A.mSize;
	
            // forward substituting 
            for (int i=0; i < np; i++)
	        {
                // as ii is a non-negative value, it will become 
                // the index of the first nonvanishing element of b
		        ll = indx[i];
		        sum = b[ll];		
                b[ll]=b[i];

		        if (ii >= 0)
                {
			        for (int j = ii; j <= i-1; j++)
                    {
				        sum = sum - A.mM[i][j] * b[j];
                    } // j
                }
		        else if (sum > 0.0)
                {
			        ii=i;
                }                   
                    
                b[i]=sum;
            } // i

            // backward substituting
	        for (int i = np-1; i >= 0; i--)
            {
		        sum = b[i];
		        if (i < np)
                {
			        for (int j = i+1; j < np; j++)
                    {
				        sum=sum - A.mM[i][j] * b[j];
                    }
                }
		        b[i]=sum / A.mM[i][i];
            }

	        return 0;
        }

        /// <summary>
        /// Generate an orthogonal matrix with 2 given vector
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <param name="X"></param>
        /// <returns></returns>               
        static public int Orthog(Vector v1, Vector v2, Matrix X)
        {

            Vector tempv, uv1, uv2;

            tempv = new Vector(v1.dimension);

            uv1 = v1.UnitVector();
            X.mM[0][0] = uv1.i;
            X.mM[1][0] = uv1.j;
            X.mM[2][0] = uv1.k;

            Vector.CrossProduct(v1, v2, tempv);

            uv2 = tempv.UnitVector();
            X.mM[0][2] = uv2.i;
            X.mM[1][2] = uv2.j;
            X.mM[2][2] = uv2.k;

            Vector.CrossProduct(uv2, uv1, tempv);

            X.mM[0][1] = tempv.i;
            X.mM[1][1] = tempv.j;
            X.mM[2][1] = tempv.k;

            return 0;
        }

        /// This MakeBasis is different from Matrix.Orthog()
        /// FIXME - Should put this method back to Matrix
        ///<source>SUBROUTINE MAKEBASIS(V1,V2,D)</source>
        ///<algorithm>
        ///3. the resultant three basis are at
        ///   a. u_i = u1 + u2 / |u1 + u2|;
        ///   b. u_k = v1 x v2 ;
        ///   c. u_j = u_k x u_i;
        ///</algorithm>
        static public void MakeBasis(Vector V1, Vector V2, Matrix D)
        {
            Vector[] vUs = new Vector[3];
            Vector vU1, vU2;
            Vector vUI = new Vector(3);                        
            Vector vUJ = new Vector(3);
            Vector vUK = new Vector(3);

            // 1. V1 X V2 / |V1 X V2|
            vU1 = V1.UnitVector();
            vU2 = V2.UnitVector();
            //Console.WriteLine("Flag Pole.... {0}, {1}, {2}", vU1.dimension, vU2.dimension, vUI.dimension);
            Vector.Plus(vU1, vU2, vUI);
            //Console.WriteLine("DB351-1  vU1 = {0}", vU1.ToString());
            //Console.WriteLine("DB351-2  vU2 = {0}", vU2.ToString());
            //Console.WriteLine("DB351-3  vUI = {0}", vUI.ToString());
            vUs[0] = vUI.UnitVector();

            // 2. u_k = v1 x v2
            Vector.CrossProduct(V1, V2, vUK);
            vUs[2] = vUK.UnitVector();

            // 3. u_j = v_k x v_i
            vUs[1] = new Vector(3);
            Vector.CrossProduct(vUs[2], vUs[0], vUs[1]);

#if false
            for (int i = 0; i < 3; i++)
                Console.WriteLine("DB319  v_{0}:  {1}", i, vUs[i]);
#endif

            // 4. Construct matrix
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    D.M[j][i] = vUs[i].V[j];
                }
            }
            return;
        }

        /// <summary>
        /// Matrix multiply Vector c = M x b
        /// </summary>
        /// <param name="M"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        static public void MultVector(Matrix M, Vector b, Vector c)
        {
            // 1. Check
            if (M.mSize != b.dimension || M.mSize != c.dimension)
            {
                throw new Exception("Matrix and Vecotr(s) with different size");
            }

            // 2. Calculate
            for (int i = 0; i < b.dimension; i++)
            {
                double subsum = 0.0;
                for (int j = 0; j < b.dimension; j++)
                {
                    subsum += M.mM[i][j] * b.Get(j);
                }
                c.Set(i, subsum);
            }

            return;
        }

        #endregion 

        #region Public methods

        /// <summary>
        /// Customerized information output
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string ms = "";
            // Console.WriteLine("size = {0}", mSize);
            for (int i = 0; i < mSize; i++)
            {
                for (int j = 0; j < mSize; j++)
                    ms += String.Format("{0,7:0.0000}   ", mM[i][j]);
                ms += "\n";
            }
            return ms;
        }

        /// <summary>
        /// Copy all the elements to another matrix
        /// </summary>
        /// <param name="Y"></param>
        public void CopyToMatrix(Matrix Y)
        {
            if (mSize != Y.mSize){
                throw new Exception("Size No Match");
            }

            for (int i = 0; i < mSize; i++)
                for (int j = 0; j < mSize; j++)
                    Y.M[i][j] = M[i][j];
            return;
        }

        /// <summary>
        /// Invert this matrix
        /// </summary>
        /// <param name="Y">Result matrix (output)</param>
        /// <param name="det">Output determination</param>
        /// <returns></returns>
        public int InvertMatrix(Matrix Y, out double det)
        {
            if (mSize <= 0)
                throw new Exception("Non-Square Matrix Inversion Has NOT Been Tested");

            int np = mSize;

            // 1. LU decomposition
            double d;
            int[] index = new int[np];

            Matrix lumatrix = this.Duplicate();
            int errorcode = LUDcomp(lumatrix, index, out d);
            if (errorcode == 1){
                // Exit case: singular matrix!
                det = 0.0;
                return -1;
            } else {
                det = 1.0;
                for (int i = 0; i < np; i++)
		        {
			        det = det * mM[i][i];
                }
            }

            // 3. Find inverse by columns (back substitution)
            double[] col = new double[np];
            for (int j = 0; j < np; j++)
	        {
                // for a unit vector at j-direction
                for (int i = 0; i < np; i++)
                {
                    col[i] = 0.0;
                }
                col[j] = 1.0;

                // back substituion
                LUBkSb(lumatrix, index, col);

                // put result to matrix to return
                Y.SetColumn(np, j, col);

            }
	 
            return 0;
        }

        /// <summary>
        /// Transpose 
        /// </summary>
        /// <param name="outmatrix"></param>
        /// <returns></returns>
        public int Transpose(Matrix outmatrix)
        {
            int np = mSize;

            if (outmatrix.rowSize != colSize ||
                outmatrix.colSize != rowSize)
                throw new Exception("Transposed Matrix Dimension Error!");

	        for (int i = 0; i < rowSize; i++)
	        {
		        for (int j = 0; j < colSize; j++)
		        {
			        outmatrix.mM[j][i] = mM[i][j];
		        }
	        }

	        return 0;
        }

        /// <summary>
        /// Set a vector to a column in matrix
        /// </summary>
        /// <param name="col"></param>
        /// <param name="vector"></param>
        public void SetColumn(int col, Vector vector)
        {
            for (int i = 0; i < vector.dimension; i++)
            {
                mM[i][col] = vector.Get(i);
            }
            return;
        }

        /// <summary>
        /// Set a vector to a column in matrix
        /// </summary>
        /// <param name="n">size of vector</param>
        /// <param name="col">number of column to set value</param>
        /// <param name="vector">array of doubles with size n</param>
        public void SetColumn(int n, int col, double[] vector)
        {
            for (int i = 0; i < n; i++)
            {
                mM[i][col] = vector[i];
            }
            return;
        }

        /// <summary>
        /// Set all elements in the matrix to zero
        /// </summary>
        public void SetToZero()
        {
            for (int i = 0; i < mSize; i++)
                for (int j = 0; j < mSize; j++)
                    mM[i][j] = 0;
            return;
        }

        /// <summary>
        /// Project a row to a vector
        /// </summary>
        /// <param name="row"></param>
        /// <returns></returns>
        public Vector Project(int row)
        {
            Vector v = new Vector(mSize);
            for (int i = 0; i < mSize; i++)
            {
                v.Set(i, mM[row][i]);
            }
            return v;
        }

        /// <summary>
        /// Duplicate this matrix
        /// </summary>
        /// <returns></returns>
        public Matrix Duplicate()
        {
            Matrix nmatrix = new Matrix(mSize);
            for (int i = 0; i < mSize; i++)
            {
                for (int j = 0; j < mSize; j++)
                {
                    nmatrix.mM[i][j] = mM[i][j];
                }
            }

            return nmatrix;

        }

        /// <summary>
        /// Print itself
        /// </summary>
        public void Print()
        {
            string ms = ToString();
            Console.WriteLine("Matrix:\n{0}", ms);
        }

        /// <summary>
        /// Print Matrix
        /// </summary>
        /// <param name="longdigit">select b/w long digit output or short</param>
        public void Print(bool longdigit)
        {
            if (longdigit == true)
            {
                string ms = "";
                // Console.WriteLine("size = {0}", mSize);
                for (int i = 0; i < mSize; i++)
                {
                    for (int j = 0; j < mSize; j++)
                        ms += String.Format("{0,15:0.000000000}   ", mM[i][j]);
                    ms += "\n";
                }
                Console.WriteLine(ms);
            }
            else
                Print();
            return;
        }


        /// <summary>
        /// c = M x b
        /// </summary>
        /// <param name="b">(in)</param>
        /// <param name="c">(out)</param>
        public void MultVector(Vector b, Vector c)
        {
            // Check vector and matrix size matching or not
            if (mSize != b.dimension)
                throw new Exception("Matrix and Vecotr with different size");
      
            // Calculate          
            for (int i = 0; i < b.dimension; i++)
            {
                double subsum = 0.0;
                for (int j = 0; j < b.dimension; j++)
                {
                    subsum += mM[i][j] * b.Get(j);
                }
                c.Set(i, subsum);
            }

            return;
        }

        /// <summary>
        /// Matrix multiply Vector
        /// </summary>
        /// <param name="b">Vector</param>
        /// <returns>Vector</returns>
        public Vector MultVector(Vector b)
        {
            Vector c = new Vector(b.dimension);
            MultVector(b, c);
            return c;
        }

        #endregion

    }
}
