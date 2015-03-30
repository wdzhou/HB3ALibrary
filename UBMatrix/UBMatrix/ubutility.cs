using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace DAS.UBMatrix
{
    /// <summary>
    /// Class to store UB matrix candidates for future processing
    /// </summary>
    public class UBStorage
    {
        #region private property
        MotorIncidentAngles[] mRawIncidentAngles;
        List<MillerIndices[]> mCalHKL;
        List<double[]> mDiff2Thetas;
        List<double[]> mDiffUphis;
        List<Matrix> mUBMatrices;
        List<int[]> mRootIncidentAngleIndices;
        double tolerance2theta;
        double toleranceuphi;
        /// <summary>
        /// UB matrix average from all added UB matrix
        /// </summary>
        Matrix avgUB;
        /// <summary>
        /// sigma (uncerntainty of all UB)
        /// </summary>
        Matrix sigUB;
        #endregion

        #region constructors
        /// <summary>
        /// Initialization
        /// </summary>
        /// <param name="inprawincangles">list of incident angles as input</param>
        public UBStorage(MotorIncidentAngles[] inprawincangles)
        {
            mRawIncidentAngles = new MotorIncidentAngles[inprawincangles.Length];
            for (int i = 0; i < inprawincangles.Length; i++)
            {
                mRawIncidentAngles[i] = new MotorIncidentAngles();
                inprawincangles[i].CopyTo(mRawIncidentAngles[i]);
            }

            mCalHKL = new List<MillerIndices[]>();
            mDiff2Thetas = new List<double[]>();
            mDiffUphis = new List<double[]>();
            mUBMatrices = new List<Matrix>();
            mRootIncidentAngleIndices = new List<int[]>();

            avgUB = new Matrix(3);
            sigUB = new Matrix(3);

            tolerance2theta = 1.0;
            toleranceuphi = 1.0;

            return;
        }
        #endregion

        #region public property
        /// <summary>
        /// tolerance on calculated and observed 2theta
        /// </summary>
        public double Tolerance2Theta{
            set { tolerance2theta = value; }
            get { return tolerance2theta; }
        }
        /// <summary>
        /// tolerance on the angle b/w calculated and observed u_phi
        /// </summary>
        public double ToleranceUphi
        {
            get { return toleranceuphi; }
            set { toleranceuphi = value; }
        }
#endregion

        #region public methods
        /// <summary>
        /// Add a reflection result to the storage
        /// </summary>
        /// <param name="calUB"></param>
        /// <param name="anglepair"></param>
        /// <param name="inpcalhkl"></param>
        /// <param name="inpdiff2theta"></param>
        /// <param name="inpdiffUphi"></param>
        public void AddUBSuite(Matrix calUB, int[] anglepair, MillerIndices[] inpcalhkl, 
            double[] inpdiff2theta, double[] inpdiffUphi)
        {
            Matrix newub = new Matrix(3);
            calUB.CopyToMatrix(newub);

            int[] newangleindexpair = new int[2];
            for (int i = 0; i < 2; i++)
                newangleindexpair[i] = anglepair[i];
     
            MillerIndices[] newcalhkl = new MillerIndices[mRawIncidentAngles.Length];
            double[] diff2theta = new double[mRawIncidentAngles.Length];
            double[] diffUphi = new double[mRawIncidentAngles.Length];
            for (int i = 0; i < inpcalhkl.Length; i++)
            {
                newcalhkl[i] = new MillerIndices();
                inpcalhkl[i].CopyTo(newcalhkl[i]);

                diff2theta[i] = inpdiff2theta[i];
                diffUphi[i] = inpdiffUphi[i];
            }

            mUBMatrices.Add(newub);
            mRootIncidentAngleIndices.Add(newangleindexpair);
            mCalHKL.Add(newcalhkl);
            mDiff2Thetas.Add(diff2theta);
            mDiffUphis.Add(diffUphi);

            return;
        }

        /// <summary>
        /// Do statistic to all the stored UB matrix
        /// </summary>
        public void Stat()
        {
            for (int i = 0; i < 3; i ++)
                for (int j = 0; j < 3; j++)
                {
                    double d1 = 0.0;
                    double d2 = 0.0;
                    //Console.WriteLine("{0}  {1}", i, j);
                    for (int k = 0; k < mUBMatrices.Count; k++)
                    {
                        //Console.WriteLine("k = {0}", k);
                        //Console.WriteLine("{0}", mUBMatrices[k].ToString());
                        d1 += mUBMatrices[k].M[i][j];
                        d2 += (mUBMatrices[k].M[i][j]) * (mUBMatrices[k].M[i][j]);
                        // Console.WriteLine("d1 = {0}    d2 = {1}", d1, d2);
                    }
                    // Console.WriteLine("d1 = {0} summed over {1} UB matrices", d1, mUBMatrices.Count);
                    avgUB.M[i][j] = d1 / mUBMatrices.Count;
                    sigUB.M[i][j] = Math.Sqrt(d2 / mUBMatrices.Count - avgUB.M[i][j] * avgUB.M[i][j]);
                }
            //avgUB.Print();
            //sigUB.Print();
            return;
        }

        /// <summary>
        /// Print customized format string
        /// </summary>
        /// <param name="filename"></param>
        public void Print(string filename)
        {
            string wbuf = "";

            // 1. Print all the information
            for (int i = 0; i < mUBMatrices.Count; i++)
            {
                wbuf += String.Format("Result {0}\n", i);
                wbuf += String.Format("\tFrom Incident Angle {0}, {1}\n",
                    mRootIncidentAngleIndices[i][0].ToString(), mRootIncidentAngleIndices[i][1].ToString());
                wbuf += String.Format("\tList of Calculate HKL:\n");
                for (int j = 0; j < mRawIncidentAngles.Length; j++)
                    wbuf += String.Format("\t\t{0}\t\tdelta(2theta) = {1},   delta(u_phi) = {2}\n",
                        mCalHKL[i][j].ToString(), mDiff2Thetas[i][j], mDiffUphis[i][j]);
                for (int j = 0; j < 3; j++)
                {
                    wbuf += "\t";
                    for (int k = 0; k < 3; k++)
                        wbuf += String.Format("{0}   ", mUBMatrices[i].M[j][k]);
                    wbuf += "\n";
                }
                wbuf += "\n\n";                                
            }

            // 2. Statistic information
            wbuf += "UB Matrix Statistic:\n";
            for (int i = 0; i < 3; i++)
            {
                wbuf += "\t";
                for (int j = 0; j < 3; j++)
                    wbuf += String.Format("{0}   ", avgUB.M[i][j]);
                wbuf += "\n";
            }
            for (int i = 0; i < 3; i++)
            {
                wbuf += "\t";
                for (int j = 0; j < 3; j++)
                    wbuf += String.Format("{0}   ", sigUB.M[i][j]);
                wbuf += "\n";
            }
                


            // 3. Print to file
            if (filename.Length < 3)
            {
                Console.WriteLine(wbuf);
            }
            else
            {

                using (StreamWriter sw = File.CreateText(filename))
                {
                    sw.WriteLine(wbuf);
                    sw.Close();
                }
            }
            return;

        }
        #endregion

    }
}
