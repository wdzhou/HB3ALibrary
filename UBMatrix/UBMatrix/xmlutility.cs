using System;
using System.Collections.Generic;
using System.Text;

namespace DAS.UBMatrix
{
    /// <summary>
    /// Utility functions 
    /// </summary>
    public class XmlUtility
    {
        /// <summary>
        /// Convert a 2D array (matrix) to string of all elements 
        /// </summary>
        /// <param name="matrix"></param>
        /// <returns></returns>
        static public string ConvertToString(double[][] matrix)
        {
            string s = "";
            for (int i = 0; i < matrix.Length; i++)
                for (int j = 0; j < matrix[i].Length; j++)
                    s += String.Format("{0} ", matrix[i][j]);

            return s;
        }

        /// <summary>
        /// Convert data in the string to a 2D double array
        /// </summary>
        /// <param name="ins"></param>
        /// <param name="matrix"></param>
        static public void ConvertTo2DDoubleArray(string ins, double[][] matrix)
        {
            // 1. Parse and clean input string
            string[] rterms = ins.Split();
            List<string> terms = new List<string>();
            for (int i = 0; i < rterms.Length; i++)
                if (rterms[i] != "")
                    terms.Add(rterms[i]);

            int matrixsize = matrix.Length * matrix[0].Length;
            if (terms.Count != matrixsize)
            {
                string errmsg = string.Format("Input {0} Contains Different Number of Elements {2} Than Supposed To Be {1}",
                    ins, matrixsize, terms.Count);
                throw new Exception(errmsg);
            }
        
            for (int i = 0; i < matrix.Length; i++)
                for (int j = 0; j < matrix[i].Length; j++)
                    matrix[i][j] = Convert.ToDouble(terms[i * matrix.Length + j]);

            return;
        }


    }
}
