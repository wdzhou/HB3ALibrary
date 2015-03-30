using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace DAS.UBMatrix
{
    public class StringUtility
    {
        /// <summary>
        /// Parse string to integers
        /// </summary>
        /// <param name="input"></param>
        /// <param name="numintegers"></param>
        /// <param name="parsedvalues"></param>
        static public void ParseInteger(string input, int numintegers, int[] parsedvalues)
        {
            string[] terms = input.Split(new char[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
            if (terms.Length < numintegers || numintegers > parsedvalues.Length)
            {
                string errormessage = String.Format(
                    "ParseInteger(): Input Error!  {0} has no enough items comparing to {1} and input array length {2}or input array length {2} is smaller than {1}",
                    input, numintegers, parsedvalues.Length);
                throw new Exception(errormessage);
            }

            try
            {
                for (int i = 0; i < numintegers; i++)
                {
                    parsedvalues[i] = Convert.ToInt32(terms[i]);
                }

            }
            catch (FormatException fe)
            {
                string errormessage = String.Format("Invalid Item To Convert To Integer In Input {0}\n{1}",
                    input, fe.Message);
                throw new FormatException(errormessage);
            }

        }

        /// <summary>
        /// Parse string to doubles
        /// </summary>
        /// <param name="input"></param>
        /// <param name="numdoubles"></param>
        /// <param name="parsedvalues"></param>
        static public void ParseDoubles(string input, int numdoubles, double[] parsedvalues)
        {
            string[] terms = input.Split(new char[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
            if (terms.Length < numdoubles || numdoubles > parsedvalues.Length)
            {
                string errormessage = String.Format(
                    "ParseInteger(): Input Error!  {0} has no enough items comparing to {1} and input array length {2}or input array length {2} is smaller than {1}",
                    input, numdoubles, parsedvalues.Length);
                throw new Exception(errormessage);
            }

            try
            {
                for (int i = 0; i < numdoubles; i++)
                {
                    parsedvalues[i] = Convert.ToDouble(terms[i]);
                }

            }
            catch (FormatException fe)
            {
                string errormessage = String.Format("Invalid Item To Convert To Integer In Input {0}\n{1}",
                    input, fe.Message);
                throw new FormatException(errormessage);
            }

        }

        /// <summary>
        /// Parse string to integers
        /// </summary>
        /// <param name="input"></param>
        /// <param name="numintegers"></param>
        /// <param name="parsedvalues"></param>
        static public int ParseIndefinateInteger(string input, int[] parsedvalues)
        {
            string[] terms = input.Split(new char[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);

            int numvalues = 0;

            try
            {
                for (int i = 0; i < terms.Length; i++)
                {
                    try
                    {
                        parsedvalues[numvalues] = Convert.ToInt32(terms[i]);
                        numvalues++;
                    }
                    catch (FormatException fe)
                    {
                        string errormessage = String.Format("Invalid Item To Convert To Integer In Input {0}\n{1}",
                            input, fe.Message);
                        throw new FormatException(errormessage);
                    }
                }
            }
            catch (IndexOutOfRangeException ioore)
            {
                string errormessage = String.Format(
                    "ParseInteger(): Input Error {3}!  {0} has no enough items comparing to {1} and input array length {2}or input array length {2} is smaller than {1}",
                    input, numvalues, parsedvalues.Length, ioore.Message);
                throw new Exception(errormessage);

            }

            return numvalues;
        }
    }

    public class UBMatrixUserInterface
    {
        /// <summary>
        /// Parse cell parameters file
        /// Requirement: all in one line in sequence of a, b, c, alpha, beta, gamma
        /// The data must be in the first line
        /// </summary>
        /// <param name="unitcellfilename"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="alpha"></param>
        /// <param name="beta"></param>
        /// <param name="gamma"></param>
        static public void ReadCellFile(string unitcellfilename,
            out double a, out double b, out double c,
            out double alpha, out double beta, out double gamma)
        {
            string cellline = "";
            try
            {
                StreamReader sr = new StreamReader(unitcellfilename);
                cellline = sr.ReadLine();
                if (cellline[0] == '#')
                    cellline = sr.ReadLine();
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Unable to read file %s due to %s.", unitcellfilename, e.Message);
                throw new Exception("UB");
            }

            double[] doubles = new double[6];
            StringUtility.ParseDoubles(cellline, 6, doubles);

            a = doubles[0];
            b = doubles[1];
            c = doubles[2];
            alpha = doubles[3];
            beta = doubles[4];
            gamma = doubles[5];

            return;
        }

    }
}
