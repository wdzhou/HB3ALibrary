using System;
using System.Collections.Generic;
using System.Text;

using DAS.UBMatrix;

namespace Testers
{
    /// <summary>
    /// System test is used to verify all the major functionalities provided by LVAdaptor for UB matrix renders
    /// correct result
    /// </summary>
    public class DebugTest
    {
        static public void Main() // _DebugTest
        {
            Console.WriteLine("********************************************************");
            Console.WriteLine("* Debug test for UB matrix and its LabVIEW adaptor. ");
            Console.WriteLine("********************************************************");

            // Initialization of test values
            //====================================================================================
            const double lambda = 1.003;
            const string latticeparameterstring = "5.4  7.6 5.3058 90 90 90";
            const string hkllistfname = "HB3A_Sample02_Reflections.txt";
            const int maxlength = 100;

            double delta_2theta = 0.5;
            double step_2theta = 0.05;
            double delta_omega = 0.5;
            double step_omega = 0.05;
            double delta_chi = 0.0;
            double step_chi = 0.05;
            //====================================================================================

            // Create instance
            LabVIEWAdaptor ubcalculator = new LabVIEWAdaptor(lambda, maxlength);

            // Set up lattice
            ubcalculator.SetLatticeParameters(latticeparameterstring);

            // Import reflections
            ubcalculator.ImportReflectionsFromFile(hkllistfname);

            //====================================================================================
            // Calculate UB Matrix  
            double[] motorangles1 = new double[4];
            int[] hkl1 = new int[3];
            ubcalculator.GetImportedReflection(0, motorangles1, hkl1);

            double[] motorangles2 = new double[4];
            int[] hkl2 = new int[3];
            ubcalculator.GetImportedReflection(1, motorangles2, hkl2);

            ubcalculator.GenerateUBMatrix2Reflections(hkl1, motorangles1, hkl2, motorangles2);

            // Calcualte (HKL) by calculated miller index
            /*
            int[] hklout1 = new int[3];
            double[] hkloutd1 = new double[3];
            ubcalculator.CalculateMillerIndexFromMotorAngle(hklout1, hkloutd1, motorangles1);

            string wbuf1 = "Reflection 1 (HKL-exp,  HKL-cal,  HKL-round)\n";
            for (int i = 0; i < 3; ++i)
            {
                wbuf1 += String.Format("{0}\t{1}\t{2}\n", hkl1[i], hkloutd1[i], hklout1[i]);
            }

            int[] hklout2 = new int[3];
            double[] hkloutd2 = new double[3];
            ubcalculator.CalculateMillerIndexFromMotorAngle(hklout2, hkloutd2, motorangles2);

            string wbuf2 = "Reflection 2 (HKL-exp,  HKL-cal,  HKL-round)\n";
            for (int i = 0; i < 3; ++i)
            {
                wbuf2 += String.Format("{0}\t{1}\t{2}\n", hkl2[i], hkloutd2[i], hklout2[i]);
            }
           
            if (!(ba0 && ba1 && ba2 && ba3 && ba4 && ba5))
            {
                Console.WriteLine("*********************************************");
                Console.WriteLine("* System Test Error In Step 2.  Test Stopped!");
                Console.WriteLine("*********************************************");
                return;
            }

            // Calculate angles between 2 reflections
            double anglehkl = ubcalculator.CalculateAngleBW2MillerIndices(hkl1, hkl2);
            double anglemotor = ubcalculator.CalculateAngleBW2MotorAngles(motorangles1, motorangles2);
            TS_ASSERT_DELTA(anglehkl, anglemotor, 0.5);
            */

            //====================================================================================  


            #region Refine UB matrix
            // Refine (linear)
            // Set up refine list
            int numhkl = ubcalculator.GetNumberOfReflections();
            bool[] refinelist = new bool[numhkl];
            for (int i = 0; i < refinelist.Length; ++i)
                refinelist[i] = true;

            // Refine
            double errorinput, erroroutput;
            bool isexec = ubcalculator.RefineUBMatrix(refinelist, out errorinput, out erroroutput);
            Console.WriteLine("Error of UB Matrix: {0} ----> {1}.", errorinput, erroroutput);

            #endregion

            #region Refine UB matrix with grid-motor-offset
            // Refine UB matrix on a grid of 2theta, omega, phi and chi, i..e, 
            // considering the uncertainty in motor angles' offset

            double[] bestmotorposoffsets = new double[3];
            double outerror = 0.0;

            bool refinegood = ubcalculator.RefineUBMatrixUnknowMotorAngleOffset(refinelist, delta_2theta, step_2theta, delta_omega, step_omega,
                delta_chi, step_chi, bestmotorposoffsets, out outerror);

            TS_ASSERT(refinegood);
            Console.WriteLine("Best shift: 2theta = {0}, omega = {1}, chi = {2}.", bestmotorposoffsets[0], bestmotorposoffsets[1],
                bestmotorposoffsets[2]);

            #endregion

            return;

            //#region Generate reflection list
            //// Set up
            //double twotheta_min = 15.0;
            //double twotheta_max = 30.0;
            //double omega_min = 6.0;
            //double omega_max = 20.0;
            //double chi_min = -2.0;
            //double chi_max = 2.0;
            //double phi_min = -2.0;
            //double phi_max = 2.0;

            //double[] minmotorpos = new double[4] { twotheta_min, omega_min, phi_min, chi_min };
            //double[] maxmotorpos = new double[4] { twotheta_max, omega_max, phi_max, chi_max };

            //int[] minhklindex = new int[3] { 0, 0, 0 };
            //int[] maxhklindex = new int[3] { 4, 3, 3 };


            //double hshift = 0.5;
            //double kshfit = 0.0;
            //double lshift = 0.0;
            //double[] hklshift = new double[3] { hshift, kshfit, lshift };

            //string condition = "(H+K)%2 == 1 || (H+L)%4 != 0";

            //// Call
            //string macro_commands, processing_results;
            //string suffix_commands = "preset time 60\nscantitle \"test\"\nscanrel k -0.1 0.1 0.01";

            //bool cangencommand = ubcalculator.GenerateSPICEMacroCommands(minhklindex, maxhklindex, minmotorpos, maxmotorpos, hklshift, condition,
            //    suffix_commands, out macro_commands, out processing_results);

            //TS_ASSERT(cangencommand);

            //// Console.WriteLine("Generated command: \n{0}", commands);

            //#endregion

            //return;

            //#region Manipulate with reflections

            //ubcalculator.ExportReflections("C:\\Temp\\exported.txt");

            //int numreflections = ubcalculator.GetNumberOfReflections();
            //bool[] selectflags = new bool[numreflections];
            //for (int i = 0; i < numreflections; ++i)
            //    selectflags[i] = true;

            //selectflags[1] = false;
            //selectflags[3] = false;

            //ubcalculator.CleanReflections(selectflags);
            //numreflections = ubcalculator.GetNumberOfReflections();
            //TS_ASSERT_EQUALS(23, numreflections);

            //ubcalculator.ClearReflections();
            //TS_ASSERT_EQUALS(0, ubcalculator.GetNumberOfReflections());

            //#endregion

            //return;
        }


        #region Testing data generator
        /// <summary>
        /// Create a reflection list file
        /// </summary>
        static private void CreateReflectionListFile(string filename)
        {
            string hklfilebuf = "";
            hklfilebuf += "# hkl angle                                           \n";
            hklfilebuf += "1 0 1     15.01575 7.53825 2.15225 4.6285             \n";
            hklfilebuf += "0 2 0     15.51 7.75325 83.20475 110.348              \n";
            hklfilebuf += "2 2 0     25.85025 12.89225 33.48375 55.311           \n";
            hklfilebuf += "2 -2 0    25.90525 12.89775 -40.24825 46.3595         \n";
            hklfilebuf += "1 4 5     66.89175 33.40775 34.37375 -25.1795         \n";
            hklfilebuf += "-1 -4 5   66.84 33.3265 -22.7015 -50.5575             \n";
            hklfilebuf += "1 -4 5    66.8255 33.40525 -23.92875 -30.3315         \n";
            hklfilebuf += "1 4  -5   66.8415 33.33025 22.88175 129.4435          \n";
            hklfilebuf += "-1 4 5    66.841 33.326 35.722 -47.854                \n";
            hklfilebuf += "-1 -4 -5  66.838 33.43625 -34.24825 154.8215          \n";
            hklfilebuf += "1 -4 -5   66.852 33.37025 -35.68425 132.1465          \n";
            hklfilebuf += "-1  4 -5  66.83975 33.39375 24.10025 149.669          \n";
            hklfilebuf += "0 2 2     26.9045 13.37175 41.408 -35.8375            \n";
            hklfilebuf += "0 -2 2    26.8475 13.35325 -29.34275 -40.9385         \n";
            hklfilebuf += "3 1 1     33.91425 17.02375 12.03225 33.2485          \n";
            hklfilebuf += "-3 1 1    33.878 16.854 18.77575 -110.5615            \n";
            hklfilebuf += "-3 -1 1   33.905 16.856 -7.96725 -108.46              \n";
            hklfilebuf += "-3 -1 -1  33.9345 17.0055 -11.7955 -146.7515          \n";
            hklfilebuf += "3 -1 -1   33.89025 16.89175 -18.57075 69.439          \n";
            hklfilebuf += "3 1 -1    33.9065 16.878 8.20225 71.54                \n";
            hklfilebuf += "-3 1 -1   33.94175 17.022 14.75375 -149.9525          \n";
            hklfilebuf += "3 -1 1    33.87525 17.017 -14.3635 30.0475            \n";
            hklfilebuf += "1 0 5     57.489 28.6855 5.2945 -27.8845              \n";
            hklfilebuf += "-1 0 5    57.49825 28.581 6.53875 -49.291             \n";
            // hklfilebuf += "1 0 -5    57.4495 28.60775 -6.39225 130.709           \n";

            System.IO.File.WriteAllText(filename, hklfilebuf);

            return;
        }

        #endregion

        #region Methods to verify

        static private bool TS_ASSERT(bool realvalue)
        {
            if (!realvalue)
            {
                Console.WriteLine("\n*************************************************************************");
                Console.WriteLine("[Assertion Failure]: Expected value is true; Real value is false.\n");
                Console.WriteLine("*************************************************************************\n");
            }
            return realvalue;
        }

        /// <summary>
        /// Verify whether 2 doubles differs with each other in a certain range
        /// </summary>
        /// <param name="?"></param>
        /// <returns></returns>
        static private bool TS_ASSERT_DELTA(double expectvalue, double realvalue, double delta)
        {
            bool returntrue = false;
            if (Math.Abs(expectvalue - realvalue) <= delta)
                returntrue = true;

            if (!returntrue)
            {
                Console.WriteLine("\n*************************************************************************");
                Console.WriteLine("* Assertion Failure]: Expecteded value is {0} +/- {1};  Real value is {2}. ",
                    expectvalue, delta, realvalue);
                Console.WriteLine("*************************************************************************\n");
            }

            return returntrue;
        }

        static private bool TS_ASSERT_EQUALS(int expectvalue, int realvalue)
        {
            bool returntrue = (expectvalue == realvalue);
            if (!returntrue)
            {
                Console.WriteLine("\n*************************************************************************");
                Console.WriteLine("* [Assertion Failure]: Expected value is {0}; Real value is {1}",
                    expectvalue, realvalue);
                Console.WriteLine("*************************************************************************\n");
            }

            return returntrue;
        }

        #endregion

    }
}
