using System;
using System.Collections.Generic;
using System.Text;
using DAS.UBMatrix;

namespace Testers
{
    /// <summary>
    /// Main application for testing without LabVIEW.  I
    /// It calls the class's methods defined for LabVIEW 
    /// </summary>
    public class StandAloneApp
    {
        static public void Main_StandAlone(string[] args) // 
        {
            Console.WriteLine("********************************************************");
            Console.WriteLine("* UB Matrix Standalone Calculator");
            Console.WriteLine("********************************************************\n");  

            // Justify input arguments
            if (args.Length < 3)
            {
                Console.WriteLine("Format: Exec-name [Input Unit Cell File] [Input Reflection File] [Wavelength]");
                Console.WriteLine("\tUnit Cell File: a  b  c  alpha beta gamma (free format)");
                Console.WriteLine("\tInput Reflection File: \n");
                Console.WriteLine("\t\t# HKL Motor-Positoin / Motor-Position HKL /Motor-Position  (motor position: 2theta omega chi phi)");
                System.Environment.Exit(0);
            }

            // Initialize UB matrix object
            const int maxreflection = 50;
            double wavelength = Convert.ToDouble(args[2]);
            if (wavelength <= 0)
            {
                Console.WriteLine("Input wavelength is unphysical.  It is equal to {0}. ", wavelength);
                System.Environment.Exit(0);
            }

            LabVIEWAdaptor lvadapt = new LabVIEWAdaptor(wavelength, maxreflection);

            // Process input arguments
            string unitcellfilename = args[0];
            string reflectfilename = args[1];

            // Parse unit cell file
            double a, b, c, alpha, beta, gamma;
            UBMatrixUserInterface.ReadCellFile(unitcellfilename, out a, out b, out c, out alpha, out beta, out gamma);
            lvadapt.SetLatticeParameters(a, b, c, alpha, beta, gamma);

            // Import reflections: output -> m_inputHKLArrays/m_inputMotorAnglesArray
            int numreflection = lvadapt.ImportReflectionsFromFile(reflectfilename);
            if (numreflection <= 0)
            {
                Console.WriteLine("Failed to import reflection file or there is no reflections in reflection file. ");
                System.Environment.Exit(-1);
            }

            // Convert input motor angles (Brian omega) to m_observedMotorAnglesArray (Busing omega)
            int[] temphkl = new int[3];
            double[] tempmotorangles = new double[4];
            for (int i = 0; i < lvadapt.NumberOfImportedReflections; i++)
            {
                lvadapt.GetImportedReflection(i, tempmotorangles, temphkl);
                // lvadapt.SetSubObservedMotorAngles(i, tempmotorangles);
            }
            bool consolecontinue = true;

            while (consolecontinue)
            {
                Console.WriteLine("Input Options: \n");
                Console.WriteLine("\t(1/c) Calculate UB Matrix              ");
                Console.WriteLine("\t(2/r) Refine UB Matrix                 ");
                Console.WriteLine("\t(3/m) Calcualte Motor Position         ");
                Console.WriteLine("\t(4/a) Print All Reflections            ");
                Console.WriteLine("\t(5/u) Print UB Matrix                  ");
                Console.WriteLine("\t(6/l) Print Lattice Parameters         ");
                Console.WriteLine("\t(7/h) Calcualte (HKL) from motor angle ");
                Console.WriteLine("\t(9/g) Generate a list of reflections from UB Matrix ");
                Console.WriteLine("\t(e) Exit");
                Console.Write("> ");
                string inputchoise = Console.ReadLine();
                switch (inputchoise[0])
                {
                    case '1':
                    case 'c':
                        CalculateUBMatrix(lvadapt);
                        break;
                    case '2':
                    case 'r':
                        RefineUBMatrix(lvadapt);
                        break;
                    case '3':
                    case 'm':
                        CalculateMotorPosition(lvadapt);
                        break;
                    case '4':
                    case 'a':
                        PrintInputReflections(lvadapt);
                        break;
                    case '5':
                    case 'u':
                        PrintUBMatrix(lvadapt);
                        break;
                    case 'e':
                        consolecontinue = false;
                        break;
                    case '6':
                    case 'l':
                        PrintLatticeParameters(lvadapt);
                        break;
                    case '7':
                    case 'h':
                        CalcualteHKLFromMotorAngle(lvadapt);
                        break;
                    case '8':
                        // Hidden option
                        CheckAllInputIndicents(lvadapt);
                        break;
                    case '9':
                        // (New 2013.12) 
                        GenerateReflections(lvadapt);
                        break;
                    default:
                        Console.WriteLine("Option {0} Is Not Accepted", inputchoise[0]);
                        break;
                }

            }

            Console.WriteLine("\nQuit UB Matrix Calculator");

            return;
        }

        /// <summary>
        /// Calculate UB Matrix from input reflections
        /// </summary>
        /// <param name="lvadp"></param>
        static private void CalculateUBMatrix(LabVIEWAdaptor lvadp)
        {
            // Prompt for inputs
            Console.WriteLine("");
            Console.WriteLine("# Calculate UB Matrix From 2 Reflections.");
            Console.Write("> Input 2 Reflection's Number: ");

            string userinput = Console.ReadLine();
            Console.WriteLine("$ User input : {0}", userinput);
            int[] reflections = new int[2];
            StringUtility.ParseInteger(userinput, 2, reflections);

            if (reflections[0] == reflections[1])
            {
                Console.WriteLine("Invalid Error!  Must be two different reflections");
                return;
            }
            else
            {
                Console.WriteLine("UB matrix is to be built on reflection #{0} and #{1}. ", reflections[0],
                    reflections[1]);
            }

            double[] ma1 = new double[4];
            int[] ml1 = new int[3];
            double[] ma2 = new double[4];
            int[] ml2 = new int[3];

            lvadp.GetImportedReflection(reflections[0], ma1, ml1);
            lvadp.GetImportedReflection(reflections[1], ma2, ml2);

            lvadp.GenerateUBMatrix2Reflections(ml1, ma1, ml2, ma2);

            double millerangle = lvadp.CalculateAngleBW2MillerIndices(ml1, ml2);
            double motorangle = lvadp.CalculateAngleBW2MotorAngles(ma1, ma2);

            Console.WriteLine("2 Reflections Angle:\n\tAngle between Miller Indices = {0}\n\tAngle between Motor Angles = {1}",
                millerangle, motorangle);

            double[] edges = new double[3];
            double[] angles = new double[3];
            lvadp.GetLatticeParameters(edges, angles);
            Console.WriteLine("Lattice parameters are ({0}, {1}, {2}) and ({3}, {4}, {5})",
                edges[0], edges[1], edges[2], angles[0], angles[1], angles[2]);

            return;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="lvadp"></param>
        static private void RefineUBMatrix(LabVIEWAdaptor lvadp)
        {
            // Prompt
            Console.WriteLine("\n> Refine UB Matrix\nInput Reflections's Number Tended To Be Excluded: ");
            string userinput = Console.ReadLine();

            int[] exclusionlist = new int[50];
            int numexclused = StringUtility.ParseIndefinateInteger(userinput, exclusionlist);
            for (int i = 0; i < numexclused; i++)
            {
                Console.WriteLine("{0}:  {1}", i, exclusionlist[i]);
            }

            // 2. Form refine including flag
            bool[] refineflags = new bool[lvadp.NumberOfImportedReflections];
            for (int i = 0; i < refineflags.Length; i++)
                refineflags[i] = true;
            for (int i = 0; i < numexclused; i++)
                refineflags[exclusionlist[i]] = false;

            // 3. Refine
            double inerror, outerror;
            // TODO Refine UB matrix in a loop with different offset to motor angles.  Find best one (outerror)
            lvadp.RefineUBMatrix(refineflags, out inerror, out outerror);

            return;
        }

        /// <summary>
        /// Calculate motor position of a reflection input from keyboard
        /// </summary>
        /// <param name="lvadp"></param>
        static private void CalculateMotorPosition(LabVIEWAdaptor lvadp)
        {
            Console.Write("\n> Calculate Motor Position From UB Matrix.  Enter (H, K, L) Below\n> ");
            string userinput = Console.ReadLine();

            int[] hkl = new int[3];
            StringUtility.ParseInteger(userinput, 3, hkl);

            double[] mangles = new double[4];
            lvadp.CalculateMotorAnglesFromHKL(hkl, mangles);

            Console.WriteLine("HKL = ({0}, {1}, {2}) ==> Motor Angles (2theta, omega, chi, phi) = {3}, {4}, {5}, {6}",
                hkl[0], hkl[1], hkl[2], mangles[0], mangles[1], mangles[2], mangles[3]);

            return;         
        }

        static private void PrintInputReflections(LabVIEWAdaptor lvadp)
        {
            Console.WriteLine("\n> User Input Reflections:  Total Number = {0}", lvadp.NumberOfImportedReflections);
            int[] millerindex = new int[3];
            double[] motorangles = new double[4];
            string buf;
            for (int i = 0; i < lvadp.NumberOfImportedReflections; i++)
            {
                lvadp.GetImportedReflection(i, motorangles, millerindex);
                buf = String.Format("[{0}]:\t", i);
                for (int j = 0; j < 4; j++)
                    buf += String.Format("{0:000.000}  ", motorangles[j]);
                buf += "\t\t";
                for (int j = 0; j < 3; j++)
                    buf += String.Format("{0} ", millerindex[j]);
                Console.WriteLine(buf);
            }
        }

        static private void PrintUBMatrix(LabVIEWAdaptor lvadp)
        {
            Console.WriteLine("\n> UB Matrix:");
            string buf = "";
            double[] ubrow = new double[3];
            for (int i = 0; i < 3; i++)
            {
                lvadp.GetUBMatrix(ubrow, i);
                for (int j = 0; j < 3; j++)
                {
                    buf += String.Format("{0:000.0000  }", ubrow[j]);
                }
                buf += "\n";
            }
            Console.WriteLine(buf);
        }

        static private void PrintLatticeParameters(LabVIEWAdaptor lvadp)
        {
            Console.WriteLine("\n> Lattice Parameters:");
            double[] edges = new double[3];
            double[] angles = new double[3];
            lvadp.GetLatticeParameters(edges, angles);
            for (int i = 0; i < 3; i++)
            {
                Console.WriteLine("{0:0.000}\t{1}", angles[i], edges[i]);
            }
            Console.WriteLine();
        }

        /// <summary>
        /// Calcualte (H, K, L) from input motor angle
        /// </summary>
        /// <param name="lvadapt"></param>
        static public void CalcualteHKLFromMotorAngle(LabVIEWAdaptor lvadapt)
        {
            Console.WriteLine("\n> Input Motor Angle (Option 1: 2theta omega phi chi; Option 2: Input Reflection Index");
            string input = Console.ReadLine();
            string[] terms = input.Split(new char[] { ' ', '\t', ',' });
            double[] importedmotorangles = new double[4];
            double[] originalhkl = new double[3];
            int[] hkl = new int[3];
            if (terms.Length == 1)
            {
                int index = Convert.ToInt32(terms[0]);
                lvadapt.GetImportedReflection(index, importedmotorangles, hkl);
                lvadapt.CalculateMillerIndexFromMotorAngle(hkl, originalhkl, importedmotorangles);
            }
            else
            {
                throw new Exception("To Be Implemented Soon For Input Motor Angles");
            }

            Console.WriteLine("Calcualte H, K, L: ({0}, {1}, {2})", hkl[0], hkl[1], hkl[2]);

            return;
        }

        /// <summary>
        /// Generate a list of relfections and their motor positions
        /// </summary>
        /// <param name="lvadapt"></param>
        static public void GenerateReflections(LabVIEWAdaptor lvadapt)
        {
            // Get input from user
            Console.Write("> Limit on (HKL: H_min, H_max, K_min, K_max, L_min, L_max): ");
            string hkllimitstr = Console.ReadLine();
            Console.Write("> Motor angles limits (2thta_min, 2theta_max, omega_min, omega_max, phi_min, phi_max, chi_min, chi_max): ");
            string motorposlimitstr = Console.ReadLine();
            Console.Write("> Conditions: ");
            string conditionstr = Console.ReadLine();

            Console.WriteLine("\n**********\nUser input HKL limit:            {0}\n           Motor position limit: {1}\n           Condition:            {2}\n**********", 
                hkllimitstr, motorposlimitstr, conditionstr);


            // Parse HKL range          
            int[] vecHKLLimits = new int[6];
            StringUtility.ParseInteger(hkllimitstr, 6, vecHKLLimits);
#if false
            MillerIndices minindex = new MillerIndices();
            MillerIndices maxindex = new MillerIndices();
            minindex.H = vecHKLLimits[0];
            maxindex.H = vecHKLLimits[1];
            minindex.K = vecHKLLimits[2];
            maxindex.K = vecHKLLimits[3];
            minindex.L = vecHKLLimits[4];
            maxindex.L = vecHKLLimits[5];
#else
            int[] minindex = new int[3];
            int[] maxindex = new int[3];
            minindex[0] = vecHKLLimits[0];
            maxindex[0] = vecHKLLimits[1];
            minindex[1] = vecHKLLimits[2];
            maxindex[1] = vecHKLLimits[3];
            minindex[2] = vecHKLLimits[4];
            maxindex[2] = vecHKLLimits[5];
#endif

            // Parse motor limit
            double[] vecMotorPosLimits = new double[8];
            StringUtility.ParseDoubles(motorposlimitstr, 8, vecMotorPosLimits);

#if false
            MotorIncidentAngles minmotorpos = new MotorIncidentAngles();
            MotorIncidentAngles maxmotorpos = new MotorIncidentAngles();
            minmotorpos.twotheta = vecMotorPosLimits[0];
            maxmotorpos.twotheta = vecMotorPosLimits[1];
            minmotorpos.omega = vecMotorPosLimits[2];
            maxmotorpos.omega = vecMotorPosLimits[3];
            minmotorpos.phi = vecMotorPosLimits[4];
            maxmotorpos.phi = vecMotorPosLimits[5];
            minmotorpos.chi = vecMotorPosLimits[6];
            maxmotorpos.chi = vecMotorPosLimits[7];
#else
            double[] minmotorpos = new double[4];
            double[] maxmotorpos = new double[4];
            minmotorpos[0] = vecMotorPosLimits[0];
            maxmotorpos[0] = vecMotorPosLimits[1];
            minmotorpos[1] = vecMotorPosLimits[2];
            maxmotorpos[1] = vecMotorPosLimits[3];
            minmotorpos[2] = vecMotorPosLimits[4];
            maxmotorpos[2] = vecMotorPosLimits[5];
            minmotorpos[3] = vecMotorPosLimits[6];
            maxmotorpos[3] = vecMotorPosLimits[7];
#endif

            // Parse condition
            if (conditionstr.Length > 0)
            {
                string error_message;
                bool conditionvalid = lvadapt.CheckConditionValidity(conditionstr, out error_message);
                if (!conditionvalid)
                {
                    Console.WriteLine("Condition {0} is not valid.  Log: {1}", conditionstr, error_message);
                    return;
                }
            }

            // Call LabVIEWAdaptor's 
            string info, finallist;
            double[] hklshift = new double[3] { 0, 0, 0 };
            lvadapt.GenerateReflectionList(minindex, maxindex, minmotorpos, maxmotorpos, hklshift, conditionstr, out finallist, out info);
            Console.WriteLine("Generated HKL/Motor Positions: \n{0}", finallist);
            Console.WriteLine("Information FYI: \n{0}", info);

            return;
        }

        /// <summary>
        /// Check against all input incident angles with calculated UB matrix and report it
        /// </summary>
        /// <param name="lvadapt"></param>
        /// <returns></returns>
        static private void CheckAllInputIndicents(LabVIEWAdaptor lvadapt)
        {
            double[] calculatedmotorangles = new double[4];
            double[] importedmotorangles = new double[4];
            double[] originalhkl = new double[3];
            int[] hkl = new int[3];
            int[] hkl2 = new int[3];
            string outputbuf = "";
            Console.WriteLine("Motor Angles (2theta, omega, chi, phi):");
            for (int i = 0; i < lvadapt.NumberOfImportedReflections; i++)
            {
                lvadapt.GetImportedReflection(i, importedmotorangles, hkl);
                lvadapt.CalculateMillerIndexFromMotorAngle(hkl2, originalhkl, importedmotorangles);
                lvadapt.CalculateMotorAnglesFromHKL(hkl, calculatedmotorangles);

                bool samehkl = true;
                for (int j = 0; j < 3; j ++)
                    if (hkl[j] != hkl2[j])
                    {
                        samehkl = false;
                        break;
                    }

                if (samehkl)
                {
                    outputbuf = String.Format("HKL = ({0,-3}, {1,-3}, {2,-3}) [Original = ({3,-4:0.00}, {4,-4:0.00}, {5,-4:0.00})]\n",
                        hkl[0], hkl[1], hkl[2], originalhkl[0], originalhkl[1], originalhkl[2]);
                }
                else
                {
                    outputbuf = String.Format("Wrong HKL Calculated!  User HKL = ({0}, {1}, {2})  vs. ({3}, {4}, {5}) Calcualted From Motor Angles \n",
                        hkl[0], hkl[1], hkl[2], hkl2[0], hkl2[1], hkl2[2]);
                }
                outputbuf += String.Format("\tObserved:  \t{0,-10:0.000} {1,-10:0.000} {2,-10:0.000} {3,-10:0.000}\t{4,-10:0.000} {5,-10:0.000} {6,-10:0.000} {7,-10:0.000}\n",
                    importedmotorangles[0], importedmotorangles[1],
                    importedmotorangles[2], importedmotorangles[3],
                    (importedmotorangles[0] - calculatedmotorangles[0]),
                    (importedmotorangles[1] - calculatedmotorangles[1]),
                    (importedmotorangles[2] - calculatedmotorangles[2]),
                    (importedmotorangles[3] - calculatedmotorangles[3]));
                outputbuf += String.Format("\tCalculated:\t{0,-10:0.000} {1,-10:0.000} {2,-10:0.000} {3,-10:0.000}",
                    calculatedmotorangles[0], calculatedmotorangles[1],
                    calculatedmotorangles[2], calculatedmotorangles[3]);
                Console.WriteLine(outputbuf);
            }

        }

    }
}
