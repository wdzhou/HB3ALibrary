using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using DAS.UBMatrix;
using System.Threading;

namespace Testers
{
    class WorkFlowInteractive
    {
        /// <summary>
        /// Test Plan:
        /// 1. Vector
        /// 2. Matrix
        /// 3. UBMatrix
        /// </summary>
        /// <param name="args"></param>
        static void MainWorkflow(string[] args)
        {
            if (args.Length < 4)
            {
                Console.WriteLine("Format:  [Work Flow Choice 1/2/3]  [Input Unit Cell File]  [Input Reflection File]  [Wavelength]");
                Console.WriteLine("Workflow 1:  Regular Interactive Test (Least Square Fitting May Not Work)");
                Console.WriteLine("         2:  Tester on Least Square Fitting (Use 2 reflections to build UB; Refine one by one for the rest");
                Console.WriteLine("         3:  Tester of Auto Indexing");
                Console.WriteLine("         4:  Test LabVIEW Adaptor");
                Console.WriteLine("         5:  Test LabVIEW Adaptor 2010.01.11  [reflection file name] [random]  [random]");
                Environment.Exit(0);
            }

            int workflowoption = Convert.ToInt16(args[0]);
            string cellfilename = args[1];
            string reflectionfilename = args[2];
            double wavelength = Convert.ToDouble(args[3]);

            if (workflowoption == 1)
                //InteractiveDetermineUB();
                WorkFlow1(wavelength, cellfilename, reflectionfilename);
            else if (workflowoption == 2)
                WorkFlow2(wavelength, cellfilename, reflectionfilename);
            else if (workflowoption == 3)
                AutoIndexing1(wavelength, cellfilename, reflectionfilename);
            else if (workflowoption == 4)
                TestLabViewAdaptor(reflectionfilename);
            else if (workflowoption == 5)
            {
                reflectionfilename = args[1];
                TestLabVIEWAdaptorVer20(reflectionfilename);
            }
            else
                Console.WriteLine("Option {0} Unsupported!", workflowoption);

            // -1: waiting for user finish
            Console.WriteLine("Press Any Key To Finish");
            Console.ReadLine();

            return;
        }

        /// <summary>
        /// Interactive workflow
        /// </summary>
        /// <param name="wavelength"></param>
        /// <param name="unitcellfilename"></param>
        /// <param name="reflectionfilename"></param>
        static public void WorkFlow1(double wavelength, string unitcellfilename, string reflectionfilename)
        {
            // 1. Get Input and Generate UB matrix
            double lambda;
            int numinputreflections;
            UnitCell myunitcell = new UnitCell();
            MotorIncidentAngles[] mincidentangles = new MotorIncidentAngles[100];
            for (int i = 0; i < 100; i++)
                mincidentangles[i] = new MotorIncidentAngles();

            ReadUnitCellFile(unitcellfilename, myunitcell);
            if (wavelength <= 0.0)
                lambda = 1.56;
            else
                lambda = wavelength;
            ReadReflectionMotorAngles(reflectionfilename, mincidentangles, out numinputreflections);

            UBMatrix myubmatrix = new UBMatrix(lambda);
            myubmatrix.SetDirectUnitCell(myunitcell);
            myubmatrix.StoreUnitCell();

            // 2. Print out input information
            Console.WriteLine("Input Unit Cell: {0}", myunitcell);
            Console.WriteLine("Lambda = {0}    Default = 1.56", lambda);
            for (int i = 0; i < numinputreflections; i++)
                Console.WriteLine("{1}:   {0}", mincidentangles[i], i);

            // 2.5 List all pairs of angles
            /*
            int numpairs = (numinputreflections - 1) * numinputreflections / 2;
            double[] ang12list = new double[numpairs];
            int[] pair1list = new int[numpairs];
            int[] pair2list = new int[numpairs];
            myubmatrix.ListPairAngles(mincidentangles, numinputreflections, ang12list, pair1list, pair2list);
            
            int numshowpair = 120;
            for (int i = numpairs - 1; i >= 0; i--)
            {
                Console.WriteLine("u_phi angle b/w {0} and {1}\t= {2:000.0000}",
                    pair1list[i], pair2list[i], ang12list[i]);
                if (numpairs - i == numshowpair)
                    break;
            }
            */

            bool continueloop = true;
            int[] refindex = new int[2];
            MillerIndices[] refhkl = new MillerIndices[2];
            for (int i = 0; i < 2; i++)
            {
                refindex[i] = -1;
                refhkl[i] = new MillerIndices();
            }


            while (continueloop)
            {

                Console.WriteLine("Type Command Type: (1) Angle b/w 2 reflections  (2) Calculate UB Matrix");
                Console.WriteLine("                   (3) Exit                     (4) Refine (may not work)");
                string stropt = Console.ReadLine().Split()[0];
                int option = Convert.ToInt16(stropt);

                switch (option)
                {
                    case 1:
                        // Calculate angles b/w 2 input
                        Console.WriteLine("Input Reflection Number");
                        string infostring = Console.ReadLine();
                        string[] terms = infostring.Split();
                        int r1 = Convert.ToInt16(terms[0]);
                        int r2 = Convert.ToInt16(terms[1]);

                        Vector uH1 = new Vector(3);
                        Vector uH2 = new Vector(3);
                        MotorIncidentAngles.UnitVectorFromAngle(mincidentangles[r1], uH1);
                        MotorIncidentAngles.UnitVectorFromAngle(mincidentangles[r2], uH2);
                        double dotp = Vector.DotProduct(uH1, uH2);
                        double ang12 = DasMath.AcosD(dotp);
                        Console.WriteLine("Angle = {0:0.000}  For Reflection {1} and {2}",
                            ang12, r1, r2);
                        break;

                    case 2:
                        // Calculate UB matrix
                        // 3. User input HKL for generating UB matrix
                        myubmatrix.RestoreUnitCell();
                        Console.WriteLine("Unit Cell Restored");

                        Console.WriteLine("Input (H, K, L) for 2 selected reflections\nFormat: Index of Motor Angle, H, K, L");

                        int count = 0;
                        while (count < 2)
                        {
                            bool valid = true;
                            try
                            {
                                string hklline = Console.ReadLine();
                                string[] term1s = hklline.Split();
                                refindex[count] = Convert.ToInt16(term1s[0]);
                                refhkl[count] = new MillerIndices();
                                for (int i = 0; i < 3; i++)
                                    refhkl[count].V[i] = Convert.ToInt16(term1s[i + 1]);
                            }
                            catch (Exception e)
                            {
                                Console.WriteLine("Encountered exception %s.", e.Message);
                                valid = false;
                            }

                            if (valid)
                                count++;
                        }

                        for (int i = 0; i < 2; i++)
                            Console.WriteLine("Reflection {0} Picked:  {1}  {2}",
                                i, mincidentangles[refindex[i]], refhkl[i]);


                        // 4. Construct UB matrix
                        int ierr;
                        myubmatrix.MakeUB(refhkl[0], refhkl[1], mincidentangles[refindex[0]], mincidentangles[refindex[1]], out ierr);
                        myubmatrix.Print();

                        // 4.1 
                        Console.WriteLine("*****************   Test 954    ***************");
                        CheckAllReflectionByUBMatrix(myubmatrix, mincidentangles, numinputreflections);
                        Console.WriteLine("************************************************");
                        Console.WriteLine("Press Any Key To Continue!");
                        Console.ReadLine();

                        break;
                    case 3:
                        // Exit
                        continueloop = false;
                        break;

                    case 4:
                        // Refine

                        // 5. Continue to refine
                        // 5.1 Construct a reflection list
                        Reflection[] mreflections = new Reflection[100];
                        int numreflections = 2;
                        for (int i = 0; i < 2; i++)
                        {
                            mreflections[i] = new Reflection();
                            mreflections[i].MotorAngles = mincidentangles[refindex[i]];
                            mreflections[i].MillerIndex = refhkl[i];
                        }

                        for (int i = 0; i < numinputreflections; i++)
                        {
                            if (i != refindex[0] && i != refindex[1])
                            {
                                // not two base reflections

                                // 1. Test reflection is good for UB matrix to calculate HKL
                                MillerIndices mi = new MillerIndices();
                                Vector uphi = new Vector(3);
                                double diff2tehta, diffmangle;
                                myubmatrix.CalculateHKLFromIncidentAngle(mincidentangles[i], mi, uphi);
                                bool hklvalid = myubmatrix.CheckError(mincidentangles[i], mi, 1.0, 1.0, out diff2tehta, out diffmangle);
                                if (hklvalid)
                                {
                                    mreflections[numreflections] = new Reflection();
                                    mreflections[numreflections].MotorAngles = mincidentangles[i];
                                    mreflections[numreflections].MillerIndex = mi;
                                    numreflections++;

                                    // b. Refine
                                    bool canfit = myubmatrix.RefineUBLinearLSF(mreflections, numreflections);
                                    if (canfit)
                                    {
                                        Console.WriteLine("*****************   Test 1217 Test Case Only!    ***************");
                                        CheckAllReflectionByUBMatrix(myubmatrix, mincidentangles, numinputreflections);
                                        Console.WriteLine("************************************************");
                                    }
                                    else
                                        numreflections--;

                                }
                                else
                                {
                                    Console.WriteLine("Bad Reflections!");
                                }



                                Console.WriteLine("Press Any Key For Next Reflection");
                                Console.ReadLine();
                            }
                        }
                        break;

                    default:
                        Console.WriteLine("Option {0} Not Supported", option);
                        break;
                }
            }

            return;
        }

        /// <summary>
        /// Refinement workflow
        /// </summary>
        /// <param name="wavelength"></param>
        /// <param name="unitcellfilename"></param>
        /// <param name="reflectionfilename"></param>
        static public void WorkFlow2(double wavelength, string unitcellfilename, string reflectionfilename)
        {
            // 1. Get Input and Generate UB matrix
            double lambda;
            int numinputreflections;
            UnitCell myunitcell = new UnitCell();
            MotorIncidentAngles[] mincidentangles = new MotorIncidentAngles[100];
            for (int i = 0; i < 100; i++)
                mincidentangles[i] = new MotorIncidentAngles();
          
            ReadUnitCellFile(unitcellfilename, myunitcell);
            if (wavelength <= 0.0)
                lambda = 1.56;
            else
                lambda = wavelength;
            ReadReflectionMotorAngles(reflectionfilename, mincidentangles, out numinputreflections);

            UBMatrix myubmatrix = new UBMatrix(lambda);
            myubmatrix.SetDirectUnitCell(myunitcell);
            myubmatrix.StoreUnitCell();

            // 2. Print out input information
            Console.WriteLine("Input Unit Cell: {0}", myunitcell);
            Console.WriteLine("Lambda = {0}    Default = 1.56", lambda);
            for (int i = 0; i < numinputreflections; i++)
                Console.WriteLine("{1}:   {0}", mincidentangles[i], i);

            // 2.5 List all pairs of angles
            int numpairs = (numinputreflections-1)*numinputreflections/2;
            double[] ang12list = new double[numpairs];
            int[] pair1list = new int[numpairs];
            int[] pair2list = new int[numpairs];
            myubmatrix.ListPairAngles(mincidentangles, numinputreflections, ang12list, pair1list, pair2list);

            int numshowpair = 80;
            for (int i = numpairs - 1; i >= 0; i--)
            {
                Console.WriteLine("u_phi angle b/w ({0}, {1})  \t = {2:000.000}",
                    pair1list[i], pair2list[i], ang12list[i]);
                if (numpairs - i == numshowpair)
                    break;
            }

            // 3. User input HKL for generating UB matrix
            Console.WriteLine("Input (H, K, L) for 2 selected reflections\nFormat: Index of Motor Angle, H, K, L");
            int[] refindex = new int[2];
            MillerIndices[] refhkl = new MillerIndices[2];

            int count = 0;
            while (count < 2)
            {
                bool valid = true;
                try
                {
                    string hklline = Console.ReadLine();
                    string[] terms = hklline.Split();
                    refindex[count] = Convert.ToInt16(terms[0]);
                    refhkl[count] = new MillerIndices();
                    for (int i = 0; i < 3; i++)
                        refhkl[count].V[i] = Convert.ToInt16(terms[i + 1]);
                }
                catch (Exception e)
                {
                    Console.WriteLine("Encounted exception %s. ", e.Message);
                    valid = false;
                }

                if (valid)
                    count++;
            }

            for (int i = 0; i < 2; i++)
                Console.WriteLine("Reflection {0} Picked:  {1}  {2}",
                    i, mincidentangles[refindex[i]], refhkl[i]);
            
            // 4. Construct UB matrix
            int ierr;
            myubmatrix.MakeUB(refhkl[0], refhkl[1], mincidentangles[0], mincidentangles[1], out ierr);
            myubmatrix.Print();

            // 5. Continue to refine
            // 5.1 Construct a reflection list
            Reflection[] mreflections = new Reflection[100];
            int numreflections = 2;
            for (int i = 0; i < 2; i ++){
                mreflections[i] = new Reflection();
                mreflections[i].MotorAngles = mincidentangles[refindex[i]];
                mreflections[i].MillerIndex = refhkl[i];
            }

            for (int i = 0; i < numinputreflections; i++)
            {
                if (i != refindex[0] && i != refindex[1])
                {
                    // not two base reflections

                    // 1. Test reflection is good for UB matrix to calculate HKL
                    MillerIndices mi = new MillerIndices();
                    Vector uphi = new Vector(3);
                    double diff2tehta, diffmangle;
                    myubmatrix.CalculateHKLFromIncidentAngle(mincidentangles[i], mi, uphi);
                    bool hklvalid = myubmatrix.CheckError(mincidentangles[i], mi, 1.0, 1.0, out diff2tehta, out diffmangle);
                    if (hklvalid)
                    {
                        mreflections[numreflections] = new Reflection();
                        mreflections[numreflections].MotorAngles = mincidentangles[i];
                        mreflections[numreflections].MillerIndex = mi;
                        numreflections++;

                        // b. Refine
                        bool canfit = myubmatrix.RefineUBLinearLSF(mreflections, numreflections);
                        if (canfit)
                        {
                            Console.WriteLine("*****************   Test 1217 Test Case Only!    ***************");
                            CheckAllReflectionByUBMatrix(myubmatrix, mincidentangles, numinputreflections);
                            Console.WriteLine("************************************************");
                        }
                        else
                        {
                            // UB1138:  Testing Whether The Difference To Put More Reflections In Making Refinement Better
                            numreflections += 0;
                            // numreflections--;
                        }
                    }
                    else
                    {
                        Console.WriteLine("Bad Reflections!");
                    }

                    

                    Console.WriteLine("Press Any Key For Next Reflection");
                    Console.ReadLine();
                }
            }

            return;
        }

        /// <summary>
        /// Test Auto Indexing 
        /// </summary>
        /// <param name="wavelength"></param>
        /// <param name="unitcellfilename"></param>
        /// <param name="reflectionfilename"></param>
        static public void AutoIndexing1(double wavelength, string unitcellfilename, string reflectionfilename)
        {
            // 1. Get Input and Generate UB matrix
            double lambda;
            int numinputreflections;
            UnitCell myunitcell = new UnitCell();
            MotorIncidentAngles[] mincidentangles = new MotorIncidentAngles[100];
            for (int i = 0; i < 100; i++)
                mincidentangles[i] = new MotorIncidentAngles();

            ReadUnitCellFile(unitcellfilename, myunitcell);
            if (wavelength <= 0.0)
                lambda = 1.56;
            else
                lambda = wavelength;
            ReadReflectionMotorAngles(reflectionfilename, mincidentangles, out numinputreflections);

            UBMatrix myubmatrix = new UBMatrix(lambda);
            myubmatrix.SetDirectUnitCell(myunitcell);
            myubmatrix.StoreUnitCell();

            // 2. Print out input information
            Console.WriteLine("Input Unit Cell: {0}", myunitcell);
            Console.WriteLine("Lambda = {0}    Default = 1.56", lambda);
            for (int i = 0; i < numinputreflections; i++)
                Console.WriteLine("{1}:   {0}", mincidentangles[i], i);

            myubmatrix.GenerateUBFromUnIndexedReflectiuons(mincidentangles, numinputreflections);

        }


        static public void TestLabViewAdaptor(string filename)
        {
            Console.WriteLine("1...");
            LabVIEWAdaptor lva = new LabVIEWAdaptor(1.56, 50);
            Console.WriteLine("2...");
            lva.SetLatticeParameters(1, 1, 1, 90, 90, 90);
            Console.WriteLine("3...");
            int numreflection = lva.ImportReflectionsFromFile(filename);
            Console.WriteLine("Number of Reflection: {0}", numreflection);
            double[] motorangles = new double[4];
            int[] millerindices = new int[3];
            for (int i = 0; i < numreflection; i++)
                lva.GetImportedReflection(i, motorangles, millerindices);
            return;
        }

        static public void TestLabVIEWAdaptorVer20(string filename)
        {
            // Init, Wavelength, Parameter
            LabVIEWAdaptor lva = new LabVIEWAdaptor(50);


            // Set Value 
#if false
            lva.SetWaveLength(1.57);
            lva.SetLatticeParameters(4.753, 5.72, 4.968, 90, 90.08, 90);
            int[] hkl1 = new int[] { 0, 1, 0 };
            double[] ma1 = new double[] { 15.71, 7.73, -7.07, -50.71 };
            int[] hkl2 = new int[] { 0, 1, 1 };
            double[] ma2 = new double[] { 23.78, 11.9, -36.12, -93.97 };
#endif

            lva.SetWaveLength(1.01);
            lva.SetLatticeParameters(6.1007, 3.5235, 11.898, 90, 106.056, 90);
            int[] hkl1 = new int[]{1, 0, 7};
            double[] ma1 = new double[] { 39.55, 19.3497, 75.03, 176.40 };
            int[] hkl2 = new int[] { 4, 0, 0 };
            double[] ma2 = new double[] { 40.2, 21.3278, 19.239, 133.34 }; 

            // Calculate UB matrix from 2 reflections
            double a2ml = lva.CalculateAngleBW2MillerIndices(hkl1, hkl2);
            double a2ma = lva.CalculateAngleBW2MotorAngles(ma1, ma2);
            Console.WriteLine("Angles Between (HKL) = {0}    Between Motor Angles = {1}", a2ml, a2ma);
            lva.GenerateUBMatrix2Reflections(hkl1, ma1, hkl2, ma2);

            Console.WriteLine("\nUB Matrix:");
            double[] row = new double[3];
            for (int i = 0; i < 3; i++)
            {
                lva.GetUBMatrix(row, i);
                for (int j = 0; j < 3; j++)
                    Console.Write("{0}\t", row[j]);
                Console.WriteLine();
            }

            // Loop on calcualte motor angles
            int[] inputhkl = new int[3];
            double[] mangles = new double[4];
            while (true)
            {
                Console.Write("Input (H, K, L): ");
                string hklstring = Console.ReadLine();

                // 1. Quit loop or not? 
                if (hklstring.Contains("quit") || hklstring.Contains("exit"))
                {
                    break;
                }

                string[] terms = hklstring.Split(); 
                for (int i = 0; i < 3; i++)
                    inputhkl[i] = Convert.ToInt32(terms[i]);
                lva.CalculateMotorAnglesFromHKL(inputhkl, mangles);
                Console.Write("Result = ");
                for (int i = 0; i < 4; i++)
                    Console.Write("{0}\t", mangles[i]);
                Console.WriteLine();
            }

            // Import
            int numref = lva.ImportReflectionsFromFile(filename);
            double[][] obsmotorangles = new double[numref][];
            for (int i = 0; i < numref; i ++){
                obsmotorangles[i] = new double[4];
            }

            // Refine
            bool[] refflags = new bool[50];
            double[][] covarmatrix = new double[3][];
            for (int i = 0; i < 3; i++)
                covarmatrix[i] = new double[3];
            int[] millerindex = new int[3];
            for (int i = 0; i < numref; i++)
            {
                refflags[i] = true;
                lva.GetImportedReflection(i, obsmotorangles[i], millerindex);
                // lva.SetSubObservedMotorAngles(i, obsmotorangles[i]);
            }
            double inerror, outerror;
            lva.RefineUBMatrix(refflags, out inerror, out outerror);
            for (int i = 0; i < 3; i++)
            {
                lva.GetCovarianceMatrix(i, covarmatrix[i]);
            }
            string matrixtoprint = "";
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    matrixtoprint += String.Format("{0:0.0000}\t", covarmatrix[i][j]);
                matrixtoprint += "\n";
            }
            Console.WriteLine("Covariance Matrix:   Chi2 = {1}\n{0}", matrixtoprint, lva.GetChi2());

            // Save
            string pathfilename = "c:\\spicetmp\\ubsaved.xml";
            Console.WriteLine("Saving to {0}", pathfilename);
            lva.SaveUBMatrix(pathfilename);

            Thread.Sleep(1000);
            lva.LoadUBMatrix(pathfilename);
            

            while (true)
            {
                Console.Write("Input (H, K, L): ");
                string hklstring = Console.ReadLine();

                // 1. Quit loop or not? 
                if (hklstring.Contains("quit") || hklstring.Contains("exit"))
                {
                    break;
                }

                string[] terms = hklstring.Split();
                for (int i = 0; i < 3; i++)
                    inputhkl[i] = Convert.ToInt32(terms[i]);
                lva.CalculateMotorAnglesFromHKL(inputhkl, mangles);
                Console.Write("Result = ");
                for (int i = 0; i < 4; i++)
                    Console.Write("{0}\t", mangles[i]);
                Console.WriteLine();
            }

#if false
            Console.Write("Input Unitcell:  ");
            string pc = Console.ReadLine();
            Console.WriteLine("\n\nInput {0}", pc);
            lva.SetLatticeParameters(pc);
#endif

            return;
        }



        /// <summary>
        /// Determine UB matrix via interactive mode
        /// </summary>
        static public void InteractiveDetermineUB()
        {

            // 1. Get Input
            UnitCell myunitcell = new UnitCell();
            double lambda;
            Reflection[] myreflections = new Reflection[100];
            for (int i = 0; i < 100; i++)
                myreflections[i] = new Reflection();

            Console.WriteLine("Input File of Unitcell in the order of a b c alpha beta gamma (default: unitcell.txt");
            string unitcellfilename = Console.ReadLine().Split()[0];
            if (unitcellfilename.Length <= 1)
                unitcellfilename = "unitcell.txt";
            ReadUnitCellFile(unitcellfilename, myunitcell);
            // Console.WriteLine("{0}   {1}", unitcellfilename, unitcellfilename.Length);
            Console.WriteLine("Input Unit Cell: {0}", myunitcell);

            Console.WriteLine("Input Wave Length.  Default = 1.56");
            string stlambda = Console.ReadLine();
            try
            {
                lambda = Convert.ToDouble(stlambda);
                if (lambda < 1.0E-7)
                    lambda = 1.56;
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed to convert %s to double with reason %s.  Use Lambda = 1.56 instead.", stlambda, e.Message);
                lambda = 1.56;
            }
            Console.WriteLine("Lambda = {0}", lambda);

            Console.WriteLine("Input File Of 2 Reflections In Order of h, k, l, 2theta, omega, chi, phi, ");
            string instring = Console.ReadLine().Split()[0];
            Console.WriteLine(instring);
            ReadReflections(instring, 2, myreflections);
            for (int i = 0; i < 2; i ++)
                Console.WriteLine("{0}", myreflections[i]);

            // 2. Calculate UB
            UBMatrix myubmatrix = new UBMatrix(lambda);
            myubmatrix.SetDirectUnitCell(myunitcell);

            Console.WriteLine("Initial Reflects To Make UB");
            for (int n = 0; n < 2; n++)
                Console.WriteLine("{0}  {1}", n, myreflections[n]);

            int ierr;
            myubmatrix.MakeUB(myreflections[0].MillerIndex, myreflections[1].MillerIndex, myreflections[0].MotorAngles,
                myreflections[1].MotorAngles, out ierr);
            myubmatrix.Print();

            Console.WriteLine("Initial Reflects After UB Made");
            for (int n = 0; n < 2; n++)
                Console.WriteLine("{0}  {1}",n, myreflections[n]);

            bool ContinueInput = true;
            int numreflection = 2;
            while (ContinueInput == true)
            {
                Console.WriteLine("Try (H, K, L):   Invalid HKL Is For Stopping");
                string stinhkl = Console.ReadLine();
                bool tostop;
                MillerIndices trymiller = new MillerIndices();
                ParseInputHKL(stinhkl, trymiller, out tostop);

                if (tostop == true)
                    ContinueInput = false;
                else
                {                    
                    Console.WriteLine("{0}", trymiller);
                    double tryq;
                    myubmatrix.Calculate2ThetaFromHKL(trymiller, out tryq);
                    int ierr2;
                    MotorIncidentAngles newmotorangle = new MotorIncidentAngles();
                    myubmatrix.CalMotorAnglesFromMillerIndex(trymiller, newmotorangle, false, out ierr2);
                    Console.WriteLine("Bisecting:  {0}", newmotorangle);

                    Console.WriteLine("Observed Motor Angle:  Invalid = Flag For Not To Refine");
                    string strinpmotorangle = Console.ReadLine();
                    bool torefine = ParseInputIncidentAngle(strinpmotorangle, newmotorangle);
                    Console.WriteLine("Measured Angle = {0}", newmotorangle);

                    if (torefine == true)
                    {
                        numreflection += 1;
                        newmotorangle.CopyTo(myreflections[numreflection - 1].MotorAngles);
                        trymiller.CopyTo(myreflections[numreflection - 1].MillerIndex);

                        myubmatrix.RefineUBLinearLSF(myreflections, numreflection);
                        myubmatrix.Print(true);
                        Console.WriteLine("");
                        myubmatrix.PrevUB.Print(true);
                    }
                }

            } // while


        }

        /// <summary>
        /// Read unit cell file
        /// </summary>
        /// <param name="filename"></param>
        static public void ReadUnitCellFile(string filename, UnitCell outcell)
        {
            TextReader tr = new StreamReader(filename);
            string infoline = tr.ReadLine();

            string[] terms = infoline.Split();
            for (int i = 0; i < terms.Length; i++)
            {
                if (i < 3)
                    outcell.Axiles[i] = Convert.ToDouble(terms[i]);
                else
                    outcell.Angles[i - 3] = Convert.ToDouble(terms[i]);
            }
            tr.Close();

            return;
        }

        /// <summary>
        /// Read reflections from the file
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="numtoread"></param>
        /// <param name="outreflections"></param>
        static public void ReadReflections(string filename, int numtoread, Reflection[] outreflections)
        {
            TextReader tr = new StreamReader(filename);
            
            for (int i = 0; i < numtoread; i++)
            {
                string infoline = tr.ReadLine();
                string[] terms = infoline.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);             

                for (int j = 0; j < 3; j++)
                    outreflections[i].MillerIndex.V[j] = Convert.ToInt32(terms[j]);
                outreflections[i].MotorAngles.twotheta = Convert.ToDouble(terms[3]);
                outreflections[i].MotorAngles.omega    = Convert.ToDouble(terms[4])-outreflections[i].MotorAngles.twotheta/2.0;
                outreflections[i].MotorAngles.chi      = Convert.ToDouble(terms[5]);
                outreflections[i].MotorAngles.phi      = Convert.ToDouble(terms[6]);
            }

            tr.Close();

            return;

        }

        /// <summary>
        /// Read reflections (motor angles only!) from the file
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="incidentangles"></param>
        /// <param name="?"></param>
        static public void ReadReflectionMotorAngles(string filename, MotorIncidentAngles[] incidentangles, out int numreflections)
        {
            numreflections = 0;
            StreamReader sr = File.OpenText(filename);
            string infoline = null;
            while ((infoline = sr.ReadLine()) != null){
                if (infoline[0] != '#')
                {
                    // Skip comment line started with #
                    string[] terms = infoline.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

                    incidentangles[numreflections].twotheta = Convert.ToDouble(terms[0]);
                    incidentangles[numreflections].omega = Convert.ToDouble(terms[1]) - incidentangles[numreflections].twotheta / 2.0;
                    incidentangles[numreflections].chi = Convert.ToDouble(terms[2]);
                    incidentangles[numreflections].phi = Convert.ToDouble(terms[3]);

                    numreflections++;
                }
            }
            sr.Close();

            return;

        }

        /// <summary>
        /// Parse input miller indices
        /// </summary>
        /// <param name="instring"></param>
        /// <param name="outmiller"></param>
        /// <param name="tostop"></param>
        static public void ParseInputHKL(string instring, MillerIndices outmiller, out bool tostop)
        {
            string[] terms = instring.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

            if (terms.Length < 3)
            {
                tostop = true;
            }
            else
            {
                tostop = false;
                try
                {
                    for (int i = 0; i < 3; i++)
                        outmiller.V[i] = Convert.ToInt32(terms[i]);
                }
                catch (Exception e)
                {
                    Console.WriteLine("Failed to convert string to integer with reason %s. ", e.Message);
                    tostop = true;
                }
            }

            return;

        }

        /// <summary>
        /// Parse input incident angles
        /// </summary>
        /// <param name="instrang"></param>
        /// <param name="inpincidengangle"></param>
        static public bool ParseInputIncidentAngle(string instrang, MotorIncidentAngles inpincidentangle)
        {
            string[] terms = instrang.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            bool toreturn = true;
            if (terms.Length < 4)
            {
                toreturn = false;
            }
            else
            {
                inpincidentangle.twotheta = Convert.ToDouble(terms[0]);
                inpincidentangle.omega = Convert.ToDouble(terms[1]) - inpincidentangle.twotheta/2.0;
                inpincidentangle.chi = Convert.ToDouble(terms[2]);
                inpincidentangle.phi = Convert.ToDouble(terms[3]);
                toreturn = true;
            }
            return toreturn;
        } // End Function

        /// <summary>
        /// Calculate the error of a UB matrix to all reflections
        /// </summary>
        /// <param name="myubmatrix"></param>
        /// <param name="mincidentangles"></param>
        /// <param name="numinputreflections"></param>
        static public void CheckAllReflectionByUBMatrix(UBMatrix myubmatrix, MotorIncidentAngles[] mincidentangles, int numinputreflections)
        {
            double sumdiff2theta = 0;
            double sumdiffmangle = 0;

            // 1. Calculate All Error
            for (int i = 0; i < numinputreflections; i++)
            {
                double diff2theta, diffmangle;

                MillerIndices mi = new MillerIndices();
                Vector uphi = new Vector(3);
                myubmatrix.CalculateHKLFromIncidentAngle(mincidentangles[i], mi, uphi);
                myubmatrix.CheckError(mincidentangles[i], mi, 1.00, 1.00, out diff2theta, out diffmangle);

                sumdiff2theta += diff2theta * diff2theta;
                sumdiffmangle += diffmangle * diffmangle;
            }

            // 2. Output Error
            sumdiff2theta = Math.Sqrt(sumdiff2theta/numinputreflections);
            sumdiffmangle = Math.Sqrt(sumdiffmangle/numinputreflections);
            Console.WriteLine("Average Error:   2theta = {0:0.000}    motor-angle = {1:0.000}",
                sumdiff2theta, sumdiffmangle);

            return;
        } // CheckAllReflectionByUBMatrix
        
    }


}