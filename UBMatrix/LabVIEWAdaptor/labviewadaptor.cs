using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace DAS.UBMatrix
{
    public class LabVIEWAdaptor
    {
        static bool _DBOUTPUT = false;

        #region private properties

        /// <summary>
        /// UB matrix object
        /// </summary>
        UBMatrix m_ubmatrix;

        /// <summary>
        /// Number of exactly imported reflections
        /// </summary>
        int m_numInputReflections;

        /// <summary>
        /// Maximum number of reflections
        /// </summary>
        int m_maxNumInputReflections;

        /// <summary>
        /// Input motor angles from reflection file
        /// </summary>
        MotorIncidentAngles[] m_inputMotorAnglesArray;

        /// <summary>
        /// Input Miller indexes from reflection file
        /// </summary>
        MillerIndices[] m_inputHKLArray;

        /// <summary>
        /// HKL array for generated Miller indexes
        /// </summary>
        MillerIndices[] m_outputHKLArray;
          
        /// <summary>
        /// Error message
        /// </summary>
        private string m_errorMessage = "";
        #endregion

        #region contructor
        /// <summary>
        /// Extended convenient constructor
        /// </summary>
        /// <param name="wavelength"></param>
        /// <param name="maxlength"></param>
        public LabVIEWAdaptor(double wavelength, int maxlength):this(maxlength)
        {
            try
            {
                SetWaveLength(wavelength);

                return;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.LabVIEWAdaptor(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Basic constructor
        /// </summary>
        /// <param name="maxlength">maximum length of reflection list</param>
        public LabVIEWAdaptor(int maxlength)
        {
            try
            {
                // Check Input's validity
                if (maxlength <= 0)
                    throw new Exception(String.Format("Max Length of Reflections List {0} Is Not Allowed", maxlength));

                // Create UBMatrix
                m_ubmatrix = new UBMatrix();
                m_ubmatrix.WaveLength = 1.01; // Default

                // Initialize class variables
                m_maxNumInputReflections = maxlength;
                m_inputMotorAnglesArray = new MotorIncidentAngles[m_maxNumInputReflections];
                m_inputHKLArray = new MillerIndices[m_maxNumInputReflections];

                for (int r = 0; r < m_maxNumInputReflections; r++)
                {
                    m_inputMotorAnglesArray[r] = new MotorIncidentAngles();
                    m_inputHKLArray[r] = new MillerIndices();
                }

                return;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.LabVIEWAdaptor(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }
        #endregion

        #region public properties
        /// <summary>
        /// Read-only error message
        /// </summary>
        public string ErrorMessage
        {
            get { return m_errorMessage; }
        }
        /// <summary>
        /// Read-only number of imported reflections (motor angle/HKL)
        /// </summary>
        public int NumberOfImportedReflections
        {
            get { return m_numInputReflections; }
        }
        #endregion

        #region methods about lattice parameters and wave length
        /// <summary>
        /// Set wavelength to the UB matrix
        /// </summary>
        /// <param name="wavelength"></param>
        public void SetWaveLength(double wavelength)
        {
            try
            {
                // Check value
                if (wavelength <= 0.0)
                    throw new Exception(String.Format("Wavelength = {0} Is Not Allowed!", wavelength));

                // Set value to UB matrix
                m_ubmatrix.WaveLength = wavelength;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.SetWaveLength(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Set lattice parameters to the adaptor
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="alpha"></param>
        /// <param name="beta"></param>
        /// <param name="gamma"></param>
        public void SetLatticeParameters(double a, double b, double c, double alpha, double beta, double gamma)
        {
            try
            {
                // Check whether input values are valid
                if (a < 1.0E-6 || b < 1.0E-6 || c < 1.0E-6 || alpha < 1.0E-6 || beta < 1.0E-6 || gamma < 1.0E-6)
                {
                    // Return with error message
                    string errmsg = String.Format("Lattice Parameters are Not Valid:  a = {0}, b = {1}, c = {2}, alpha = {3}, beta = {4}, gamma = {5}",
                        a, b, c, alpha, beta, gamma);
                    // throw new Exception(errmsg);
                    m_errorMessage = String.Format("{0}\n{1}", errmsg, m_errorMessage);
                    Log.Write(Log.Info, "Exception in DAS.UBMatrix.LabVIEWAdaptor.SetLatticeParameters(). " + errmsg, "Error");
                }
                else
                {
                    // Build unit cell
                    string logmsg;
                    UnitCell myunitcell = new UnitCell();
                    myunitcell.a = a;
                    myunitcell.b = b;
                    myunitcell.c = c;
                    myunitcell.alpha = alpha;
                    myunitcell.beta = beta;
                    myunitcell.gamma = gamma;
                    myunitcell.CalcualteVolume(out logmsg);
                    Log.Write(Log.Info, "In DAS.UBMatrix.LabVIEWAdaptor.SetLatticeParameters(). " + logmsg, "Info");

                    // c. Set cell
                    m_ubmatrix.SetDirectUnitCell(myunitcell);
                    m_ubmatrix.StoreUnitCell();
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.SetLatticeParameters(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Set lattice parameters by input a string with order a, b, c, alpha, beta, gamma with space only
        /// </summary>
        /// <param name="ps"></param>
        public bool SetLatticeParameters(string ps)
        {
            try
            {
                // Parse string to an array of doubles
                string[] terms = ps.Split();
                double[] values = new double[6];
                int count = 0;
                for (int i = 0; i < terms.Length; i++)
                {
                    // Empty string
                    if (terms[i].Length == 0)
                        continue;

                    // Try to convert to double
                    try
                    {
                        // Convert
                        double v = Convert.ToDouble(terms[i]);
                        values[count] = v;

                        if (count == 5)
                        {
                            // Only first 6 values are needed
                            break;
                        }
                        else
                        {
                            ++count;
                        }
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("Warning.  Unable to convert {0}-th lattice parameter '{1}' to double from string due to {2}",
                            i, terms[i], e.Message);
                    }
                }

                // Check number of valid input
                if (count < 5)
                {
                    string errmsg = String.Format("Error! There are only {0} valid inputs for lattice paramter, while 6 are required. ",
                        count + 1);
                    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.SetLatticeParameters(). " + errmsg, "Error");
                    Console.WriteLine(errmsg);
                    return false;
                }

                // Set unitcell
                UnitCell uc = new UnitCell();
                for (int i = 0; i < 3; i++)
                {
                    uc.Axiles[i] = values[i];
                    uc.Angles[i] = values[i + 3];
                }
                m_ubmatrix.SetDirectUnitCell(uc);

                Console.WriteLine("Set Unit Cell (Adaptor): {0}", uc);

                return true;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.SetLatticeParameters(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Get lattice parameters.
        /// If UB matrix is refined, the lattice parameters have been refined
        /// </summary>
        /// <param name="edges"></param>
        /// <param name="angles"></param>
        public void GetLatticeParameters(double[] edges, double[] angles)
        {
            try
            {
                if (edges.Length < 3 || angles.Length < 3)
                    throw new Exception(String.Format("Input/Output Double Array's Length Not Allowed (< 3): Edges.L = {0}, Angles.L = {1}",
                        edges.Length, angles.Length));

                string logmessage = String.Format("Get Lattice (Debug):  Input = {0}, {1}", edges, angles);

                m_ubmatrix.GetDirectUnitCell(edges, angles);

                logmessage += String.Format("\t\tOutput = {0}, {1}", edges, angles);
                Log.Write(Log.Info, "In DAS.UBMatrix.LabVIEWAdaptor.GetLatticeParameters(). " + logmessage, "Info");
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GetLatticeParameters(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Get lattice parameters
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="alpha"></param>
        /// <param name="beta"></param>
        /// <param name="gamma"></param>
        public void GetLatticeParameters(out double a, out double b, out double c, out double alpha,
            out double beta, out double gamma)
        {
            try
            {
                a = m_ubmatrix.CrystallCell.a;
                b = m_ubmatrix.CrystallCell.b;
                c = m_ubmatrix.CrystallCell.c;
                alpha = m_ubmatrix.CrystallCell.alpha;
                beta = m_ubmatrix.CrystallCell.beta;
                gamma = m_ubmatrix.CrystallCell.gamma;

                return;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GetLatticeParameters(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Get wave length of the UB matrix
        /// </summary>
        /// <returns></returns>
        public double GetWaveLength()
        {
            return m_ubmatrix.WaveLength;
        }

        #endregion

        #region Methods to do UB matrix calculation
        /// <summary>
        /// Generate UB matrix from 2 reflections.  Result is stored to class variable m_UBMatrix
        /// Note: omega in the input value is of Bryan's definition (not Busing's)
        /// </summary>
        /// <param name="hkl1"></param>
        /// <param name="motorangle1">HFIR 4-circle motor position</param>
        /// <param name="hkl2"></param>
        /// <param name="motorangle2">HFIR 4-circle motor position</param>
        public void GenerateUBMatrix2Reflections(int[] hkl1, double[] motorangle1, int[] hkl2, double[] motorangle2)
        {
            try
            {
                // Check
                if (hkl1.Length < 3 || motorangle1.Length < 4 || hkl2.Length < 3 || motorangle2.Length < 4)
                {
                    // Input arrays do not have valid length
                    string errmsg = "Input hkl, motor angles' array size not correct";
                    m_errorMessage = String.Format("{0}\n{1}", errmsg, m_errorMessage);
                    Log.Write(Log.Info, "Error in DAS.UBMatrix.LabVIEWAdaptor.GenerateUBMatrix2Reflections(). " + m_errorMessage, "Error");

                    return;
                    // throw new Exception(errmsg);
                }
                else
                {
                    Console.WriteLine("Generate UB Matrix from reflection ({0}, {1}, {2}) and ({3}, {4}, {5}). ",
                        hkl1[0], hkl1[1], hkl1[2], hkl2[0], hkl2[1], hkl2[2]);
                }

                // Convert input array to MillerIndices
                MillerIndices refhkl1 = new MillerIndices();
                MillerIndices refhkl2 = new MillerIndices();
                for (int i = 0; i < 3; i++)
                {
                    refhkl1.V[i] = hkl1[i];
                    refhkl2.V[i] = hkl2[i];
                }

                // Convert input array to incident angles
                MotorIncidentAngles mangle1 = new MotorIncidentAngles();
                ConvertMotorAnglesFromBryanToBusing(motorangle1, mangle1);

                MotorIncidentAngles mangle2 = new MotorIncidentAngles();
                ConvertMotorAnglesFromBryanToBusing(motorangle2, mangle2);

                // Build UB matrix
                int ierr;
                m_ubmatrix.MakeUB(refhkl1, refhkl2, mangle1, mangle2, out ierr);
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GenerateUBMatrix2Reflections(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Pass value of UB matrix out (follow
        /// </summary>
        /// <param name="outUBmatrix"></param>
        public void GetUBMatrix(double[,] outUBmatrix)
        {
            try
            {
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        outUBmatrix[i, j] = m_ubmatrix.M[i][j];
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GetUBMatrix(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Output one row of UB matrix
        /// </summary>
        /// <param name="outSubUBmatrix"></param>
        /// <param name="row"></param>
        public void GetUBMatrix(double[] outSubUBmatrix, int row)
        {
            try
            {
                if (outSubUBmatrix.Length < 3)
                    throw new Exception(String.Format("Input Array's Length = {0} < 3", outSubUBmatrix.Length));

                try
                {
                    for (int j = 0; j < 3; j++)
                        outSubUBmatrix[j] = m_ubmatrix.M[row][j];
                }
                catch (IndexOutOfRangeException ioe)
                {
                    throw new Exception(String.Format("Index Out of Boundary!  Input Index row = {0}\n{1}\n",
                        row, ioe));
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GetUBMatrix(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Get a row value of covariant matrix
        /// </summary>
        /// <param name="rownumber"></param>
        /// <param name="rowitems"></param>
        public void GetCovarianceMatrix(int rownumber, double[] rowitems)
        {
            try
            {
                if (rownumber >= 0 && rownumber < 3)
                {
                    for (int col = 0; col < 3; col++)
                        rowitems[col] = m_ubmatrix.CovarianceMatrix.M[rownumber][col];
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GetCovarianceMatrix(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Get Chi2 of the UB matrix refined
        /// </summary>
        public double GetChi2()
        {
            return m_ubmatrix.Chi2;
        }

        #endregion  

        #region Methods to calculate motor angles

        /// <summary>
        /// Calcualte motor angles from floating point HKL
        /// Using bi-sect algorithm to calculate motor angle, i.e, omega = 0
        /// Note: returned motor angles must be Bryan/4-circle's but NOT busing's 
        /// </summary>
        /// <param name="hkl">In: [h, k, l]</param>
        /// <param name="motorangles">In/Out: [2theta, omega, chi, phi]</param>
        /// <param name="usebisec">Flag to use bisect algorithm or fixed phi. default in bisection mode</param>
        public bool CalculateMotorAnglesFromExactHKL(double[] hkl, double[] motorangles, bool usebisec = true)
        {
            try
            {
                MillerIndices mhkl = new MillerIndices();
                MotorIncidentAngles mangles = new MotorIncidentAngles();
                bool fixed_phi = !usebisec;  // bisect algorithm
                int ierror;

                // Check inputs
                if (hkl.Length < 3)
                {
                    throw new Exception(String.Format("Input Integer Array For (HKL) With Lenght = {0} < 3", hkl.Length));
                }
                if (motorangles.Length < 4)
                {
                    throw new Exception(String.Format("Input/Outpout Double Array For Motor Angles With Lenght = {0} < 4", motorangles.Length));
                }

                // Convert array to Miller Indices
                for (int i = 0; i < 3; i++)
                    mhkl.V[i] = hkl[i];

                // Calculate motor angles
                m_ubmatrix.CalMotorAnglesFromMillerIndex(mhkl, mangles, fixed_phi, out ierror);

                // Convert motor angles to double array
                ConvertMotorAnglesFromBusingToBryan(mangles, motorangles);

                return true;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CalculateMotorAnglesFromExactHKL(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }


        /// <summary>
        /// Calcualte motor angles from HKL
        /// Using bi-sect algorithm to calculate motor angle, i.e, omega = 0
        /// Note: returned motor angles must be Bryan/4-circle's but NOT busing's 
        /// </summary>
        /// <param name="hkl">In: [h, k, l]</param>
        /// <param name="motorangles">In/Out: [2theta, omega, chi, phi]</param>
        /// <param name="usebisec">Flag to use bisect algorithm or fixed phi. default in bisection mode</param>
        public bool CalculateMotorAnglesFromHKL(int[] hkl, double[] motorangles, bool usebisec=true)
        {
            try
            {
                MillerIndices mhkl = new MillerIndices();
                MotorIncidentAngles mangles = new MotorIncidentAngles();
                bool fixed_phi = !usebisec;  // bisect algorithm
                int ierror;

                // Check inputs
                if (hkl.Length < 3)
                {
                    throw new Exception(String.Format("Input Integer Array For (HKL) With Lenght = {0} < 3", hkl.Length));
                }
                if (motorangles.Length < 4)
                {
                    throw new Exception(String.Format("Input/Outpout Double Array For Motor Angles With Lenght = {0} < 4",
                            motorangles.Length));
                }

                // Convert integer array to Miller Indices
                for (int i = 0; i < 3; i++)
                    mhkl.V[i] = hkl[i];

                // Calculate motor angles
                m_ubmatrix.CalMotorAnglesFromMillerIndex(mhkl, mangles, fixed_phi, out ierror);

                // Convert motor angles to double array
                ConvertMotorAnglesFromBusingToBryan(mangles, motorangles);

                return true;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CalculateMotorAnglesFromHKL(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        #endregion

        #region Methods to handle reflections list including importing and exporting
     
        /// <summary>
        /// Import the reflections given in a text file to this object
        /// Format:
        /// (1) 1st line: # hkl angle
        ///     2+  line: h k l 2theta   omega   chi   phi
        /// (2) 1st line: # angle hkl
        ///     2+  line: 2theta omega chi phi h k l
        /// Class variables affected:
        /// - m_numImportedReflections
        /// - m_inputMotorAngleArray
        /// - m_inputHKLArray
        /// </summary>
        /// <param name="filename"></param>
        /// <returns>Int: Number of reflections imported</returns>
        public int ImportReflectionsFromFile(string filename)
        {
            int numreflections = 0;

            try
            {
                string logmessage = "";

                // Read text file
                StreamReader sr;
                try
                {
                    int filestyle;

                    sr = File.OpenText(filename);
                    try
                    {
                        Console.WriteLine("Opened reflection file {0}. ", filename);
                        string infoline = null;

                        // Format information line
                        string line0 = sr.ReadLine();

                        string[] terms = line0.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        if (terms[0] != "#")
                        {
                            throw new Exception(
                                String.Format("First Line Must Start With #.  Now It is '{0}'", line0));
                        }

                        if (terms[1] == "hkl" && terms[2] == "angle")
                        {
                            // hkl angle
                            filestyle = 1;
                        }
                        else if (terms[1] == "angle" && terms[2] == "hkl")
                        {
                            // angle hkl
                            filestyle = 2;
                        }
                        else
                        {
                            string errormessage =
                                String.Format("Line 0 of reflection file {0} must have 'angle' and 'hkl' (case sensitive).  Current line 0 is \"{1}\". ",
                                filename, line0);
                            throw new Exception(errormessage);
                        }

                        // Parse hkl and motor angles in the rest lines
                        numreflections = 0;
                        int refindex = m_numInputReflections;
                        double[] tempmotorpos = new double[4];
                        while ((infoline = sr.ReadLine()) != null)
                        {
                            // Array boundary check
                            if (refindex >= m_maxNumInputReflections)
                            {
                                // Break loop if out of maximum number of allowed reflections
                                throw new Exception(String.Format("Unable to import {0}-th input reflections and above due to exceeding limit of {1} reflections.",
                                    numreflections + 1, m_maxNumInputReflections));
                            }

                            if (infoline.Length > 0 && infoline[0] != '#')
                            {
                                // Skip comment line started with #
                                terms = infoline.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

                                if (terms.Length >= 7)
                                {
                                    try
                                    {
                                        if (filestyle == 1)
                                        {
                                            // Allow nonintegral hkl in reflection file. rdg; 1/8/14
                                            //m_inputHKLArray[refindex].H = Convert.ToInt32(terms[0]);
                                            //m_inputHKLArray[refindex].K = Convert.ToInt32(terms[1]);
                                            //m_inputHKLArray[refindex].L = Convert.ToInt32(terms[2]);
                                            m_inputHKLArray[refindex].H = Convert.ToDouble(terms[0]);
                                            m_inputHKLArray[refindex].K = Convert.ToDouble(terms[1]);
                                            m_inputHKLArray[refindex].L = Convert.ToDouble(terms[2]);

                                            // Import and Export functions keep motor angles in Bryan coordinate system because they are only meant
                                            // for passing data between SPICE and reflection files. rdg; 1/17/14
                                            m_inputMotorAnglesArray[refindex].twotheta = Convert.ToDouble(terms[3]);
                                            m_inputMotorAnglesArray[refindex].omega = Convert.ToDouble(terms[4]);
                                            m_inputMotorAnglesArray[refindex].chi = Convert.ToDouble(terms[5]);
                                            m_inputMotorAnglesArray[refindex].phi = Convert.ToDouble(terms[6]);
                                        }
                                        else if (filestyle == 2)
                                        {
                                            // Import and Export functions keep motor angles in Bryan coordinate system because they are only meant
                                            // for passing data between SPICE and reflection files. rdg; 1/17/14
                                            m_inputMotorAnglesArray[refindex].twotheta = Convert.ToDouble(terms[0]);
                                            m_inputMotorAnglesArray[refindex].omega = Convert.ToDouble(terms[1]);
                                            m_inputMotorAnglesArray[refindex].chi = Convert.ToDouble(terms[2]);
                                            m_inputMotorAnglesArray[refindex].phi = Convert.ToDouble(terms[3]);

                                            // Allow nonintegral hkl in reflection file. rdg; 1/8/14
                                            //m_inputHKLArray[refindex].H = Convert.ToInt32(terms[4]);
                                            //m_inputHKLArray[refindex].K = Convert.ToInt32(terms[5]);
                                            //m_inputHKLArray[refindex].L = Convert.ToInt32(terms[6]);
                                            m_inputHKLArray[refindex].H = Convert.ToDouble(terms[4]);
                                            m_inputHKLArray[refindex].K = Convert.ToDouble(terms[5]);
                                            m_inputHKLArray[refindex].L = Convert.ToDouble(terms[6]);
                                        }

                                        ++numreflections;
                                        ++refindex;
                                    }
                                    catch (FormatException fe)
                                    {
                                        Console.WriteLine("Line {0}: '{1}' Format Error due to {2}!",
                                            numreflections, infoline, fe.Message);
                                    }
                                }
                                else
                                {
                                    logmessage += String.Format("Input Data {0} Cannot Be Understood\n", infoline);
                                }
                            }
                        }
                    }
                    catch
                    {
                        throw;
                    }
                    finally
                    {
                        sr.Close();
                    }
                }
                catch (Exception e)
                {
                    throw new Exception(String.Format("Fail To Import File \"{0}\" due to error: \n{1}\n", filename, e.Message));
                }

                m_numInputReflections += numreflections;

                if (logmessage != "")
                {
                    Log.Write(Log.Info, "In DAS.UBMatrix.LabVIEWAdaptor.ImportReflectionsFromFile(). " + logmessage, "Info");
                    Console.WriteLine(logmessage);
                }

                return numreflections;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.ImportReflectionsFromFile(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Reads HKL Intensity File.
        /// </summary>
        /// <param name="filename"></param>
        /// <returns></returns>
        public double[,] ReadHKLIntensityFile(string filename)
        {
            try
            {
                List<double[]> hkl_intensity_values = new List<double[]>();

                string logmessage = "";

                // Read text file
                StreamReader sr;
                try
                {
                    sr = File.OpenText(filename);
                    try
                    {
                        Console.WriteLine("Opened hkl-intensity file {0}. ", filename);
                        string infoline = null;

                        // Parse hkl and intensity in the rest lines
                        double[] temp_hkl_intensity = new double[4];
                        while ((infoline = sr.ReadLine()) != null)
                        {
                            if (infoline.Length > 0 && infoline[0] != '#')
                            {
                                // Skip comment line started with #
                                string[] terms = infoline.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

                                if (terms.Length >= 4)
                                {
                                    try
                                    {
                                        double[] hkl_intensity = new double[4];

                                        // h
                                        hkl_intensity[0] = Convert.ToDouble(terms[0]);
                                        // k
                                        hkl_intensity[1] = Convert.ToDouble(terms[1]);
                                        // l
                                        hkl_intensity[2] = Convert.ToDouble(terms[2]);
                                        // intensity
                                        hkl_intensity[3] = Convert.ToDouble(terms[3]);

                                        hkl_intensity_values.Add(hkl_intensity);
                                    }
                                    catch (FormatException fe)
                                    {
                                        Console.WriteLine("Line {0}: '{1}' Format Error due to {2}!",
                                            hkl_intensity_values.Count, infoline, fe.Message);
                                    }
                                }
                                else
                                {
                                    logmessage += String.Format("Input Data {0} Cannot Be Understood\n", infoline);
                                }
                            }
                        }
                    }
                    catch
                    {
                        throw;
                    }
                    finally
                    {
                        sr.Close();
                    }
                }
                catch (Exception e)
                {
                    throw new Exception(String.Format("Fail To Import File \"{0}\" due to error: \n{1}\n", filename, e.Message));
                }

                // Convert to 2D array for return.
                double[,] rtn_value = new double[hkl_intensity_values.Count, 4];

                for (int i = 0; i < hkl_intensity_values.Count; ++i)
                {
                    rtn_value[i, 0] = hkl_intensity_values[i][0];
                    rtn_value[i, 1] = hkl_intensity_values[i][1];
                    rtn_value[i, 2] = hkl_intensity_values[i][2];
                    rtn_value[i, 3] = hkl_intensity_values[i][3];
                }

                return rtn_value;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.ReadHKLIntensityFile(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// ImportReflections()
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="motor_angles"></param>
        /// <param name="hkl"></param>
        /// <returns></returns>
        public Reflection[] ImportReflections(string filename)
        {
            try
            {
                ClearReflections();
                int numreflections = ImportReflectionsFromFile(filename);

                return GetAllCurrentReflections();
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.ImportReflections(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }
  
        /// <summary>
        /// Append a reflection to imported motor angles and hkl list
        /// </summary>
        /// <param name="motorangles">array containing motor angles in order of 2theta/omega/chi/phi</param>
        /// <param name="message">error message if return false</param>
        /// <returns>boolean to indicate whether appending is successful </returns>
        public bool AppendReflection(double[] motorangles, out string message)
        {
            try
            {
                int refindex = m_numInputReflections;
                if (refindex >= m_maxNumInputReflections)
                {
                    // Exceeding limit of number of reflections
                    throw new Exception(String.Format("Input reflections have already reached maximum number limit."));
                }

#if false
            m_inputMotorAnglesArray[refindex].twotheta = motorangles[0];
            m_inputMotorAnglesArray[refindex].omega = motorangles[1];
            m_inputMotorAnglesArray[refindex].chi = motorangles[2];
            m_inputMotorAnglesArray[refindex].phi = motorangles[3];
#else
                ConvertMotorAnglesFromBryanToBusing(motorangles, m_inputMotorAnglesArray[refindex]);
#endif

                m_inputHKLArray[refindex].H = 0;
                m_inputHKLArray[refindex].K = 0;
                m_inputHKLArray[refindex].L = 0;

                ++m_numInputReflections;

                message = "";

                return true;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.AppendReflection(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                message = log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// AppendReflection()
        /// </summary>
        /// <param name="motorangles"></param>
        /// <param name="hkl"></param>
        public void AppendReflection(double[] motorangles, double[] hkl)
        {
            try
            {
                int refindex = m_numInputReflections;

#if false
            m_inputMotorAnglesArray[refindex].twotheta = motorangles[0];
            m_inputMotorAnglesArray[refindex].omega = motorangles[1];
            m_inputMotorAnglesArray[refindex].chi = motorangles[2];
            m_inputMotorAnglesArray[refindex].phi = motorangles[3];
#else
                ConvertMotorAnglesFromBryanToBusing(motorangles, m_inputMotorAnglesArray[refindex]);
#endif

                m_inputHKLArray[refindex].H = hkl[0];
                m_inputHKLArray[refindex].K = hkl[1];
                m_inputHKLArray[refindex].L = hkl[2];

                ++m_numInputReflections;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.AppendReflection(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }


        /// <summary>
        /// Get the number of imported reflections
        /// </summary>
        /// <returns></returns>
        public int GetNumberOfReflections()
        {
            return m_numInputReflections;
        }

        /// <summary>
        /// Get one reflection
        /// </summary>
        /// <param name="index"></param>
        /// <param name="motorangles">out: motor angles [2theta, omega, chi, phi]</param>
        public void GetImportedReflection(int index, double[] motorangles, int[] millerindex)
        {
            try
            {
                if (index > m_inputMotorAnglesArray.Length)
                {
                    throw new Exception(
                        String.Format("Out of boundary  input index = {0}, Max Length = {1}",
                        index, m_inputMotorAnglesArray.Length));
                }

                // Motor angles
#if false
            reflection[0] = m_inputMotorAnglesArray[index].twotheta;
            reflection[1] = m_inputMotorAnglesArray[index].omega;
            reflection[2] = m_inputMotorAnglesArray[index].chi;
            reflection[3] = m_inputMotorAnglesArray[index].phi;
#else
                ConvertMotorAnglesFromBusingToBryan(m_inputMotorAnglesArray[index], motorangles);
#endif

                // HKL 
                for (int i = 0; i < 3; i++)
                    millerindex[i] = Convert.ToInt32(m_inputHKLArray[index].V[i]);

                return;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GetImportedReflection(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Export the reflections to file as 
        /// Output format: 
        ///    1st line: # angle hkl
        ///    2+  line: 2theta omega chi phi h k l
        /// </summary>
        /// <param name="filename"></param>
        public void ExportReflections(string filename)
        {
            try
            {
                // Buffer string
                StringBuilder sb = new StringBuilder();
                sb.AppendLine("# angle hkl\n");

                double[] tempmotors = new double[4];

                for (int reflection_index = 0; reflection_index < m_numInputReflections; ++reflection_index)
                {
                    // Import and Export functions keep store motor angles in Bryan coordinate system because they are only meant
                    // for passing data between SPICE and reflection files.
                    // We convert here because the Append() method, which was used
                    // to move the reflections from SPICE to this DLL converts
                    // from Bryan to Busing. rdg; 1/17/14
                    ConvertMotorAnglesFromBusingToBryan(m_inputMotorAnglesArray[reflection_index], tempmotors);

                    string tempstr = String.Format("{0:F4} \t{1:F4} \t{2:F4} \t{3:F4} \t{4:F4} \t{5:F4} \t{6:F4}",
                        tempmotors[0], tempmotors[1], tempmotors[2], tempmotors[3],
                        m_inputHKLArray[reflection_index].H, m_inputHKLArray[reflection_index].K, m_inputHKLArray[reflection_index].L);
                    sb.AppendLine(tempstr);
                }

                // Export to file
                using (StreamWriter outfile = new StreamWriter(filename))
                {
                    outfile.Write(sb.ToString());
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.ExportReflections(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Clear all the imported reflections
        /// </summary>
        public void ClearReflections()
        {
            m_numInputReflections = 0;
        }

        /// <summary>
        /// Remove the reflections that are not wanted
        /// </summary>
        /// <param name="selectflags"></param>
        public void CleanReflections(bool[] selectflags)
        {
            try
            {
                // Determine max
                int numref2proc = Math.Min(m_numInputReflections, selectflags.Length);

                // Loop around
                int index = 0;
                for (int i = 0; i < numref2proc; ++i)
                {
                    // Skip 
                    if (!selectflags[i])
                    {
                        // Unselected
                        continue;
                    }
                    else if (index == i)
                    {
                        // Skip same: no need to rewrite
                        ++index;
                        continue;
                    }

                    // Move
                    m_inputMotorAnglesArray[i].CopyTo(m_inputMotorAnglesArray[index]);
                    m_inputHKLArray[i].CopyTo(m_inputHKLArray[index]);

                    ++index;
                }

                // Reset num reflection
                m_numInputReflections = index;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CleanReflections(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        #endregion

        #region private methods about refining UB matrix
        /// <summary>
        /// Create a list of Reflections from user-specified refinement flag of each 
        /// </summary>
        /// <returns></returns>
        private Reflection[] CreateReflectionListToRefine(bool[] reflectionincludedflags)
        {
            try
            {
                // Determine the number of reflections for refining UB matrix'
                int numobservedreflection = Math.Min(m_numInputReflections, reflectionincludedflags.Length);

                // Find out how many reflections to be involved in refinement
                int numreflections = 0;
                for (int i = 0; i < numobservedreflection; ++i)
                {
                    if (reflectionincludedflags[i])
                    {
                        ++numreflections;
                    }
                }

                // Build the array of reflections (HKL+motor anges)
                Reflection[] mreflectionslist = new Reflection[numreflections];
                int refindex = 0;
                for (int i = 0; i < numobservedreflection; ++i)
                {
                    // Console.WriteLine("Input reflection {0} ", m_inputHKLArray[i]);
                    if (reflectionincludedflags[i])
                    {
                        // Check range of refindex
                        if (refindex >= numreflections)
                        {
                            m_errorMessage = "Logic Error!  Reflection index is larger than number of reflections.";
                            throw new Exception(m_errorMessage);
                        }

                        // To refine with this refection                    
#if false
                    MotorIncidentAngles tempangles = new MotorIncidentAngles();
                    tempangles.twotheta = m_inputMotorAnglesArray[i].twotheta;
                    tempangles.omega = m_inputMotorAnglesArray[i].omega; // -0.5 * tempangles.twotheta;
                    tempangles.chi = m_inputMotorAnglesArray[i].chi;
                    tempangles.phi = m_inputMotorAnglesArray[i].phi;
#endif

                        mreflectionslist[refindex] = new Reflection();
                        mreflectionslist[refindex].MotorAngles = m_inputMotorAnglesArray[i];
                        mreflectionslist[refindex].MillerIndex = m_inputHKLArray[i];

                        ++refindex;
                    }
                }

                return mreflectionslist;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CreateReflectionListToRefine(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// GetAllCurrentReflections()
        /// </summary>
        /// <returns></returns>
        private Reflection[] GetAllCurrentReflections()
        {
            try
            {
                // Build the array of reflections (HKL+motor anges)
                Reflection[] reflection_list = new Reflection[m_numInputReflections];
                for (int ref_index = 0; ref_index < m_numInputReflections; ++ref_index)
                {
                    // To refine with this refection
                    reflection_list[ref_index] = new Reflection();
                    m_inputMotorAnglesArray[ref_index].CopyTo(reflection_list[ref_index].MotorAngles);
                    m_inputHKLArray[ref_index].CopyTo(reflection_list[ref_index].MillerIndex);
                }

                return reflection_list;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GetAllCurrentReflections(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }


        #endregion

        #region public methods about refining UB Matrix
        /// <summary>
        /// Refine UB matrix
        /// </summary>
        /// <param name="reflectionincludedflags">array of boolean to indicate whether the reflection is used in refining</param>
        public bool RefineUBMatrix(bool[] reflectionincludedflags, out double inerror, out double outerror)
        {
            try
            {
                string logmessage = "LV Adaptor: Starting to Refine UB Matrix\n";
                logmessage += String.Format("ReflectionIncludedFlag# = {0}\n", reflectionincludedflags.Length);

                Reflection[] mreflectionslist = CreateReflectionListToRefine(reflectionincludedflags);
                int numreflections = mreflectionslist.Length;

                inerror = m_ubmatrix.CalculateUBMatrixError(mreflectionslist);

                // Logging
                logmessage += string.Format("Total {0} Reflections/Motor Angles Go To Least Square Fitting", numreflections);
                Console.WriteLine("DB: Total {0} Reflections/Motor Angles Go To Least Square Fitting", numreflections);
                Log.Write(Log.Info, "In DAS.UBMatrix.LabVIEWAdaptor.RefineUBMatrix(). " + logmessage, "Info");

                // Refine
                bool refinegood = m_ubmatrix.RefineUBLinearLSF(mreflectionslist, numreflections);
                if (!refinegood)
                {
                    outerror = 1.234E10;
                    return false;
                }

                // Error and return
                outerror = m_ubmatrix.CalculateUBMatrixError(mreflectionslist);

                return true;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.RefineUBMatrix(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Refine UB matrix
        /// </summary>
        /// <param name="inerror"></param>
        /// <param name="outerror"></param>
        /// <returns></returns>
        public bool RefineUBMatrix(out double inerror, out double outerror)
        {
            try
            {
                bool[] use_reflection = new bool[m_numInputReflections];

                for (int i = 0; i < m_numInputReflections; ++i) use_reflection[i] = true;
                return RefineUBMatrix(use_reflection, out inerror, out outerror);
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.RefineUBMatrix(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }
        
        /// <summary>
        /// Refine UB matrix considering uncertainty of motor offset
        /// Algorithm:  do a grid-search on 2theta, omega and chi; return the best refinement result 
        ///             i.e., one solution with smallest error
        /// </summary>
        /// <param name="reflectionincludedflags"></param>
        /// <param name="delta_2theta"></param>
        /// <param name="step_2theta"></param>
        /// <param name="delta_omega"></param>
        /// <param name="step_omega"></param>
        /// <param name="delta_chi"></param>
        /// <param name="step_chi"></param>
        /// <param name="MotorPosOffsets">Output: motor position offset yielding best UB matrix</param>
        /// <returns></returns>
        public bool RefineUBMatrixUnknowMotorAngleOffset(bool[] reflectionincludedflags,  
            double delta_2theta, double step_2theta,
            double delta_omega, double step_omega, double delta_chi, double step_chi,
            double[] MotorPosOffsets, out double outerror)
        {
            try
            {
                outerror = 0.0;

                // Check
                if (MotorPosOffsets.Length < 3)
                {
                    throw new Exception ("Input MotorPosOffset array's length is shorter than 3. ");
                }

                // Construct motor angles range
                MotorIncidentAngles motoranglerange = new MotorIncidentAngles();
                motoranglerange.twotheta = delta_2theta;
                motoranglerange.omega = delta_omega;
                motoranglerange.chi = delta_chi;
                motoranglerange.phi = 0;

                // Construct motor angle step
                MotorIncidentAngles motoranglestep = new MotorIncidentAngles();
                motoranglestep.twotheta = step_2theta;
                motoranglestep.omega = step_omega;
                motoranglestep.chi = step_chi;
                motoranglestep.phi = 1.0;

                // Reflections
                Reflection[] refineubreflections = CreateReflectionListToRefine(reflectionincludedflags);

                // Call 
                MotorIncidentAngles bestmotorangleoffset = new MotorIncidentAngles();
                bool refinestatus = m_ubmatrix.GridRefineUBMatrix(refineubreflections, motoranglerange, motoranglestep,
                    bestmotorangleoffset);

                // Set up return
                // Convert motor angles to double array
                double[] tempmotors = new double[4];
                ConvertMotorAnglesFromBusingToBryan(bestmotorangleoffset, tempmotors);
                MotorPosOffsets[0] = tempmotors[0];
                MotorPosOffsets[1] = tempmotors[1];
                MotorPosOffsets[2] = tempmotors[2];
                if (MotorPosOffsets.Length >= 4)
                    MotorPosOffsets[3] = 0;

                // Error and return
                // Need to calculate error based on shifted reflections.
                // rdg; 1/28/14
                Reflection[] vecShiftedReflects = m_ubmatrix.GenerateReflectionOffsetMotorPositions(refineubreflections, bestmotorangleoffset.twotheta,
                    bestmotorangleoffset.omega, bestmotorangleoffset.chi, bestmotorangleoffset.phi);

                outerror = m_ubmatrix.CalculateUBMatrixError(vecShiftedReflects);

                return refinestatus;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.RefineUBMatrixUnknowMotorAngleOffset(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Refine UB matrix considering uncertainty of motor offset
        /// Algorithm:  do a grid-search on 2theta, omega and chi; return the best refinement result 
        ///             i.e., one solution with smallest error
        /// </summary>
        /// <param name="delta_2theta"></param>
        /// <param name="step_2theta"></param>
        /// <param name="delta_omega"></param>
        /// <param name="step_omega"></param>
        /// <param name="delta_chi"></param>
        /// <param name="step_chi"></param>
        /// <param name="MotorPosOffsets"></param>
        /// <param name="outerror"></param>
        /// <returns></returns>
        public bool RefineUBMatrixUnknowMotorAngleOffset(double delta_2theta, double step_2theta,
            double delta_omega, double step_omega, double delta_chi, double step_chi,
            double[] MotorPosOffsets, out double outerror)
        {
            try
            {
                bool[] use_reflection = new bool[m_numInputReflections];

                for (int i = 0; i < m_numInputReflections; ++i) use_reflection[i] = true;

                return RefineUBMatrixUnknowMotorAngleOffset(use_reflection, delta_2theta, step_2theta, delta_omega, step_omega,
                                                             delta_chi, step_chi, MotorPosOffsets, out outerror);
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.RefineUBMatrixUnknowMotorAngleOffset(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }


        #endregion

        #region methods to calculate motor angles or Miller index
        /// <summary>
        /// Calculate 2theta according to Miller index (HKL)
        /// </summary>
        /// <param name="hkl"></param>
        /// <returns>2theta</returns>
        public double Calculate2ThetaFromMillerIndex(int[] hkl)
        {
            try
            {
                MillerIndices mi = new MillerIndices();
                for (int i = 0; i < 3; i++)
                    mi.V[i] = hkl[i];

                double d;
                double twotheta = m_ubmatrix.Calculate2ThetaFromHKL(mi, out d);

                return twotheta;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.Calculate2ThetaFromMillerIndex(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Calculate the angle between 2 miller indices
        /// </summary>
        /// <param name="hkl1"></param>
        /// <param name="hkl2"></param>
        /// <returns></returns>
        public double CalculateAngleBW2MillerIndices(int[] hkl1, int[] hkl2)
        {
            try
            {
                MillerIndices v1 = new MillerIndices();
                MillerIndices v2 = new MillerIndices();
                for (int i = 0; i < 3; i++)
                {
                    v1.V[i] = hkl1[i];
                    v2.V[i] = hkl2[i];
                }

                double angle = m_ubmatrix.CalculateAngle2MillerIndicies(v1, v2);

                return angle;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CalculateAngleBW2MillerIndices(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Calculate the angle between u_phis calculated from 2 motor angles 
        /// (i.e., incident angles)
        /// </summary>
        /// <param name="mangle1"></param>
        /// <param name="mangle2"></param>
        /// <returns></returns>
        public double CalculateAngleBW2MotorAngles(double[] mangle1, double[] mangle2)
        {
            try
            {
                MotorIncidentAngles[] iangles = new MotorIncidentAngles[2];
                Vector[] uphis = new Vector[2];

                for (int i = 0; i < 2; i++)
                {
                    iangles[i] = new MotorIncidentAngles();
                    double[] mangle;
                    if (i == 0)
                        mangle = mangle1;
                    else
                        mangle = mangle2;
#if false
                iangles[i].twotheta = mangle[0];
                iangles[i].omega = mangle[1] - 0.5*mangle[0];  // Conversion from HFIR 4cirlcle to Busing 4circle
                iangles[i].chi = mangle[2];
                iangles[i].phi = mangle[3];
#else
                    ConvertMotorAnglesFromBryanToBusing(mangle, iangles[i]);
#endif

                    uphis[i] = new Vector(3);
                    MotorIncidentAngles.UnitVectorFromAngle(iangles[i], uphis[i]);
                }

                double dotp = Vector.DotProduct(uphis[0], uphis[1]);
                if (dotp > 1.0)
                    dotp = 1.0;

                double diffangle = DasMath.AcosD(dotp);

                return diffangle;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CalculateAngleBW2MotorAngles(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Calculate Miller index from motor angles by UB matrix
        /// Input motor angle's omega will be re-defined to Busing paper's definition
        /// </summary>
        /// <param name="hkl">Output HKL in integer</param>
        /// <param name="ohkl"> Output HKL in double</param>
        /// <param name="motorangles">Input motor angle w/ definition of 4-circle</param>
        public void CalculateMillerIndexFromMotorAngle(int[] hkl, double[] ohkl, double[] motorangles)
        {
            try
            {
                // Check input arrays' lengthes
                if (hkl.Length < 3)
                {
                    throw new Exception(String.Format("Input Array (HKL) Length = {0} < 3.  Motor Angle = {1}", hkl.Length, motorangles));
                }

                // Convert double arrays to motor angle and take care of Brian and Busing's difference in definition
                MotorIncidentAngles angles = new MotorIncidentAngles();
#if false
            angles.twotheta = motorangles[0];
            angles.omega = motorangles[1] - angles.twotheta * 0.5;
            angles.chi = motorangles[2];
            angles.phi = motorangles[3];
#else
                ConvertMotorAnglesFromBryanToBusing(motorangles, angles);
#endif

                // Calculate miller index from motor angles and 
                MillerIndices mhkl = new MillerIndices();
                Vector uphi = new Vector(3);

                m_ubmatrix.CalculateHKLFromIncidentAngle(angles, mhkl, uphi);

                // Output
                for (int i = 0; i < 3; i++)
                {
                    hkl[i] = Convert.ToInt32(mhkl.V[i]);
                    ohkl[i] = mhkl.OriginalHKL[i];
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CalculateMillerIndexFromMotorAngle(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        #endregion

        #region method about I/O

        /// <summary>
        /// Save UB matrix, Unit Cell and wavelength information to a file 
        /// in XML format
        /// </summary>
        /// <param name="pathfilename"></param>
        public void SaveUBMatrix(string pathfilename)
        {
            try
            {
                string logmsg = String.Format("Save To File {0}", pathfilename);

                try
                {
                    m_ubmatrix.SaveToXMLFile(pathfilename);
                    Log.Write(Log.Info, "In DAS.UBMatrix.LabVIEWAdaptor.SaveUBMatrix(). " + logmsg, "Info");
                }
                catch (System.ArgumentException ae)
                {
                    throw new Exception(logmsg + String.Format(" Failed. Invalid Path: {0}\n", ae));
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CalculateMillerIndexFromMotorAngle(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Load UB matrix, unit cell and wave length from a file
        /// in XML format
        /// </summary>
        /// <param name="pathfilename"></param>
        public void LoadUBMatrix(string pathfilename)
        {
            try
            {
                try
                {
                    m_ubmatrix.LoadFromXMLFile(pathfilename);
                }
                catch (System.ArgumentException ae)
                {
                    throw new Exception(String.Format("Invalid Path {0}\n({1})\n", pathfilename, ae));
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.LoadUBMatrix(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Convert motor angles to Busing's definition from Bryan's
        /// Motor angle should be 2theta, omega, chi and phi
        /// </summary>
        /// <param name="motorangle"></param>
        /// <param name="mangle"></param>
        private void ConvertMotorAnglesFromBryanToBusing(double[] motorangle, MotorIncidentAngles mangle)
        {
            try
            {
                if (motorangle.Length != 4)
                    throw new Exception("Input motor angle array has a wrong array size.");

                mangle.twotheta = motorangle[0];
                mangle.omega = motorangle[1] - mangle.twotheta * 0.5;
                mangle.chi = motorangle[2];
                mangle.phi = motorangle[3];

            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.ConvertMotorAnglesFromBryanToBusing(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Convert motor angles to Bryan's definition from Busing's
        /// Motor angle should be 2theta, omega, chi and phi
        /// </summary>
        /// <param name="mangles"></param>
        /// <param name="motorangles"></param>
        private void ConvertMotorAnglesFromBusingToBryan(MotorIncidentAngles mangles, double[] motorangles)
        {
            try
            {
                if (motorangles.Length != 4)
                    throw new Exception("Input motor angle array has a wrong array size.");

                motorangles[0] = mangles.twotheta;
                motorangles[1] = mangles.omega + mangles.twotheta * 0.5;  // Convert to Brian's omega
                motorangles[2] = mangles.chi;
                motorangles[3] = mangles.phi;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.ConvertMotorAnglesFromBusingToBryan(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        #endregion

        #region (public/private) methods about generate HKL list

        /// <summary>
        /// Check whether an input condition is valid boolean expression for h, k, and l
        /// Its return should be a boolean string "ture / false"
        /// Formula: all in lower case
        ///          variable: h, k, l
        ///          a standard C# math expression
        /// </summary>
        /// <param name="formula"></param>
        /// <returns>valid or invalid</returns>
        public bool CheckConditionValidity(string condition, out string error_message)
        {
            try
            {
                error_message = "";

                // Force to lower case
                string formula = condition.Replace(" ", "").ToLower();

                if (formula.Length > 0)
                {
                    string[] variablenames = new string[3] { "h", "k", "l" };
                    int[] variablevalues = new int[3] { 1, 1, 1 };

                    string valuestr;
                    bool valid = ExpressioneValuator.Eval(formula, variablenames, variablevalues, out valuestr);

                    if (!valid)
                    {
                        error_message = "Invalid condition: " + valuestr;
                        Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.CheckConditionValidity(). " + error_message, "Error");
                        return false;
                    }
                    else if (valuestr != "true" && valuestr != "false")
                    {
                        error_message = "Condition " + condition + " does not yield a boolean result, but " + valuestr;
                        Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.CheckConditionValidity(). " + error_message, "Error");
                        return false;
                    }
                }

                return true;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CheckConditionValidity(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// Generate a list of reflections and their motor positions within by motor position limits
        /// Note: this is an 'adaptor-only' method as a tool 
        /// </summary>
        /// <param name="minhklindex"></param>
        /// <param name="maxhklindex"></param>
        /// <param name="minmotorpos">order: 2theta, omega, phi, chi</param>
        /// <param name="maxmotorpos">order: 2theta, omega, phi, chi</param>
        /// <param name="hklshift">shift on </param>
        /// <param name="condition"></param>
        public bool GenerateReflectionList(int[] minhklindex, int[] maxhklindex, double[] minmotorpos, double[] maxmotorpos, double[] hklshift, string condition,
            out string finalist, out string processing_info)
        {
            try
            {
                finalist = "";
                processing_info = "";

                // Check array validity 
                if (minhklindex.Length != 3 || maxhklindex.Length != 3)
                {
                    processing_info = "Error: Unable to process due to invalid number of HKL.";
                    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateReflectionList(). " + processing_info, "Error");
                    return false;
                }

                if (minmotorpos.Length != 4 || maxmotorpos.Length != 4)
                {
                    processing_info = "Error: Unable to process due to invalid number of motor positions.";
                    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateReflectionList(). " + processing_info, "Error");
                    return false;
                }

                if (hklshift.Length != 3)
                {
                    processing_info = "Error: Unable to process due to invalid number of HKL shift.";
                    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateReflectionList(). " + processing_info, "Error");
                    return false;
                }

                if (_DBOUTPUT)
                {
                    Console.WriteLine("Generate a list of reflections with range H = ({0}, {1}), K = ({2}, {3}), L = ({4}, {5}). ",
                        minhklindex[0], maxhklindex[0], minhklindex[1], maxhklindex[1], minhklindex[2], maxhklindex[2]);
                    Console.WriteLine("  Motor position is limited as 2theta = ({0}, {1}), Omega = ({2}, {3}), Phi = ({4}, {5}), Chi = ({6}, {7}).",
                        minmotorpos[0], maxmotorpos[0], minmotorpos[1], maxmotorpos[1],
                        minmotorpos[2], maxmotorpos[2], minmotorpos[3], maxmotorpos[3]);
                }

                // NOTE: No longer restricting reflection list to nonnegative hkl. rdg; 2/4/14
                //// Check input value in range     
                //if (minhklindex[0] < 0 || minhklindex[1] < 0 || minhklindex[2] < 0)
                //{
                //    // Output error message and return
                //    processing_info = "Error: Input HKL limits is unphysical.";
                //    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateReflectionList(). " + processing_info, "Error");
                //    return false;
                //}

                // Process condition
                bool checkcondition;
                condition = condition.Replace(" ", "");
                condition = condition.ToLower();
                ExpressioneValuator conditionevaluator = null;
                if (condition.Length > 0)
                {
                    // Check condition
                    conditionevaluator = new ExpressioneValuator();
                    string[] parameternames = new string[3] { "h", "k", "l" };
                    conditionevaluator.SetupIntExpression(condition, parameternames);
                    string errormessage;
                    bool expvalid = conditionevaluator.IsReturnBoolean(out errormessage);
                    if (!expvalid)
                    {
                        processing_info = String.Format("Error: Unable to evaluate boolean expression {0} caused by {1}", condition, errormessage);
                        Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateReflectionList(). " + processing_info, "Error");
                        Console.WriteLine(processing_info);
                        return false;
                    }
                    checkcondition = true;
                }
                else
                {
                    // Empty condition. Allow.
                    checkcondition = false;
                }

                int numh = (int)maxhklindex[0] - (int)minhklindex[0] + 1;
                int numk = (int)maxhklindex[1] - (int)minhklindex[1] + 1;
                int numl = (int)maxhklindex[2] - (int)minhklindex[2] + 1;

                if (numh < 1 || numk < 1 || numl < 1)
                {
                    // Output error message and return
                    processing_info = "Error: Allowed number of H/K/L is zero.";
                    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateReflectionList(). " + processing_info, "Error");
                    return false;
                }

                // Data strcutre
                List<MillerIndices> list_genHKL = new List<MillerIndices>();

                MillerIndices hkl = new MillerIndices();
                for (int h = (int)minhklindex[0]; h <= (int)maxhklindex[0]; ++h)
                {
                    hkl.H = (double)h;
                    for (int k = (int)minhklindex[1]; k <= (int)maxhklindex[1]; ++k)
                    {
                        hkl.K = (double)k;
                        for (int l = (int)minhklindex[2]; l <= (int)maxhklindex[2]; ++l)
                        {
                            hkl.L = (double)l;

                            // Filter by condition
                            if (h + k + l <= 0)
                                continue;
                            if (checkcondition)
                            {
                                bool meetcondition = CheckMillerIndexWithCondition(conditionevaluator, hkl);
                                if (!meetcondition)
                                {
                                    // HKL does not meet condition
                                    // Console.WriteLine("HKL = ({0}, {1}, {2}) does not meet reflection condition {3}",
                                    //    hkl.H, hkl.K, hkl.L, condition);
                                    processing_info += String.Format("HKL = ({0}, {1}, {2}) does not meet reflection condition {3}\n",
                                        hkl.H, hkl.K, hkl.L, condition);
                                    continue;
                                }
                            }

                            // Shift motor angles
                            MillerIndices shiftedhkl = new MillerIndices();
                            for (int im = 0; im < 3; ++im)
                                shiftedhkl.V[im] = hklshift[im] + hkl.V[im];

                            // Calculate motor angles
                            int ierror;
                            MotorIncidentAngles mangles = new MotorIncidentAngles();
                            m_ubmatrix.CalMotorAnglesFromMillerIndex(shiftedhkl, mangles, false, out ierror);

                            double[] motor_angles = new double[4];
                            ConvertMotorAnglesFromBusingToBryan(mangles, motor_angles);

                            // Justify whether the motor angles are within request
                            processing_info += String.Format("HKL = ({0}, {1}, {2}) => (2theta, omega, chi, phi) = {3}, {4}, {5}, {6}",
                                shiftedhkl.H, shiftedhkl.K, shiftedhkl.L,
                                motor_angles[0], motor_angles[1], motor_angles[2], motor_angles[3]);

                            bool inmotorrange = AreMotorAnglesWithinLimits(motor_angles, minmotorpos, maxmotorpos);
                            if (inmotorrange)
                            {
                                list_genHKL.Add(shiftedhkl);

                                processing_info += "\n";
                                finalist += String.Format("{0}\t{1}\t{2}\t\t{3} \t{4} \t{5} \t{6}\n",
                                    shiftedhkl.H, shiftedhkl.K, shiftedhkl.L,
                                    mangles.twotheta, mangles.omega, mangles.chi, mangles.phi);
                            }
                            else
                            {
                                processing_info += "... ... Discarded due to out of boundary.\n";
                            }

                        } // L
                    } // K
                } // H

                // Convert List to array
                m_outputHKLArray = list_genHKL.ToArray();

                return true;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GenerateReflectionList(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// GenerateIntensityReflectionList()
        /// </summary>
        /// <param name="hkl_intensity"></param>
        /// <param name="minmotorpos"></param>
        /// <param name="maxmotorpos"></param>
        /// <param name="condition"></param>
        /// <param name="finalist"></param>
        /// <param name="processing_info"></param>
        /// <returns></returns>
        public List<double[]> GenerateIntensityReflectionList(double[,] hkl_intensity, double[] minmotorpos, double[] maxmotorpos, string condition,
            out string finalist, out string processing_info)
        {
            try
            {
                finalist = "";
                processing_info = "";
                List<double[]> hkl_intensity_list = new List<double[]>();

                // Check array validity 
                if (hkl_intensity.GetLength(0) < 1 || hkl_intensity.GetLength(1) != 4)
                {
                    processing_info = "Error: Unable to process due to invalid size of hkl_instensity 2D array.";
                    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateIntensityReflectionList(). " + processing_info, "Error");
                    return null;
                }

                if (minmotorpos.Length != 4 || maxmotorpos.Length != 4)
                {
                    processing_info = "Error: Unable to process due to invalid number of motor positions.";
                    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateIntensityReflectionList(). " + processing_info, "Error");
                    return null;
                }

                if (_DBOUTPUT)
                {
                    Console.WriteLine("  Motor position is limited as 2theta = ({0}, {1}), Omega = ({2}, {3}), Phi = ({4}, {5}), Chi = ({6}, {7}).",
                        minmotorpos[0], maxmotorpos[0], minmotorpos[1], maxmotorpos[1],
                        minmotorpos[2], maxmotorpos[2], minmotorpos[3], maxmotorpos[3]);
                }

                // Process condition
                bool checkcondition;
                condition = condition.Replace(" ", "");
                condition = condition.ToLower();
                ExpressioneValuator conditionevaluator = null;
                if (condition.Length > 0)
                {
                    // Check condition
                    conditionevaluator = new ExpressioneValuator();
                    string[] parameternames = new string[3] { "h", "k", "l" };
                    conditionevaluator.SetupIntExpression(condition, parameternames);
                    string errormessage;
                    bool expvalid = conditionevaluator.IsReturnBoolean(out errormessage);
                    if (!expvalid)
                    {
                        processing_info = String.Format("Error: Unable to evaluate boolean expression {0} caused by {1}", condition, errormessage);
                        Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.GenerateIntensityReflectionList(). " + processing_info, "Error");
                        Console.WriteLine(processing_info);
                        return null;
                    }
                    checkcondition = true;
                }
                else
                {
                    // Empty condition. Allow.
                    checkcondition = false;
                }


                MillerIndices hkl = new MillerIndices();
                int h, k, l;

                for (int i = 0; i < hkl_intensity.GetLength(0); ++i)
                {
                    hkl.H = (double)hkl_intensity[i,0];
                    hkl.K = (double)hkl_intensity[i, 1];
                    hkl.L = (double)hkl_intensity[i, 2];

                    h = Convert.ToInt32(hkl.H);
                    k = Convert.ToInt32(hkl.K);
                    l = Convert.ToInt32(hkl.L);

                    // Filter by condition
                    if (h + k + l <= 0)
                        continue;
                    if (checkcondition)
                    {
                        bool meetcondition = CheckMillerIndexWithCondition(conditionevaluator, hkl);
                        if (!meetcondition)
                        {
                            // HKL does not meet condition
                            // Console.WriteLine("HKL = ({0}, {1}, {2}) does not meet reflection condition {3}",
                            //    hkl.H, hkl.K, hkl.L, condition);
                            processing_info += String.Format("HKL = ({0}, {1}, {2}) does not meet reflection condition {3}\n",
                                hkl.H, hkl.K, hkl.L, condition);
                            continue;
                        }
                    }

                    // Calculate motor angles
                    int ierror;
                    MotorIncidentAngles mangles = new MotorIncidentAngles();
                    m_ubmatrix.CalMotorAnglesFromMillerIndex(hkl, mangles, false, out ierror);

                    double[] motor_angles = new double[4];
                    ConvertMotorAnglesFromBusingToBryan(mangles, motor_angles);

                    // Justify whether the motor angles are within request
                    processing_info += String.Format("HKL = ({0}, {1}, {2}) => (2theta, omega, chi, phi) = {3}, {4}, {5}, {6}",
                        hkl.H, hkl.K, hkl.L,
                       motor_angles[0], motor_angles[1], motor_angles[2], motor_angles[3]);

                    bool inmotorrange = AreMotorAnglesWithinLimits(motor_angles, minmotorpos, maxmotorpos);
                    if (inmotorrange)
                    {
                        hkl_intensity_list.Add(new double[4] { hkl_intensity[i, 0], hkl_intensity[i, 1], hkl_intensity[i, 2], hkl_intensity[i, 3] });

                        processing_info += "\n";
                        finalist += String.Format("{0}\t{1}\t{2}\t\t{3} \t{4} \t{5} \t{6}\n",
                            hkl.H, hkl.K, hkl.L,
                            mangles.twotheta, mangles.omega, mangles.chi, mangles.phi);
                    }
                    else
                    {
                        processing_info += "... ... Discarded due to out of boundary.\n";
                    }
                }

                return hkl_intensity_list;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GenerateIntensityReflectionList(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }


        /// <summary>
        /// GenerateSPICEMacroCommands() generates SPICE macro commands based on hkl, motor limit, and condition settings.
        /// </summary>
        /// <param name="minhklindex"></param>
        /// <param name="maxhklindex"></param>
        /// <param name="minmotorpos"></param>
        /// <param name="maxmotorpos"></param>
        /// <param name="hklshift"></param>
        /// <param name="condition"></param>
        /// <param name="suffix_commands"></param>
        /// <param name="macro_commands"></param>
        /// <param name="processing_results"></param>
        /// <param name="number_of_generated_macro_points"></param>
        /// <returns></returns>
        public bool GenerateSPICEMacroCommands(int[] minhklindex, int[] maxhklindex, double[] minmotorpos, double[] maxmotorpos, double[] hklshift, string condition,
            string suffix_commands, out string macro_commands, out string processing_results, out int number_of_generated_macro_points)
        {
            try
            {
                string error_message;
                number_of_generated_macro_points = 0;

                if (CheckConditionValidity(condition, out error_message))
                {
                    string final_list;

                    bool generation_successful = GenerateReflectionList(minhklindex, maxhklindex, minmotorpos, maxmotorpos, hklshift, condition, out final_list, out processing_results);

                    if (generation_successful)
                    {
                        macro_commands = "";

                        // Check output HKL list length
                        if (m_outputHKLArray.Length == 0)
                        {
                            error_message = "Unable to generate commands because there is no generated HKL list.";
                            m_errorMessage += error_message;
                            Console.WriteLine(error_message);
                            processing_results += "\n" + error_message;
                            return false;
                        }

                        // Generate
                        for (int i = 0; i < m_outputHKLArray.Length; ++i)
                        {
                            string tempcommand = String.Format("br {0} {1} {2}\n", m_outputHKLArray[i].H,
                                    m_outputHKLArray[i].K, m_outputHKLArray[i].L);

                            if ((suffix_commands != null) && (suffix_commands.Length > 0))
                            {
                                tempcommand += suffix_commands + "\n";
                            }

                            ++number_of_generated_macro_points;
                            macro_commands += tempcommand;
                        }

                        return true;
                    }
                    else
                    {
                        // Failed to generate reflections.
                        macro_commands = "";
                        processing_results = "Unable to generate reflection list: " + processing_results;
                        return false;
                    }
                }
                else
                {
                    // Supplied condition is not valid.
                    macro_commands = "";
                    processing_results = "Condition invalid: " + error_message;
                    return false;
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GenerateSPICEMacroCommands(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        /// <summary>
        /// GenerateSPICEIntensityMacroCommands()
        /// </summary>
        /// <param name="hkl_intensity"></param>
        /// <param name="minmotorpos"></param>
        /// <param name="maxmotorpos"></param>
        /// <param name="condition"></param>
        /// <param name="suffix_commands"></param>
        /// <param name="scan_time_for_intensity_1_seconds"></param>
        /// <param name="max_scan_time_seconds"></param>
        /// <param name="macro_commands"></param>
        /// <param name="processing_results"></param>
        /// <param name="total_counting_time"></param>
        /// <returns></returns>
        public bool GenerateSPICEIntensityMacroCommands(double[,] hkl_intensity, double[] minmotorpos, double[] maxmotorpos, string condition,
            string suffix_commands, double scan_time_for_intensity_1_seconds, double max_scan_time_seconds, 
            out string macro_commands, out string processing_results,
            out double total_counting_time, out int number_of_generated_macro_points, out int max_scan_time_points)
        {
            try
            {
                string error_message;
                total_counting_time = 0.0;
                number_of_generated_macro_points = 0;
                max_scan_time_points = 0;

                if (CheckConditionValidity(condition, out error_message))
                {
                    string final_list;

                    List<double[]> filtered_list = GenerateIntensityReflectionList(hkl_intensity, minmotorpos, maxmotorpos, condition, out final_list, out processing_results);

                    if (filtered_list != null)
                    {
                        macro_commands = "";

                        // Check output HKL list length
                        if (filtered_list.Count == 0)
                        {
                            error_message = "Unable to generate commands because there is no generated HKL list.";
                            m_errorMessage += error_message;
                            Console.WriteLine(error_message);
                            processing_results += "\n" + error_message;
                            return false;
                        }

                        // Generate
                        for (int i = 0; i < filtered_list.Count; ++i)
                        {
                            string tempcommand = String.Format("br {0} {1} {2}\n", filtered_list[i][0],
                                    filtered_list[i][1], filtered_list[i][2]);

                            double intensity_to_use = Math.Abs(filtered_list[i][3]);

                            // Protect against intensity being near zero.
                            double scan_time_seconds = 0.0;
                            if (scan_time_for_intensity_1_seconds <= intensity_to_use * max_scan_time_seconds)
                            {
                                // Default calcualted scan time will be sufficiently short.
                                scan_time_seconds = scan_time_for_intensity_1_seconds / intensity_to_use;
                            }
                            else
                            {
                                // Default calculated scan time would be longer than the specified maximum.
                                scan_time_seconds = max_scan_time_seconds;
                                ++ max_scan_time_points;
                            }

                            tempcommand += String.Format("preset time {0}\n", scan_time_seconds);
                            ++number_of_generated_macro_points;
                            total_counting_time += scan_time_seconds;

                            if ((suffix_commands != null) && (suffix_commands.Length > 0))
                            {
                                tempcommand += suffix_commands + "\n";
                            }

                            macro_commands += tempcommand;
                        }

                        return true;
                    }
                    else
                    {
                        // Failed to generate reflections.
                        macro_commands = "";
                        processing_results = "Unable to generate reflection list: " + processing_results;
                        return false;
                    }
                }
                else
                {
                    // Supplied condition is not valid.
                    macro_commands = "";
                    processing_results = "Condition invalid: " + error_message;
                    return false;
                }
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.GenerateSPICEIntensityMacroCommands(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }


        /// <summary>
        /// Check whether input Miller index meets condition
        /// </summary>
        /// <param name="hkl"></param>
        /// <param name="condition"></param>
        /// <returns></returns>
        public bool CheckMillerIndexWithCondition(ExpressioneValuator evaluator, MillerIndices hkl)
        {
            try
            {
                // Construct variables
                string[] variablenames = new string[] { "h", "k", "l" };
                int[] variablevalues = new int[3];
                variablevalues[0] = (int)hkl.H;
                variablevalues[1] = (int)hkl.K;
                variablevalues[2] = (int)hkl.L;

                string valuestr;
                bool valid = evaluator.Calculate(variablevalues, out valuestr);
                if (!valid)
                {
                    Log.Write(Log.Error, "In DAS.UBMatrix.LabVIEWAdaptor.CheckMillerIndexWithCondition(). " + "Condition is not valid! \n" + valuestr, "Error");
                    return false;
                }

                if (valuestr == "true")
                {
                    return true;
                }

                return false;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.CheckMillerIndexWithCondition(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

      
        /// <summary>
        /// Filter motor angle/positions by predefined boundaries
        /// </summary>
        /// <param name="motor_angles"></param>
        /// <param name="minmotorpos">Lower boundary: 2theta, omega, phi and chi are not correlated. </param>
        /// <param name="maxmotorpos">Upper boundary: 2theta, omega, phi and chi are not correlated. (storage only)</param>
        /// <returns></returns>
        private bool AreMotorAnglesWithinLimits(double[] motor_angles, double[] minmotorpos, double[] maxmotorpos)
        {
            try
            {
                bool inrange = true;

                if (motor_angles[0] < minmotorpos[0] || motor_angles[0] > maxmotorpos[0])
                {
                    inrange = false;
                }
                else if (motor_angles[1] < minmotorpos[1] || motor_angles[1] > maxmotorpos[1])
                {
                    inrange = false;
                }
                else if (motor_angles[2] < minmotorpos[2] || motor_angles[2] > maxmotorpos[2])
                {
                    inrange = false;
                }
                else if (motor_angles[3] < minmotorpos[3] || motor_angles[3] > maxmotorpos[3])
                {
                    inrange = false;
                }

                return inrange;
            }
            catch (Exception exception)
            {
                string log_message = "Exception in DAS.UBMatrix.LabVIEWAdaptor.AreMotorAnglesWithinLimits(). " + exception.Message;
                Log.Write(Log.Info, log_message, "Error");
                m_errorMessage += "\n" + log_message;
                throw new Exception(log_message);
            }
        }

        #endregion

    }
}
