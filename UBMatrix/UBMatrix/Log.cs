using System;
using System.Text;
using System.Diagnostics;
using System.IO;


namespace DAS.UBMatrix
{
 
    /// <summary>Log singleton</summary>
    public sealed class Log
    {
        /// <summary>
        /// Static constructor to initiate logging.
        /// </summary>
        static Log()
        {
            string log_file_directory = "C:\\spicetmp\\UBMatrixLog";
            string log_file_name = log_file_directory + "\\UBMatrixLog";

            // Assure log directory exists.
            System.IO.Directory.CreateDirectory(log_file_directory);

            string log_level_string = "Verbose";
        
            Log.TraceOn(log_file_name, log_level_string);
            Log.Write(true, string.Format("In DAS.UBMatrix.Log() static constructor. Starting UBMatrix log. log_file_name: {0} log_level: {1}", log_file_name, log_level_string), "Info");
            Log.Write(Log.Switch.TraceInfo, "In DAS.UBMatrix.Log() static constructor.", "Off");
        }

        /// <summary>Trace switch</summary>
        public static readonly TraceSwitch Switch = new TraceSwitch("SNS.DAS", "SNS.DAS Trace Switch");

        /// <summary>Verbose log</summary>
        public static bool Verbose { get { return Switch.TraceVerbose; } }
        /// <summary>Info log</summary>
        public static bool Info { get { return Switch.TraceInfo; } }
        /// <summary>Warning log</summary>
        public static bool Warning { get { return Switch.TraceWarning; } }
        /// <summary>Error log</summary>
        public static bool Error { get { return Switch.TraceError; } }

        /// <summary>Log a message</summary>
        /// <param name="log">log flag (if true the message will be logged)</param>
        /// <param name="message">message to log</param>
        /// <param name="category">message category</param>
        public static void Write(bool log, string message, string category)
        {
            if (!log) return;
            Trace.WriteLine(String.Format("{0:yyyy-MM-dd HH:mm:ss} {1}", System.DateTime.Now, message), category);                            
        }

        /// <summary>Exception log</summary>
        /// <param name="ex">exception to log</param>
        public static void Exception(Exception ex)
        {
            if (ex != null)
            {
#if DEBUG
                Trace.WriteLine(String.Format("{0:yyyy-MM-dd HH:mm:ss} {1}", System.DateTime.Now, ex)        , "Exception");                
#else
                Trace.WriteLine(String.Format("{0:yyyy-MM-dd HH:mm:ss} {1}", System.DateTime.Now, ex.Message), "Exception");				
#endif
            }
        }

        /// <summary>Exception log</summary>
        /// <param name="message">exception message</param>
        public static void Exception(string message)
        {
            Trace.WriteLine(String.Format("{0:yyyy-MM-dd HH:mm:ss} {1}", System.DateTime.Now, message), "Exception");
        }

        /// <summary>
        /// Turn on logging
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="log_level_string"></param>
        public static void TraceOn(string filename, string log_level_string)
        {
            try
            {
                TraceLevel log_level;

                int log_level_int;

                // First, try to parse log level as a number.
                if (int.TryParse(log_level_string, out log_level_int))
                {
                    if (log_level_int < 0)
                    {
                        log_level_int = 0;
                    }
                    else if (log_level_int > 4)
                    {
                        log_level_int = 4;
                    }
                    log_level = (TraceLevel)log_level_int;
                }
                else
                {
                    // Log level not given as a number. Try the level name.

                    //rdg; Since we are coding to .NET Framework 3.5 (max level LabVIEW 2012 works with by default),
                    // and since .NET 3.5 does not have an Enu.TryParse() we punt here by catching an exception for
                    // a failed parse.
                    try
                    {
                        log_level = (TraceLevel)Enum.Parse(typeof(TraceLevel), log_level_string, true);
                    }
                    catch (Exception)
                    {
                        // Unable to parse.
                        log_level = TraceLevel.Off;
                    }
                }

 
                // Create a file for output named TestFile.txt.
                // Add timestamp.
                string full_filename = string.Format("{0} ({1:yyyy-MM-dd_hh-mm-ss-tt}).log", filename, DateTime.Now);
                Stream log_file = new FileStream(full_filename, FileMode.Create, FileAccess.Write, FileShare.Delete | FileShare.ReadWrite);

                TextWriterTraceListener defaultListener = new TextWriterTraceListener(log_file);

                Trace.Listeners.Add(defaultListener);
                Trace.AutoFlush = true;
                Log.Switch.Level = log_level;
            }
            catch (Exception ex)
            {
                throw new Exception("Exception in Log.TraceOn().", ex);
            }
        }

        /// <summary>Hide this</summary>
        private Log() { }
    }
}
