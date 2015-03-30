using System;
using System.Collections.Generic;
using System.Text;

using Microsoft.CSharp;
using System.CodeDom.Compiler;
using System.Reflection;

namespace DAS.UBMatrix
{
    public class ExpressioneValuator
    {
        #region private class variables
        /// <summary>
        /// Condition string
        /// </summary>
        string m_condition;
        /// <summary>
        /// Flag whether all arguments are integer
        /// </summary>
        //bool m_isIntExpression;
        /// <summary>
        /// Error message
        /// </summary>
        string m_errormessage;
        /// <summary>
        /// Whether is a valid expression
        /// </summary>
        bool m_isValidExpression;
        /// <summary>
        /// Number of parameters
        /// </summary>
        int m_numParameters;
        /// <summary>
        /// Assembly instance 
        /// </summary>
        object m_assemblyInstance;
        /// <summary>
        /// MethodInfo to evalulate
        /// </summary>
        MethodInfo m_methodInfo;
        #endregion

        #region constructor
        /// <summary>
        /// Constructor
        /// </summary>
        public ExpressioneValuator()
        {
            m_condition = "";
            m_errormessage = "";
            m_isValidExpression = false;
            m_numParameters = 0;
            m_methodInfo = null;
            m_assemblyInstance = null;
            //m_isIntExpression = false;
        }
        #endregion

        #region public non-static methods
        public void SetupIntExpression(string condition, string[] parameternames)
        {
            // Check condition: do striping
            string formula = condition.Replace(" ", "").Replace("\t", "").Replace("\n", "");
            if (formula.Length == 0)
            {
                m_isValidExpression = false;
                return;
            }
            else
            {
                m_condition = formula;
            }

            // Generate source code
            string sourcecode = GenerateCSharpIntEvalCode(m_condition, parameternames);
            m_numParameters = parameternames.Length;

            // Compile 
            CSharpCodeProvider c = new CSharpCodeProvider();
            ICodeCompiler icc = c.CreateCompiler();
            CompilerParameters cp = new CompilerParameters();

            cp.ReferencedAssemblies.Add("system.dll");
            cp.ReferencedAssemblies.Add("system.xml.dll");
            cp.ReferencedAssemblies.Add("system.data.dll");
            cp.ReferencedAssemblies.Add("system.windows.forms.dll");
            cp.ReferencedAssemblies.Add("system.drawing.dll");

            cp.CompilerOptions = "/t:library";
            cp.GenerateInMemory = true;

            CompilerResults cr = icc.CompileAssemblyFromSource(cp, sourcecode);

            if (cr.Errors.Count > 0)
            {
                // If compiling error, return false to indicate invalid formula and source code for reference  
                m_isValidExpression = false;
                m_errormessage = "Condition caused compilation error: \n";
                for (int i = 0; i < cr.Output.Count; ++i)
                {
                    m_errormessage += cr.Output[i] + "\n";
                }
                m_errormessage += "Source codes:\n" + sourcecode;
                return;
            }
            m_isValidExpression = true;
            //m_isIntExpression = true;

            System.Reflection.Assembly a = cr.CompiledAssembly;
            m_assemblyInstance = a.CreateInstance("CSCodeEvaler.CSCodeEvaler");
            Type t = m_assemblyInstance.GetType();
            m_methodInfo = t.GetMethod("EvalCode");

            return;                
        }

        /// <summary>
        /// Calculate the condition-expression with input
        /// </summary>
        /// <param name="parametervalues"></param>
        /// <param name="?"></param>
        /// <returns></returns>
        public bool Calculate(int[] parametervalues, out string valuestr)
        {
            // Check
            if (!m_isValidExpression)
            {
                valuestr = String.Format("Condition a valid expression:  {0}", m_errormessage);
                return false;
            }
            else if (parametervalues.Length != m_numParameters)
            {
                m_errormessage = String.Format("Condition expression has {0} input parameters, while input number of parameter values is {1}",
                    m_numParameters, parametervalues.Length);
                valuestr = m_errormessage;
                return false;
            }

            // Construct parameter object array            
            object[] parameters = new object[parametervalues.Length];
            for (int i = 0; i < m_numParameters; ++i)
                parameters[i] = parametervalues[i];

            // Invoke
            object s = m_methodInfo.Invoke(m_assemblyInstance, parameters);

            // This is the answer
            // Console.WriteLine("{0}", s.ToString());
            valuestr = s.ToString().ToLower();

            return true;
        }

        /// <summary>
        /// Justify whether return is boolean
        /// </summary>
        /// <param name="errormessage"></param>
        /// <returns></returns>
        public bool IsReturnBoolean(out string errormessage)
        {
            int[] testvalues = new int[m_numParameters];
            Random rnd = new Random();
            for (int i = 0; i < m_numParameters; ++i)
            {
                testvalues[i] = rnd.Next(1, 100);
            }

            string valuestr;
            bool valid = Calculate(testvalues, out valuestr);
            if (!valid)
            {
                errormessage = valuestr;
                return false;
            }
            else if (valuestr != "true" && valuestr != "false")
            {
                errormessage = String.Format("Condition's returned value is not boolean but {0}.", valuestr);
                return false;
            }

            errormessage = "";
            return true;
        }

        #endregion

        #region static static public

        /// <summary>
        /// Evalulate a formula/expression
        /// </summary>
        /// <param name="formula"></param>
        /// <param name="variablenames"></param>
        /// <param name="variablevalues"></param>
        /// <param name="valuestr">value of formula if valid; error message if invalid</param>
        /// <returns></returns>
        static public bool Eval(string formula, string[] variablenames, int[] variablevalues, out string valuestr)
        {
            // Check inputs
            if (variablenames.Length != variablevalues.Length)
            {
                //  mLog.Record("Input number of variable names is not equal to variable values.  Unable to evaluate the function. ");
                valuestr = "";
                return false;
            }

            // Create a string of code
            string sourcecode = GenerateCSharpIntEvalCode(formula, variablenames);

            // Compile the source code
            CSharpCodeProvider c = new CSharpCodeProvider();
            ICodeCompiler icc = c.CreateCompiler();
            CompilerParameters cp = new CompilerParameters();

            cp.ReferencedAssemblies.Add("system.dll");
            cp.ReferencedAssemblies.Add("system.xml.dll");
            cp.ReferencedAssemblies.Add("system.data.dll");
            cp.ReferencedAssemblies.Add("system.windows.forms.dll");
            cp.ReferencedAssemblies.Add("system.drawing.dll");

            cp.CompilerOptions = "/t:library";
            cp.GenerateInMemory = true;

            CompilerResults cr = icc.CompileAssemblyFromSource(cp, sourcecode);

            if (cr.Errors.Count > 0)
            {
                // If compiling error, return false to indicate invalid formula and source code for reference   
                valuestr = "";
                for (int i = 0; i < cr.Output.Count; ++i)
                {
                    valuestr += cr.Output[i] + "\n";
                }
                valuestr += "Source codes:\n" + sourcecode;
                return false;
            }

            System.Reflection.Assembly a = cr.CompiledAssembly;
            object o = a.CreateInstance("CSCodeEvaler.CSCodeEvaler");
            Type t = o.GetType();
            MethodInfo mi = t.GetMethod("EvalCode");
            // object s = mi.Invoke(o, null);
            object[] parameters = new object[variablevalues.Length];
            for (int i = 0; i < variablevalues.Length; ++i)
                parameters[i] = variablevalues[i];

            object s = mi.Invoke(o, parameters);

            // This is the answer
            // Console.WriteLine("{0}", s.ToString());
            valuestr = s.ToString().ToLower();

            return true;
        }

        #endregion  

        #region

        /// <summary>
        /// Generate C# program to evaluate HKL's condition
        /// </summary>
        /// <param name="formula"></param>
        /// <param name="variablenames"></param>
        /// <param name="variablevalues"></param>
        /// <returns></returns>
        static private string GenerateCSharpIntEvalCode(string formula, string[] variablenames)
        {
            StringBuilder sb = new StringBuilder("");

            // Library to use
            sb.Append("using System;\n");
            sb.Append("using System.Xml;\n");
            sb.Append("using System.Data;\n");
            sb.Append("using System.Data.SqlClient;\n");
            sb.Append("using System.Windows.Forms;\n");
            sb.Append("using System.Drawing;\n");

            // Name space
            sb.Append("namespace CSCodeEvaler{ \n");

            // Class
            sb.Append("public class CSCodeEvaler{ \n");
            // Function name
            sb.Append("  public object EvalCode(");
            // Define parameters
            int numvariables = variablenames.Length;
            for (int i = 0; i < numvariables; ++i)
            {
                sb.Append("int " + variablenames[i]);
                if (i != numvariables-1)
                    sb.Append(", ");
            }
            sb.Append("){\n");

            /* This is for general case.  But the speed is low for multiple calling
            for (int i = 0; i < numvariables; ++i)
            {
                sb.Append("  int " + variablenames[i] + " = " +  variablevalues[i].ToString() + ";\n");
            }
             */

            // Formular
            sb.Append("    return " + formula + "; \n");
            sb.Append("  } \n");
            sb.Append("} \n");
            sb.Append("}\n");

            // Console.WriteLine("Code: \n{0}", sb.ToString());

            return sb.ToString();
        }

        #endregion

    }
}
