using System;
using System.Collections.Generic;
using System.Text;
using DAS.UBMatrix;
using Microsoft.CSharp;
using System.CodeDom.Compiler;
using System.Reflection;

namespace Testers
{
    class Tester
    {
        static public void MainTest(string[] args)
        {
            string condition = "h == 2 && l == 3 || h - k + l == 4 || (h+k)%2 == 0";
            // condition = "h + 2";
            Console.WriteLine("Test input: \"{0}\"", condition);

            string[] varnames = new string[3] { "h", "k", "l" };
            int[] varvalues = new int[3] { 1, 1, 3 };
            string valuestr;
            bool valid = ExpressioneValuator.Eval(condition, varnames, varvalues, out valuestr);
            if (valid)
            {
                Console.WriteLine("Expression is valid.  Result = {0}", valuestr);
            }
            else
            {
                Console.WriteLine("Expression is invalid.  Error Message:{0}", valuestr);
            }

#if false
            LabVIEWAdaptor adp = new LabVIEWAdaptor(50);
            adp.ParseCondition(condition);

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

            StringBuilder sb = new StringBuilder("");
            sb.Append("using System;\n");
            sb.Append("using System.Xml;\n");
            sb.Append("using System.Data;\n");
            sb.Append("using System.Data.SqlClient;\n");
            sb.Append("using System.Windows.Forms;\n");
            sb.Append("using System.Drawing;\n");

            sb.Append("namespace CSCodeEvaler{ \n");
            sb.Append("public class CSCodeEvaler{ \n");
            sb.Append("public object EvalCode(){\n");
            sb.Append("return " + "1 + 5%2 == 3" + "; \n");
            sb.Append("} \n");
            sb.Append("} \n");
            sb.Append("}\n");

            Console.WriteLine("Code: \n{0}", sb.ToString());

            CompilerResults cr = icc.CompileAssemblyFromSource(cp, sb.ToString());
            if (cr.Errors.Count > 0)
            {
                return;
            }

            System.Reflection.Assembly a = cr.CompiledAssembly;
            object o = a.CreateInstance("CSCodeEvaler.CSCodeEvaler");

            Type t = o.GetType();
            MethodInfo mi = t.GetMethod("EvalCode");

            object s = mi.Invoke(o, null);

            // This is the answer
            Console.WriteLine("{0}", s.ToString());
#endif 
            
        }

        /// <summary>
        /// Test Plan:
        /// 1. Vector
        /// 2. Matrix
        /// 3. UBMatrix
        /// </summary>
        /// <param name="args"></param>
        static void Main01(string[] args)
        {
            
            TestUBMatrix();

            // -1: waiting for user finish
            Console.ReadLine();
        }

        /// <summary>
        /// Test UBMatrix Generator
        /// </summary>
        static void TestUBMatrix()
        {
            // 1. set up 2 reflection
            List<Reflection> reflist = Reflection.ImportFromFileDAS("ubmatrix.txt");
            for (int i = 0; i < reflist.Count; i++)
            {
                reflist[i].Print();
            }

#if false
            List<MotorIncidentAngles> anglelist = new List<MotorIncidentAngles>();
            for (int i = 0; i < reflist.Count; i++)
            {
                anglelist.Add(reflist[i].motor_angles);
            }
#endif
            MotorIncidentAngles[] anglelist = new MotorIncidentAngles[reflist.Count];
            for (int i = 0; i < reflist.Count; i++)
                anglelist[i] = reflist[i].MotorAngles;

            // 2. set unit cell
            UnitCell sicell = new UnitCell();
#if false
            sicell.a = sicell.b = sicell.c = 5.430940;
            sicell.alpha = sicell.beta = sicell.gamma = 90;
            UBMatrix ubmatrix = new UBMatrix(1.01);
#endif
            // FeWO4
            sicell.a = 4.753;
            sicell.b = 5.72;
            sicell.c = 4.968;
            sicell.alpha = sicell.gamma = 90.0;
            sicell.beta = 90.08;
            UBMatrix ubmatrix = new UBMatrix(1.56);
            

            // sicell.CalcualteVolume();
            // throw new Exception("Debug Stop");
            
            
            ubmatrix.SetDirectUnitCell(sicell);
          
            ubmatrix.Tolerance_2Theta = 0.5;
            ubmatrix.Tolerance_Euler = 5.0;
#if false
            ubmatrix.CalBMatrix();
            ubmatrix.B.Print();           
            //ubmatrix.CalUBMatrix(r1, r2);
#endif
            ubmatrix.GenerateUBFromUnIndexedReflectiuons(anglelist, -1);

            return;

        }

        /// <summary>
        /// Test for matrix
        /// </summary>
        static void TestMatrix()
        {
#if false
           Matrix tmatrix = new Matrix(2);
            Vector v1 = new Vector(2);
            Vector v2 = new Vector(2);
            v1.Set(0, 1);
            v1.Set(1, 2);
            v2.Set(0, 4);
            v2.Set(1, -3);
            tmatrix.SetColumn(0, v1);
            tmatrix.SetColumn(1, v2);

            tmatrix.Print();
            Matrix tm2 = new Matrix(2);
            tmatrix.Transpose(tm2);
            tm2.Print();
#endif
            Matrix tmatrix = new Matrix(3);
            Vector v1 = new Vector(3);
            Vector v2 = new Vector(3);
            Vector v3 = new Vector(3);

            Vector[] uvec = new Vector[3];
            for (int i = 0; i < 3; i++)
            {
                uvec[i] = new Vector(3);
                for (int j = 0; j < 3; j++)
                {
                    if (i == j)
                    {
                        uvec[i].Set(j, i + 1);
                    }
                    else
                    {
                        uvec[i].Set(j, 0);
                    }
                }
            }


            v1.Set(0, 1);
            v1.Set(1, 2);
            v1.Set(2, 4);
            v2.Set(0, 4);
            v2.Set(1, -3);
            v2.Set(2, -2);
            v3.Set(0, 10);
            v3.Set(1, -11);
            v3.Set(1, 15);
            tmatrix.SetColumn(0, v1);
            tmatrix.SetColumn(1, v2);
            tmatrix.SetColumn(2, v3);

            tmatrix.Print();

            Matrix tm2 = new Matrix(3);
            tmatrix.Transpose(tm2);
            tm2.Print();

            for (int i = 0; i < 3; i++)
            {
                tm2.SetColumn(i, uvec[i]);
            }
            tm2.Print();

            Matrix tm3 = new Matrix(3);
            Matrix.MulMatrix(tmatrix, tm2, tm3);
            tm3.Print();

            for (int i = 0; i < 3; i++)
            {
                Vector p1 = tm3.Project(i);
                p1.Print();
            }

            Matrix tm4 = tmatrix.Duplicate();
            tm4.Print();

            int[] indexlist = new int[3];
            double d;
            Matrix.LUDcomp(tm4, indexlist, out d);
            tm4.Print();


            Matrix inmatrix = new Matrix(3);
            double det;
            tmatrix.InvertMatrix(inmatrix, out det);
            Console.WriteLine("Inversed Matrix");
            inmatrix.Print();

            Matrix unitmatrix = new Matrix(3);
            Matrix.MulMatrix(tmatrix, inmatrix, unitmatrix);
            Console.WriteLine("Check");
            unitmatrix.Print();

            return;
        }

        /// <summary>
        /// Generate a random matrix, do the inversion
        /// and test
        /// </summary>
        static void TestMatrix2()
        {
            int matrixsize = 6;
            Random rand = new Random(0);

            for (int i = 0; i < 10; i++)
            {
                Matrix tM = new Matrix(matrixsize);
                // randomly assign value to matrix
                for (int j = 0; j < matrixsize; j++)
                {
                    double[] vec = new double[matrixsize];
                    for (int k = 0; k < matrixsize; k++)
                    {
                        vec[k] = rand.Next(0, 10000) / 100.00 - 50;
                    }
                    tM.SetColumn(matrixsize,j, vec);
                }
                Console.WriteLine("Test Matrix:");
                tM.Print();
                // invert
                Matrix iM = new Matrix(matrixsize);
                double det;
                tM.InvertMatrix(iM, out det);
                Console.WriteLine("Inversed Test Matrix:");
                iM.Print();
                // Check
                Matrix uM = new Matrix(matrixsize);
                Matrix.MulMatrix(iM, tM, uM);
                Console.WriteLine("M^-1 * M:");
                uM.Print();
                Matrix.MulMatrix(tM, iM, uM);
                Console.WriteLine("M * M^-1:");
                uM.Print();

                Console.WriteLine("Test {0} Over", i);
                Console.ReadLine();
            }
        }

        /// <summary>
        /// Test for vectors
        /// </summary>
        static public void TestVector()
        {
            // 1. Part 1: Vector operation
            Vector a, b, c;

            a = new Vector(3);
            a.i = 1;
            a.j = 2;
            a.k = 3;
            a.Print();
            a.Set(2, -1);
            a.Print();
            Console.WriteLine("{0}, {1}, {2}", a.Get(0), a.Get(1), a.Get(2));
            b = a.UnitVector();
            b.Print();
            Console.WriteLine("Unit Vector b's Length = {0}", b.Length());

            double x = Vector.DotProduct(a, b);
            Console.WriteLine("a dot b = {0}", x);

            Vector d11 = new Vector(3);
            d11.Set(0, 1);
            d11.Set(1, 1);
            d11 = d11.UnitVector();
            Console.WriteLine("Unit Vector @ (110)");
            d11.Print();
            double l11 = a.Project(0, 1);
            Console.WriteLine("Vector a 's length projected at (110) = {0}", l11);

            // 2. Cross production
            b.i = 2;
            b.j = -1;
            b.k = 5;

            Console.WriteLine("Test Cross Production");
            c = new Vector(3);
            Vector.CrossProduct(a, b, c);
            a.Print();
            b.Print();
            c.Print();


        }
    }
}
