using System;
using System.Collections.Generic;
using System.Text;

namespace DAS.UBMatrix
{
    /// <summary>
    /// some math methods used in DAS.UBMatrix Package
    /// </summary>
    public class DasMath
    {
        #region constants
        /// <summary>
        /// constant, rad = pi/180 (unit in /degree)
        /// </summary>
        public const double rad = Math.PI / 180.0;
        /// <summary>
        /// constant, inverse of rad (unit in degree)
        /// </summary>
        public const double radin = 180.0 / Math.PI;
        #endregion

        #region static public methods
        /// <summary>
        /// sin(degree)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        static public double SinD(double x)
        {
            double sind = Math.Sin(x * rad);
            return sind;
        }

        /// <summary>
        /// Cos(degree)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        static public double CosD(double x)
        {
            double cosd = Math.Cos(x * rad);
            return cosd;
        }

        /// <summary>
        /// cos^-1(degree)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        static public double AcosD(double x)
        {
            double acosd = radin * Math.Acos(x);
            // Console.WriteLine("Fun1109Math   x = {0}, cos^-1 = {1}", x, Math.Acos(x));
            return acosd;
        }
        /// <summary>
        /// sin^-1(degree)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        static public double AsinD(double x)
        {
            double asind = radin * Math.Asin(x);
            return asind;
        }

        /// <summary>
        /// Return the nearest integer of a double x
        /// </summary>
        /// <param name="x">double, to be converted</param>
        /// <returns></returns>
        static public int NearestInteger(double x)
        {
            int rt;
            int ni = Convert.ToInt32(x);
            double nx = Convert.ToDouble(ni);
            if (nx > x)
            {
                if ((nx - x) > (x - nx + 1))
                {
                    rt = ni - 1;
                }
                else
                {
                    rt = ni;
                }
            }
            else
            {
                if ((nx + 1 - x) > (x - nx))
                {
                    rt = ni;
                }
                else
                {
                    rt = ni + 1;
                }
            }

            return rt;
        }

        /// <summary>
        /// tan(x to y)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        static public double Atan2(double x, double y)
        {
            //double c = 0.0;
            throw new Exception("Implementing ASAP");
            //return c;
        }

        /// <summary>
        /// Sorting n arrays together, with one main arrays and unlimit arrays sorted
        /// according to the main array
        /// </summary>
        /// <param name="mainarray"></param>
        /// <param name="islavearrays"></param>
        /// <param name="dslavearrays"></param>
        static public void Sort(double[] mainarray, int[][] islavearrays, double[][] dslavearrays)
        {
            // 1. Check
            bool valid = true;
            for (int i = 0; i < islavearrays.Length; i++)
                if (mainarray.Length != islavearrays[i].Length)
                    valid = false;
            for (int i = 0; i < dslavearrays.Length; i++)
                if (mainarray.Length != dslavearrays[i].Length)
                    valid = false;
            if (valid == false)
                throw new Exception("Sort(): Input array length not equal");

#if _DBOUTPUT
            Console.WriteLine("Before Sorting:  DB 313");
            for (int i = 0; i < mainarray.Length; i++)
            {
                string buf = String.Format("{0}  ", mainarray[i]);
                for (int j = 0; j < islavearrays.Length; j++)
                    buf += String.Format("{0}  ", islavearrays[j][i]);
                for (int j = 0; j < dslavearrays.Length; j++)
                    buf += String.Format("{0}  ", dslavearrays[j][i]);
                Console.WriteLine(buf);
            }
#endif                           
            // 2. Sorting
            int counts = mainarray.Length;
            HeapSort(mainarray, islavearrays, dslavearrays, counts);

            return;
        }

        /// <summary>
        /// Swap a series of arrays
        /// </summary>
        /// <param name="a"></param>
        /// <param name="iaslave"></param>
        /// <param name="daslave"></param>
        /// <param name="s"></param>
        /// <param name="t"></param>
        static public void Swap(double[] a, int[][] iaslave, double[][] daslave, int s, int t)
        {
            Swap(a, s, t);
            for (int i = 0; i < iaslave.Length; i++)
                Swap(iaslave[i], s, t);
            for (int i = 0; i < daslave.Length; i++)
                Swap(daslave[i], s, t);
            return;
        }


        /// <summary>
        /// Swap two elements in the array
        /// </summary>
        /// <param name="a"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        static public void Swap(double[] a, int i, int j)
        {
            double temp = a[i];
            a[i] = a[j];
            a[j] = temp;
            return;
        }

        /// <summary>
        /// Swap two elements in the array
        /// </summary>
        /// <param name="a"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        static public void Swap(int[] a, int i, int j)
        {
            int temp = a[i];
            a[i] = a[j];
            a[j] = temp;
            return;
        }

        /// <summary>
        /// Heap sorting an array with length known
        /// </summary>
        /// <param name="a"></param>
        /// <param name="count"></param>
        static public void HeapSort(double[] a, int count)
        {
            // input:  an unordered array a of length count
            // first place a in max-heap order)

            Heapify(a, count);

            int end = count - 1;
            while (end > 0)
            {
                // (swap the root(maximum value) of the heap with the last element of the heap)
                Swap(a, end, 0);
                // (put the heap back in max-heap order)
                SiftDown(a, 0, end - 1);
                // (decrease the size of the heap by one so that the previous max value will
                // stay in its proper placement)
                end = end - 1;
            }

            return;
        }

        /// <summary>
        /// Main function of heap sorting
        /// </summary>
        /// <param name="a"></param>
        /// <param name="iasalve"></param>
        /// <param name="daslave"></param>
        /// <param name="count"></param>
        static public void HeapSort(double[] a, int[][] iasalve, double[][] daslave, int count)
        {
            // input:  an unordered array a of length count
            // first place a in max-heap order)

            Heapify(a, iasalve, daslave, count);

            int end = count - 1;
            while (end > 0)
            {
                // (swap the root(maximum value) of the heap with the last element of the heap)
                Swap(a, iasalve, daslave, end, 0);
                // (put the heap back in max-heap order)
                SiftDown(a, iasalve, daslave, 0, end - 1);
                // (decrease the size of the heap by one so that the previous max value will
                // stay in its proper placement)
                end = end - 1;
            }

            return;
        }

        /// <summary>
        /// Make it to heap
        /// </summary>
        /// <param name="a"></param>
        /// <param name="count"></param>
        static public void Heapify(double[] a, int count)
        {
            // (start is assigned the index in a of the last parent node)
            int start = (count - 2) / 2;

            while (start >= 0)
            {
                // (sift down the node at index start to the proper place such that all nodes below
                //  the start index are in heap order)
                SiftDown(a, start, count - 1);
                start = start - 1;
                // (after sifting down the root all nodes/elements are in heap order)
            }

            return;
        }

        /// <summary>
        /// Heapify
        /// </summary>
        /// <param name="a"></param>
        /// <param name="iaslave"></param>
        /// <param name="daslave"></param>
        /// <param name="count"></param>
        static public void Heapify(double[] a, int[][] iaslave, double[][] daslave, int count)
        {
            // (start is assigned the index in a of the last parent node)
            int start = (count - 2) / 2;

            while (start >= 0)
            {
                // (sift down the node at index start to the proper place such that all nodes below
                //  the start index are in heap order)
                SiftDown(a, iaslave, daslave, start, count - 1);
                start = start - 1;
                // (after sifting down the root all nodes/elements are in heap order)
            }

            return;
        }

        /// <summary>
        /// Shift down the sequence
        /// </summary>
        /// <param name="a"></param>
        /// <param name="start"></param>
        /// <param name="end"></param>
        static public void SiftDown(double[] a, int start, int end)
        {
            // input:  end represents the limit of how far down the heap
            //               to sift.

            int root = start;

            // (While the root has at least one child)
            while (root * 2 + 1 <= end)
            {
                // (root*2+1 points to the left child)
                int child = root * 2 + 1;
                // (If the child has a sibling and the child's value is less than its sibling's...)
                if (child + 1 <= end && a[child] < a[child + 1])
                {
                    child = child + 1;
                    // (... then point to the right child instead)
                }
                if (a[root] < a[child])
                {
                    // then (out of max-heap order)
                    Swap(a, root, child);
                    root = child;
                    // (repeat to continue sifting down the child now)
                }
                else
                {
                    return;
                }
            }
            return;
        }

        /// <summary>
        /// Sift down
        /// </summary>
        /// <param name="a"></param>
        /// <param name="iaslave"></param>
        /// <param name="daslave"></param>
        /// <param name="start"></param>
        /// <param name="end"></param>
        static public void SiftDown(double[] a, int[][] iaslave, double[][] daslave, int start, int end)
        {
            // input:  end represents the limit of how far down the heap
            //               to sift.

            int root = start;

            // (While the root has at least one child)
            while (root * 2 + 1 <= end)
            {
                // (root*2+1 points to the left child)
                int child = root * 2 + 1;
                // (If the child has a sibling and the child's value is less than its sibling's...)
                if (child + 1 <= end && a[child] < a[child + 1])
                {
                    child = child + 1;
                    // (... then point to the right child instead)
                }
                if (a[root] < a[child])
                {
                    // then (out of max-heap order)
                    Swap(a, iaslave, daslave, root, child);
                    root = child;
                    // (repeat to continue sifting down the child now)
                }
                else
                {
                    return;
                }
            }
            return;
        }
        #endregion
    }
}
