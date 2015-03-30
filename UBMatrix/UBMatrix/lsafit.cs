using System;
using System.Collections.Generic;
using System.Text;

///<note>The naming convention is that (1) use v for Vector and (2) m for Matrix</note>

namespace DAS.UBMatrix
{
    /// <summary>
    /// This is a temporary class.  
    /// Converted from UNITCELL.FOR
    /// </summary>
    public class GenUBLSF
    {


        #region constructor
        /// <summary>
        /// Constructor
        /// </summary>
        public GenUBLSF(double wavelength){
               // set the default value

            return;
        }
        #endregion

        #region public properties

        #endregion

        #region private methods

        #endregion 

        #region public methods












        /// <summary>
        /// Calculate 2theta from Miller indicies
        /// </summary>
        /// <param name="hkl"></param>
        /// <param name="mGI"></param>
        /// <param name="wavelength"></param>
        /// <returns></returns>
        public double Calculate2Theta(MillerIndices hkl, Matrix mGI, double wavelength, out double d1){
            d1 = Vector.DotProduct(hkl, mGI.MultVector(hkl));
		    double td = Math.Sqrt(wavelength*wavelength*d1/4.0);      				      			
            double tt = 2.0*DasMath.AsinD(td);  
      
            return tt;
        }





 


#if false
        /// <summary>
        /// This routine uses the unit-cell parameters, a list of observed 
	    /// diffractometer angles and an angular tolerance to find an orientation
	    /// matrix that can index the entire list.  Indices for two reflections
	    /// are chosen arbitrarily.  Original routine coded March 1995 by R. T.
	    /// Downs, Geophysical Laboratory.  The routine was converted to
	    /// use with 'SINGLE' by L.W. Finger, March 1995.
        /// </summary>
        static public void IndxC(UnitCell cell){

      		// real*4 GI(3,3),TOL(2)
      		// integer*4 ipoint
      		// common /useref/ipoint(listlen)
		
		    // 1. Get reciprocal cell metric
		    if (iutility == 0) {

			    // 1. Get list of Reflections to use
      			ncount = 0;
			    for (int i = 0; i < nlst; i ++){
				    if (iuse[i] == 0){
					    ncount ++;
				    }
                }

			    if (ncount <= 3){
				    Consonel.WriteLine("Number of Reflections is too small for INDC");
				    throw Excpetion();
				    return;
                }

			    Metr(out cell, out gi, out vol2);
				
			    iutility = 31;

			    Console.WriteLine("Tolerance for 2Theta and Eulerian angles");
            
            } else {
			    // Read in tolerance
			    double[] Tolerance  = new double[2];
			    for (int i = 0; i < 2; i ++){
				    Tolerance[i] = 1.0;
			    }
			    ReadTolerance();

			    // Generate CALCULATED AND OBSERVED ANGLES TO GET UB-MATRIX			
			    GenUB(gi, ncount, tol);
            }

		    return;
        }
#endif

        #endregion


    }

}
