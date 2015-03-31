New features:
 1. Expose the calculation of estimated error for *refined* UB matrix; 
  * Equation: :math:`\sigma` 
 2. 



I. Test

Test Plan for UBMatrix library

1. Vector
2. Matrix simple operation (multiply, transpose, and etc)
3. Matrix advanced operation (LU decomposition, matrix inverse)
4. UB matrix
5. Unit test

II. Expansion from SINGLE-Fortran 1st Quarter
1. Can use UB matrix to determine the motor position of given (H, K, L)

III. Expansion from SINGLE-Fortran 2nd Quarter
1. Least square fitting works
2. Refinement of UB matrix
