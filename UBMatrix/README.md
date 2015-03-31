New features implemented in March/April 2015
 1. Expose the calculation of estimated error for *refined* UB matrix; 



Features:
 1. UB matrix calculation
  * Calculate UB matrix from 2 reflections;
  * Calculate UB matrix from 3 reflections;
  * Refine UB matrix by least square fitting
 2. Experiment plan from known UB matrix
  * Calculate motor positions from Miller index (HKL)
 3. Testing: 
  * Functional test (DebugTest.cs)
  * System test (SystemTest.cs)
  * Command line interactive application to test without GUI (Workflow.cs)riv102er
  * 
  
Error estimation of refined UB matrix
 1. MillerIndices::SetFromCalculatedValue(newvalue): $\sigma=\sqrt{(H-H_0)^2+(K-K_0)^2+(L-L_0)^2}
 2. UBMatrix::CalculateHKLFromIncidentAngle()
 3. UBMatrix::CalculateUBMatrixError(reflectionlist): $\sigma=\sqrt(\sum_{r}^{reflection}(HKL_m-HKL_{cal})^2/N)$, where N is the number of reflections.
 
