MODULE basic_constants
USE numeric_kinds

  REAL (kind=dp),PARAMETER :: pi=3.14159265358979323846264338_dp
  real (kind=dp),parameter :: e=2.718281828459045235360287471_dp

  REAL (kind=dp), PARAMETER :: golden_mean =1.618033988749894902525738871_dp 
  COMPLEX  (kind=dp), PARAMETER :: imag=CMPLX(0._dp,1._dp)

END MODULE basic_constants
