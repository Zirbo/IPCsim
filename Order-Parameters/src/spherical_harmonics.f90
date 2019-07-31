MODULE spherical_harmonics 
   
USE numeric_kinds
USE basic_constants
 
IMPLICIT NONE
   
  PRIVATE :: faktorial, LegendreP4, LegendreP6 !gser, gcf, get_gamma, get_gamma_p 
  PUBLIC :: Y6m, Y4m, LegendreP !gamma_ln, gamma_q, gamma_p, inc_gamma, gamma 
 
 
CONTAINS 
   
!************************************************************************! 
  FUNCTION Y6m(m, theta, phi) 
!************************************************************************! 
    COMPLEX(kind=dp) :: Y6m 
    REAL(kind=dp), INTENT(in) :: theta, phi 
    INTEGER(kind=di) :: m 
! 
!    IF(theta.EQ.0._dp) THEN 
!      Y6m = (0._dp, 0._dp) 
!    ELSE    
    IF(m.GE.0) THEN 
      Y6m = SQRT(13._dp/(4._dp*pi)*faktorial(6-m)/faktorial(6+m))* & 
           & LegendreP(6,m,COS(theta)) & 
           & *EXP(imag * CMPLX(REAL(m,dp),0._dp,dc) * CMPLX(phi,0._dp,dc)) 
       
    ELSE 
      Y6m = (-1)**ABS(m)* SQRT(13._dp/(4._dp*pi)*faktorial(6+m)/faktorial(6-m))*& 
           & LegendreP(6,-m,COS(theta)) & 
           & *EXP(-imag * CMPLX(REAL(-m,dp),0._dp,dc) * CMPLX(phi,0._dp,dc)) 
    END IF 
!    END IF 
     
  END FUNCTION Y6m 
 
!************************************************************************! 
  FUNCTION Y4m(m, theta, phi) 
!************************************************************************! 
    COMPLEX(kind=dp) :: Y4m 
    REAL(kind=dp), INTENT(in) :: theta, phi 
    INTEGER(kind=di) :: m 
 
!    IF(theta.EQ.0._dp) THEN 
!      Y4m = (0._dp, 0._dp) 
!    ELSE  
    IF(m.GE.0) THEN 
      Y4m = SQRT(9._dp/(4._dp*pi)*faktorial(4-m)/faktorial(4+m))* & 
           & LegendreP(4,m,COS(theta)) & 
           & *EXP(imag * CMPLX(REAL(m,dp),0._dp,dc) * CMPLX(phi,0._dp,dc)) 
       
    ELSE 
      Y4m = (-1)**ABS(m)* SQRT(9._dp/(4._dp*pi)*faktorial(4+m)/faktorial(4-m))*& 
           & LegendreP(4,-m,COS(theta)) & 
           & *EXP(-imag * CMPLX(REAL(-m,dp),0._dp,dc) * CMPLX(phi,0._dp,dc)) 
    END IF 
!    END IF 
  END FUNCTION Y4m 
 
!************************************************************************! 
  RECURSIVE FUNCTION faktorial(n) RESULT (res) 
!************************************************************************! 
    REAL(kind=dp) :: res 
    INTEGER(kind=di), INTENT(in) :: n 
     
    IF(n.GT.22) THEN  
      PRINT*, "legendre_polynomials reports: function fak called with & 
           & argument > 22." 
      STOP 
    END IF 
    IF(n.EQ.0) THEN 
      res = 1._dp 
      return 
    END IF 
 
    IF(n.NE.1) THEN 
      res = faktorial(n-1)*REAL(n,dp) 
    ELSE 
      res = 1._dp 
    END IF 
     
  END FUNCTION faktorial 
 
 
!************************************************************************! 
  FUNCTION LegendreP(l,m,x) 
!************************************************************************! 
    INTEGER(kind=di), INTENT(in):: l,m 
    REAL(kind=dp) :: LegendreP, x 
    INTEGER(kind=di) :: i, ll 
    REAL(kind=dp) :: fact, pll, pmm, pmmp1, somx2 
    IF(m.LT.0.OR.m.GT.l.OR.ABS(x).GT.1) THEN 
      PRINT*, "LegendreP called with wrong arguments" 
      STOP 
    END IF 
 
    pmm=1._dp 
    IF(m.GT.0) THEN 
      somx2 = SQRT((1._dp-x)*(1._dp+x)) 
      fact=1._dp 
      DO i = 1,m 
        pmm=-pmm*fact*somx2 
        fact = fact + 2._dp 
      END DO 
    END IF 
    IF(l.EQ.m) THEN 
      LegendreP = pmm 
    ELSE 
      pmmp1 = x*(2._dp*REAL(m,dp)+1._dp)*pmm 
      IF(l.EQ.m+1) THEN 
        LegendreP = pmmp1 
      ELSE 
        DO ll = m+2,l 
          pll = (x*(2._dp*REAL(ll,dp)-1._dp)*pmmp1-& 
               & (REAL(ll+m,dp)-1._dp)*pmm)/REAL(ll-m,dp) 
          pmm = pmmp1 
          pmmp1= pll 
        END DO 
        LegendreP = pll 
      END IF 
    END IF 
    RETURN 
  END FUNCTION LegendreP 
!************************************************************************! 
  FUNCTION LegendreP6(m,x) 
!************************************************************************! 
    REAL (KIND=dp) :: LegendreP6 
    REAL (KIND=dp), INTENT(IN) :: x 
    INTEGER (KIND=di) :: m 
     
    SELECT CASE (m) 
    CASE(0) 
      LegendreP6 = (-5._dp + 105._dp*x**2 - 315._dp*x**4 + 231._dp*x**6)/16._dp 
    CASE(1) 
      LegendreP6 = 21._dp/8._dp*SQRT(x-1._dp)*SQRT(x+1._dp)*x*& 
           & (33._dp*x**4-30._dp*x**2+5._dp)  
!21._dp/8._dp*SQRT((-1._dp-x)/(-1._dp+x))*(x-1._dp)* & 
!      &(5._dp*x - 30._dp*x**3 + 33._dp*x**5) 
    CASE(2) 
      LegendreP6 = 105._dp/8._dp*(x**2-1)*(1._dp + 33._dp*x**4 - 18._dp*x**2) 
!-105._dp/8._dp*(x-1._dp)*(x+1._dp)*& 
! & (1._dp - 18._dp*x**2 + 33._dp*x**4) 
    CASE(3) 
      LegendreP6 = 315._dp/2._dp*(x-1._dp)**(3._dp/2._dp)*(x+1._dp)**(3._dp/2._dp)*& 
           & x*(11._dp*x**2-3._dp) 
! -315._dp/2._dp*SQRT((-1._dp-x)/(-1._dp+x))*(x-1._dp)**2*& 
! & (x+1._dp)*(-3._dp*x + 11._dp*x**3) 
    CASE(4) 
      LegendreP6 = 945._dp/2._dp*(x-1._dp)**2*(x+1._dp)**2*(-1._dp+11._dp*x**2) 
! 945._dp/2._dp*(x-1._dp)**2*(x+1._dp)**2*& 
!  & (-1._dp + 11._dp*x**2) 
    CASE(5) 
      LegendreP6 = 10395._dp*(x-1._dp)**(5._dp/2._dp)*(x+1._dp)**(5._dp/2._dp)*x 
! 10395._dp*SQRT((-1._dp - x)/(-1._dp + x))*(-1._dp + x)**3*& 
! & x*(1._dp + x)**2 
    CASE(6)  
      LegendreP6 = 10395._dp*(x-1._dp)**3*(x+1._dp)**3 
! -10395._dp*(-1._dp + x)**3*(1._dp + x)**3 
    CASE default 
      PRINT*, "legendre_polynomials reports: function LegendreP6 & 
          & called with wrong argument." 
      STOP 
    END SELECT 
     
  END FUNCTION LegendreP6 
 
!************************************************************************! 
  FUNCTION LegendreP4(m,x) 
!************************************************************************! 
    REAL (KIND=dp) :: LegendreP4 
    REAL (KIND=dp), INTENT(IN) :: x 
    INTEGER (KIND=di) :: m 
     
    SELECT CASE (m) 
    CASE(0) 
      LegendreP4 = 3._dp/8._dp - (15._dp*x**2)/4._dp + (35._dp*x**4)/8._dp 
    CASE(1) 
      LegendreP4 = 5._dp/2._dp * SQRT(x-1._dp)*SQRT(x+1._dp)*x*& 
           & (7._dp*x**2-3._dp) 
! (5._dp*SQRT((-1._dp - x)/(-1._dp + x))*(-1._dp + x)*& 
! &(-3._dp*x + 7._dp*x**3))/2._dp 
    CASE(2) 
      LegendreP4 = 15._dp/2._dp*(x**2-1._dp)*(-1._dp+7._dp*x**2) 
! (-15._dp*(-1._dp + x)*(1._dp + x)*& 
!  & (-1._dp + 7._dp*x**2))/2._dp 
    CASE(3) 
      LegendreP4 = 105._dp*(x-1._dp)**(3._dp/2._dp)*(x+1._dp)**(3._dp/2._dp)*x 
! -105._dp*Sqrt((-1._dp - x)/(-1._dp + x))*& 
! & (-1._dp + x)**2*x*(1._dp + x) 
    CASE(4) 
      LegendreP4 = 105._dp*(x-1._dp)**2*(x+1._dp)**2 
! 105._dp*(-1._dp + x)**2._dp*(1._dp + x)**2 
    CASE default 
      PRINT*, "legendre_polynomials reports: function LegendreP4 & 
           & called with wrong argument." 
      STOP 
    END SELECT 
     
  END FUNCTION LegendreP4 
 
END MODULE spherical_harmonics 

