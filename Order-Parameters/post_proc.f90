MODULE post_proc
 USE numeric_kinds
 USE basic_constants
 USE global_variables
 USE spherical_harmonics

IMPLICIT NONE

CONTAINS


SUBROUTINE calc_bop()
INTEGER (kind=di) :: m, NumNb, m1, m2, m3, ii, kk, jj, nn
COMPLEX (kind=dp), POINTER, DIMENSION(:,:) :: Q6m_bar, Q4m_bar
COMPLEX (kind=dp), POINTER, DIMENSION(:) :: W4, W6
REAL (kind=dp), POINTER, DIMENSION(:) :: W3j6, W3j4, Q6, Q4, W6_hat, W4_hat
REAL (kind=dp) :: r_q, r_q_fac, c, theta, phi, q4h, q6h, w4hut, w6hut
INTEGER (kind=di) :: i, j, na, nb, nc, n_max_c, estatus, nr_cells
REAL (kind=dp) :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, xj, yj, zj, rsq, r, histr(20), maxr
INTEGER (kind=di) :: neighbors(global%nr_particles,100), histogram(20), testv, testv2, histmax, histvalue(20), maxvalue
INTEGER (kind=di) :: tn, assigned(20)

! Structure analysis... bond order parameters

ax=global%ra(1)
ay=global%ra(2)
az=global%ra(3)
bx=global%rb(1)
by=global%rb(2)
bz=global%rb(3)
cx=global%rc(1)
cy=global%rc(2)
cz=global%rc(3)

ALLOCATE(Q6m_bar(1:global%nr_particles,-6:6), Q4m_bar(1:global%nr_particles,-4:4), W3j6(0:21), W3j4(0:10),&
   & Q6(1:global%nr_particles), Q4(1:global%nr_particles), W6(1:global%nr_particles), W4(1:global%nr_particles),&
   & W6_hat(1:global%nr_particles), W4_hat(1:global%nr_particles))

! Wigner 3j Symbols [(l,i)(l,j)(l,k)] 
!    where |i|>=|j|>=|k|

! ... for l == 6

! mapping: index = floor[(49*|i|+7*|j|+|k|)/8]
    
W3j6(0) = -0.093059500211290752118
W3j6(1) = 0.0465297501056453760589409285321
W3j6(7) = -0.0953576116969903179997867378635
W3j6(2) = 0.0511827251162099136648350213853
W3j6(8) = 0.0452320867794605134208355889657
W3j6(3) = -0.100038962727137558526722996344
W3j6(14) = -0.104459029912520160537187594731
W3j6(9)= 0.0688184284781714582592259403831
W3j6(4) = 0.0186119000422581504235763714128
W3j6(15) = 0.0408296900131293004050049089685
W3j6(10) = -0.106078646340041293368101671370
W3j6(5) = 0.127956812790524784162087553463
W3j6(21) = -0.141438195120055057824135561827
W3j6(16) = 0.129114816600118769292308582415
W3j6(11) = -0.0957541107531434804924219618366
W3j6(6) = 0.0511827251162099136648350213853
    
! ... for l == 4

! mapping: index = floor[(25*|i|+4*|j|+|k|)/6]

W3j4(0) = 0.134097046880302259764357335044
W3j4(1) = -0.0670485234401511298821786675222
W3j4(5) = 0.141350698548043900425691223634
W3j4(2) = -0.0819481953157402698559961491938
W3j4(6) = -0.0623297993338971451873100081774
W3j4(3) = 0.156446554693685969725083557552
W3j4(10) = 0.186989398001691435561930024532
W3j4(7) = -0.164909148306051217163306427573
W3j4(4) = 0.104297703129123979816722371701
    

! r_q defines "neighbourhood"
! check number of nearest neighbours

neighbors(:,:)=0


DO i=1,global%nr_particles
  DO nc=-global%nvals,global%nvals
    DO na=-global%nvals,global%nvals
      DO nb=-global%nvals,global%nvals
        DO j=1,global%nr_particles
          dx=lattice%cm(i,1)-(lattice%cm(j,1)+na*global%ra(1)+nb*global%rb(1)+nc*global%rc(1))
          dy=lattice%cm(i,2)-(lattice%cm(j,2)+na*global%ra(2)+nb*global%rb(2)+nc*global%rc(2))
          dz=lattice%cm(i,3)-(lattice%cm(j,3)+na*global%ra(3)+nb*global%rb(3)+nc*global%rc(3))
          rsq=dx**2+dy**2+dz**2

          DO testv=1,100
            r_q=1._dp+testv/100._dp
            IF((rsq.LE.r_q**2).and.(rsq.NE.0._dp)) THEN        
              neighbors(i,testv) = neighbors(i,testv) + 1
              DO testv2=testv+1,100
                neighbors(i,testv2) = neighbors(i,testv2) + 1
              END DO
              EXIT
            END IF
          END DO

        END DO
      END DO
    END DO
  END DO
END DO


histogram(:)=0
assigned(:)=0
histvalue(:)=0
histr(:)=0._dp

DO i=1,global%nr_particles
  DO testv=1,100
!    print*,'neighbors',i,testv,neighbors(i,testv)
    DO testv2=1,20
      IF (neighbors(i,testv).EQ.testv2) THEN
        IF (assigned(testv2).EQ.0) THEN
          tn=0
          DO j=1,global%nr_particles
            IF (neighbors(i,testv).NE.neighbors(j,testv)) THEN
              tn=1
            END IF
          END DO
          IF (tn.EQ.0) THEN
            histvalue(testv2)=neighbors(i,testv)
            histr(testv2)=1._dp+testv/100._dp
            assigned(testv2)=1
          END IF
        END IF
        histogram(testv2)=histogram(testv2)+1
      END IF
    END DO
  END DO
END DO


histmax=0
maxvalue=0
maxr=0._dp
DO testv=10,20
!  print*, 'hist',testv,histogram(testv)
  IF (histogram(testv).GE.histmax) THEN
    histmax=histogram(testv)
    maxvalue=histvalue(testv)
    maxr=histr(testv)
  END IF
END DO


IF (maxvalue.EQ.0) THEN

  histogram(:)=0
  assigned(:)=0
  histvalue(:)=0
  histr(:)=0._dp

  DO i=1,global%nr_particles
    DO testv=1,100
      DO testv2=1,20
        IF (neighbors(i,testv).EQ.testv2) THEN
          IF (assigned(testv2).EQ.0) THEN
            histvalue(testv2)=neighbors(i,testv)
            histr(testv2)=1._dp+testv/100._dp
            assigned(testv2)=1
          END IF
          histogram(testv2)=histogram(testv2)+1
        END IF
      END DO
    END DO
  END DO


  histmax=0
  maxvalue=0
  maxr=0._dp
  DO testv=1,20
    IF (histogram(testv).GE.histmax) THEN
      histmax=histogram(testv)
      maxvalue=histvalue(testv)
      maxr=histr(testv)
    END IF
  END DO

END IF


!print*,histmax,'times',maxvalue,'neighbors at r=',maxr
r_q =maxr 

! calculating starts...     
! 1) evaluate \bar{Q_lm} for all clusters
! Q6m_bar and Q4m_bar

Q6m_bar(:,:) = (0._dp,0._dp)
Q4m_bar(:,:) = (0._dp,0._dp)

individual%nn=0

DO i=1,global%nr_particles
! NumNb = number of neighbours
  NumNb = 0
  DO nc=-global%nvals,global%nvals
    DO na=-global%nvals,global%nvals
      DO nb=-global%nvals,global%nvals
        DO j=1,global%nr_particles
          dx=lattice%cm(i,1)-(lattice%cm(j,1)+na*global%ra(1)+nb*global%rb(1)+nc*global%rc(1))
          dy=lattice%cm(i,2)-(lattice%cm(j,2)+na*global%ra(2)+nb*global%rb(2)+nc*global%rc(2))
          dz=lattice%cm(i,3)-(lattice%cm(j,3)+na*global%ra(3)+nb*global%rb(3)+nc*global%rc(3))
          rsq=dx**2+dy**2+dz**2

          IF((rsq.LE.r_q**2).and.(rsq.NE.0._dp)) THEN
! calculate \theta and \phi for r_ij
            theta = ACOS(dz/SQRT(dx**2+dy**2+dz**2))
            phi = ATAN2(dy,dx)

!            if (m.eq.0) then
!              print*,i,sqrt(rsq)
!              print*,nc,na,nb,j
!              print*,dx,dy,dz
!              print*,theta,phi
!            end if     

            NumNb = NumNb+1

            DO m=-6,6

! spherical harmonics
              Q6m_bar(i,m) = Q6m_bar(i,m) + Y6m(m,theta,phi)
              IF(ABS(m).LE.4) THEN
                Q4m_bar(i,m) = Q4m_bar(i,m) + Y4m(m,theta,phi)
              END IF

            END DO

          END IF

        END DO
      END DO
    END DO
  END DO

!  print*,r_q, i, 'neighbours: ',NumNb

  individual%nn=individual%nn+NumNb

  IF(NumNb.EQ.0) THEN
!    PRINT*, "NumNb=0, increase r_q", r_q
!    DO testv=1,20
!      print*, 'hist',testv,histogram(testv)
!    END DO
!    print*, 'histmax,maxvalue,maxr',histmax,maxvalue,maxr
    RETURN
  END IF


  DO m=-6,6

    Q6m_bar(i,m) = Q6m_bar(i,m)/REAL(NumNb,dp)
    IF(ABS(m).LE.4) THEN
       Q4m_bar(i,m) = Q4m_bar(i,m)/REAL(NumNb,dp)
    END IF
  END DO

!print*,m,'Q6m_bar',Q6m_bar(i,m)


END DO

individual%nn=(1._dp*individual%nn)/global%nr_particles
    
    
    DO i=1,global%nr_particles
      Q6(i) = 0._dp
      Q4(i) = 0._dp
      DO m=-6,6
        Q6(i) = Q6(i) + (ABS(Q6m_bar(i,m)))**2
        IF(ABS(m).LE.4) THEN
          Q4(i) = Q4(i) + (ABS(Q4m_bar(i,m)))**2
        END IF
      END DO
      Q6(i) = SQRT(Q6(i)*4._dp*pi/13._dp)
      Q4(i) = SQRT(Q4(i)*4._dp*pi/9._dp)
    END DO
    
    
    
! W6 
    
    DO i=1,global%nr_particles
      W6(i) = (0._dp,0._dp)
      DO m1=-6,6
        DO m2=-6,6
          DO m3=-6,6
            IF((m1+m2+m3).EQ.0) THEN
              ii = ABS(MAX(ABS(m1),ABS(m2),ABS(m3)))
              kk = ABS(MIN(ABS(m1),ABS(m2), ABS(m3)))
              jj = ABS(ABS(m1)+ABS(m2)+ABS(m3)-ii-kk)
              nn = FLOOR(((49._dp*REAL(kk,dp)+7._dp*REAL(jj,dp)+REAL(ii,dp))/8._dp))
              W6(i) = W6(i) + W3j6(nn)*Q6m_bar(i,m1)*Q6m_bar(i,m2)*Q6m_bar(i,m3)
            END IF
          END DO
        END DO
      END DO
    END DO
    
! W4 
    
    DO i=1,global%nr_particles
      W4(i) = (0._dp,0._dp)
      DO m1=-4,4
        DO m2=-4,4
          DO m3=-4,4
            IF((m1+m2+m3).EQ.0) THEN
              ii = ABS(MAX(ABS(m1),ABS(m2),ABS(m3)))
              kk = ABS(MIN(ABS(m1),ABS(m2), ABS(m3)))
              jj = ABS(ABS(m1)+ABS(m2)+ABS(m3)-ii-kk)
              nn = FLOOR(((25._dp*REAL(kk,dp)+5._dp*REAL(jj,dp)+REAL(ii,dp))/6._dp))
              W4(i) = W4(i) + W3j4(nn)*Q4m_bar(i,m1)*Q4m_bar(i,m2)*Q4m_bar(i,m3)
            END IF
          END DO
        END DO
      END DO
    END DO
    
    DO i = 1,global%nr_particles
      W6_hat(i) = REAL(W6(i)/(Q6(i)**3*(13._dp/(4._dp*pi))**(3._dp/2._dp)))
      W4_hat(i) = REAL(W4(i)/(Q4(i)**3*(9._dp/(4._dp*pi))**(3._dp/2._dp)))
    END DO


w6hut=0._dp
w4hut=0._dp
q4h=0._dp
q6h=0._dp
do i=1,global%nr_particles
  w6hut=w6hut+W6_hat(i)
  w4hut=w4hut+W4_hat(i)
  q4h=q4h+Q4(i)
  q6h=q6h+Q6(i)
end do

w6hut=w6hut/global%nr_particles
w4hut=w4hut/global%nr_particles
q4h=q4h/global%nr_particles
q6h=q6h/global%nr_particles

individual%w6hat=w6hut
individual%w4hat=w4hut
individual%q6hat=q6h
individual%q4hat=q4h
    
    
DEALLOCATE(Q6,Q4,Q6m_bar,Q4m_bar,W3j6,W3j4, W4, W6, W4_hat, W6_hat)

END SUBROUTINE calc_bop

SUBROUTINE compare_bop(individual1,individual2,same)
TYPE(t_individual) :: individual1, individual2
LOGICAL :: same
REAL (kind=dp) :: tolq4, tolq6, tolw4, tolw6

tolq4=0.1_dp
tolq6=0.1_dp
tolw4=0.01_dp
tolw6=0.01_dp

same=.FALSE.

IF ((individual1%q4hat.GE.individual2%q4hat-tolq4).AND.(individual1%q4hat.LE.individual2%q4hat+tolq4)) THEN
  IF((individual1%q6hat.GE.individual2%q6hat-tolq6).AND.(individual1%q6hat.LE.individual2%q6hat+tolq6)) THEN
    IF((individual1%w4hat.GE.individual2%w4hat-tolw4).AND.(individual1%w4hat.LE.individual2%w4hat+tolw4)) THEN
      IF((individual1%w6hat.GE.individual2%w6hat-tolw6).AND.(individual1%w6hat.LE.individual2%w6hat+tolw6)) THEN
        same=.TRUE.
      END IF
    END IF
  END IF
END IF

!print*,same
!WRITE(*,fmt='(F5.2,1X,F5.2,1X,F5.2,1X,F5.2)') individual1%q4hat-individual2%q4hat,individual1%q6hat-individual2%q6hat,individual1%w4hat-individual2%w4hat,individual1%w6hat-individual2%w6hat

END SUBROUTINE compare_bop


END MODULE
