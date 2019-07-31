PROGRAM mckern
USE global_variables
USE numeric_kinds
USE basic_constants
USE spherical_harmonics
USE post_proc

IMPLICIT NONE

CALL main()

CONTAINS

SUBROUTINE main()
REAL (kind=dp)    :: mindist 

CALL potential()
CALL configuration()
CALL shortestlength(mindist)
global%nvals=int(potential_parameter%range/mindist)+1

CALL calc_bop()
print*,'   Q4_hat, Q6_hat, W4_hat, W6_hat'
WRITE(*,fmt='(A,F7.3,1X,F7.3,1X,F7.3,1X,F7.3,1X,A,F5.2,A)') '  ', &
   individual%q4hat,individual%q6hat,individual%w4hat,individual%w6hat, &
   ' (',individual%nn,')'

CALL close_lattice(lattice)
CALL close_global(global)

END SUBROUTINE main

SUBROUTINE potential()
INTEGER (kind=di) :: i, stat
REAL    (kind=dp) :: dummyreal
CHARACTER (kind=ch, len=60) :: dummystring

OPEN(9,file='mc-inputfile',status='old')
DO i=1,100
  READ(9,*,iostat=stat) dummystring, dummyreal
    IF (dummystring.EQ.'DIMENSION') THEN
      global%dimension = dummyreal
    END IF
    IF (dummystring.EQ.'NR_PATCHES') THEN
      global%nr_patches = dummyreal
    END IF
    IF (dummystring.EQ.'ATTRACTION_RANGE') THEN
      global%attraction_range = dummyreal
    END IF
END DO
CLOSE(9)

potential_parameter%sigma=1._dp                                                    
potential_parameter%rangesite=global%attraction_range*potential_parameter%sigma                                            
potential_parameter%range=potential_parameter%sigma+potential_parameter%rangesite

END SUBROUTINE potential

SUBROUTINE configuration()
INTEGER (kind=di) :: ip,jp
INTEGER (kind=di) :: nn,ntot
REAL (kind=dp) :: o1,o2,o3
REAL (kind=dp) :: ra1,ra2,ra3
REAL (kind=dp) :: rb1,rb2,rb3
REAL (kind=dp) :: rc1,rc2,rc3
REAL (kind=dp) :: cmx,cmy,cmz
REAL (kind=dp) :: rnx,rny,rnz
REAL (kind=dp) :: dummyVelX, dummyVelY, dummyVelZ
CHARACTER :: dummyType
REAL (kind=dp) :: dummyTime,simulationBoxSide

OPEN(3,file='startingstate.xyz',status='old')
 READ(3,*) nn
 global%nr_particles=nn/(global%nr_patches+1)
CLOSE(3)

global%ra => null()
global%rb => null()
global%rc => null()
lattice%cm => null()

ALLOCATE (global%rb(global%dimension))
ALLOCATE (global%ra(global%dimension))
ALLOCATE (global%rc(global%dimension))
ALLOCATE (lattice%cm(global%nr_particles,global%dimension))

OPEN(3,file='startingstate.xyz',status='old')
 READ(3,*) nn
 READ(3,*) simulationBoxSide,dummyTime
 global%ra(1)=simulationBoxSide
 global%ra(2)=0._dp
 global%ra(3)=0._dp
 global%rb(1)=0._dp
 global%rb(2)=simulationBoxSide
 global%rb(3)=0._dp
 global%rc(1)=0._dp
 global%rc(2)=0._dp
 global%rc(3)=simulationBoxSide

WRITE(*,*) global%nr_particles

 DO ip=1,global%nr_particles
    READ(3,*) dummyType,cmx,cmy,cmz,dummyVelX,dummyVelY,dummyVelZ
    lattice%cm(ip,1)=cmx*simulationBoxSide
    lattice%cm(ip,2)=cmy*simulationBoxSide
    lattice%cm(ip,3)=cmz*simulationBoxSide
    DO jp=1,global%nr_patches
        READ(3,*) dummyType,rnx,rny,rnz,dummyVelX,dummyVelY,dummyVelZ
    ENDDO

 ENDDO
CLOSE(3)

END SUBROUTINE configuration

SUBROUTINE shortestlength(short)
INTEGER (kind=di) :: i,kd
REAL (kind=dp) :: short,short1,short2,short3,short4
REAL (kind=dp),DIMENSION(:) :: la(global%dimension),lb(global%dimension),lc(global%dimension)
REAL (kind=dp) :: lam,lbm,lcm
REAL (kind=dp),DIMENSION(:) :: shortest(4)

global%boxa=0._dp
global%boxb=0._dp
global%boxc=0._dp
DO i=1,global%dimension
   global%boxa=global%boxa+global%ra(i)**2._dp
   global%boxb=global%boxb+global%rb(i)**2._dp
   global%boxc=global%boxc+global%rc(i)**2._dp
ENDDO
global%boxa=SQRT(global%boxa)
global%boxb=SQRT(global%boxb)
global%boxc=SQRT(global%boxc)

DO kd=1,global%dimension
   la(kd)=global%ra(kd)
   lb(kd)=global%rb(kd)
   lc(kd)=global%rc(kd)
ENDDO
lam=global%boxa
lbm=global%boxb
lcm=global%boxc

CALL shortestl2d(la,lb,lam,lbm,short1)
CALL shortestl2d(la,lc,lam,lcm,short2)
CALL shortestl2d(lb,lc,lbm,lcm,short3)

CALL shortestl3d(la,lb,lc,short4)

shortest(1)=short1
shortest(2)=short2
shortest(3)=short3
shortest(4)=short4

short=shortest(1)
DO i=2,4
   IF(shortest(i).LT.short) short=shortest(i)
ENDDO       

END SUBROUTINE shortestlength

SUBROUTINE shortestl2d(r1,r2,r1l,r2l,shortestl)
REAL (kind=dp), INTENT(IN) :: r1(global%dimension)
REAL (kind=dp), INTENT(IN) :: r2(global%dimension)
REAL (kind=dp), INTENT(IN):: r1l
REAL (kind=dp), INTENT(IN):: r2l
REAL (kind=dp), INTENT(OUT):: shortestl
REAL (kind=dp) :: theta,height1,height2
REAL (kind=dp) :: dotproduct12

dotproduct12=r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3)

theta=ACOS(dotproduct12/(r1l*r2l))
height1=r1l*SIN(theta)
height2=r2l*SIN(theta)

shortestl=height1
IF(height2.LT.shortestl) shortestl= height2

END SUBROUTINE shortestl2d

SUBROUTINE shortestl3d(r1,r2,r3,shortestl)
INTEGER (kind=di) :: i
REAL (kind=dp), INTENT(IN) :: r1(global%dimension)
REAL (kind=dp), INTENT(IN) :: r2(global%dimension)
REAL (kind=dp), INTENT(IN) :: r3(global%dimension)
REAL (kind=dp), INTENT(OUT):: shortestl
REAL (kind=dp) :: height1,height2,height3,d1
REAL (kind=dp), DIMENSION(:) :: normal1(global%dimension),normal2(global%dimension),normal3(global%dimension)
REAL (kind=dp), DIMENSION(:) :: crossproduct(global%dimension)

crossproduct(1)=r1(2)*r2(3)-r1(3)*r2(2)
crossproduct(2)=-1._dp*(r1(1)*r2(3)-r1(3)*r2(1))
crossproduct(3)=r1(1)*r2(2)-r1(2)*r2(1)
DO i=1,global%dimension
   normal1(i)=crossproduct(i)
ENDDO

crossproduct(1)=r1(2)*r3(3)-r1(3)*r3(2)
crossproduct(2)=-1._dp*(r1(1)*r3(3)-r1(3)*r3(1))
crossproduct(3)=r1(1)*r3(2)-r1(2)*r3(1)
DO i=1,global%dimension
   normal2(i)=crossproduct(i)
ENDDO

crossproduct(1)=r3(2)*r2(3)-r3(3)*r2(2)
crossproduct(2)=-1._dp*(r3(1)*r2(3)-r3(3)*r2(1))
crossproduct(3)=r3(1)*r2(2)-r3(2)*r2(1) 
DO i=1,global%dimension
   normal3(i)=crossproduct(i)
ENDDO

height1 = ABS(distancefromplane(normal1,r1,r3))
height2 = ABS(distancefromplane(normal2,r3,r2))
height3 = ABS(distancefromplane(normal3,r2,r1))

shortestl = height1
IF (height2.LT.shortestl) shortestl = height2
IF (height3.LT.shortestl) shortestl = height3

END SUBROUTINE shortestl3d

FUNCTION distancefromplane(nn,r0,rr)
REAL (kind=dp) :: distancefromplane
INTEGER (kind=di) :: i
REAL (kind=dp), INTENT(IN) :: rr(global%dimension)
REAL (kind=dp), INTENT(IN) :: r0(global%dimension)
REAL (kind=dp), INTENT(IN) :: nn(global%dimension)
REAL (kind=dp), DIMENSION(:) :: bb(global%dimension)
REAL (kind=dp) :: dotproductbn,dotproductnn

DO i=1,3
   bb(i)=rr(i)-r0(i)
ENDDO

dotproductbn=bb(1)*nn(1)+bb(2)*nn(2)+bb(3)*nn(3)
dotproductnn=nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3) 

distancefromplane= dotproductbn/SQRT(dotproductnn)

END FUNCTION distancefromplane

SUBROUTINE close_lattice(lattice)
TYPE (t_lattice) :: lattice

DEALLOCATE(lattice%cm)

END SUBROUTINE close_lattice

SUBROUTINE close_global(global)
TYPE (t_global) :: global

DEALLOCATE(global%ra)
DEALLOCATE(global%rb)
DEALLOCATE(global%rc)

END SUBROUTINE close_global

END PROGRAM
