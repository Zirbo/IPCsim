MODULE global_variables
USE numeric_kinds

TYPE t_global
  INTEGER (kind=di) :: nvals
  INTEGER (kind=di) :: nr_particles,nr_patches,dimension
  REAL (kind=dp)    :: attraction_range
  REAL (kind=dp),DIMENSION(:),POINTER :: ra,rb,rc
  REAL (kind=dp) :: boxa,boxb,boxc
END TYPE t_global
TYPE (t_global) :: global

TYPE t_potential_parameter
  REAL (kind=dp) :: sigma, rangesite, range
END TYPE t_potential_parameter
TYPE (t_potential_parameter) :: potential_parameter

TYPE t_lattice
  REAL (kind=dp),DIMENSION(:,:),POINTER :: cm               
END TYPE t_lattice
TYPE (t_lattice) :: lattice

TYPE t_individual
  REAL (kind=dp) ::  w6hat, w4hat, q6hat, q4hat, nn
END TYPE t_individual
TYPE (t_individual) :: individual


END MODULE global_variables
