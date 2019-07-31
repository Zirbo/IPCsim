### compiler options: ###########################################################################
### optimisation: **************************
#	-O<n>			optimisation, default=2
#	-arch <processor>	pn1, pn2, pn3, pn4, SSE, SSE1
#	-ax<processor>		directs the compiler to generate processor-specific code
#	-x<processor>		directs the compiler to generate specialized code for the Intel® processor
#	-mtune=<processor>	performs optimizations for a specified CPU
### specify output form: *******************
#	-static 		including all libraries in executable
#	-c			create object file *.o (compile only)
#	-shared			produce a dynamic shared object (DSO)
### debugging: *****************************
#	-g			compile for debugger
#	-debug extended		additionally to -g
#	-warn [no]interfaces	check interfaces
#	-syntax-only		check only for correct syntax
#	-check all		checks run-time conditions
#	-warn all		Enables all warning messages
#	-stand<keyword>		issue compile-time messages for nonstandard language elements (f90, f95)
#	-traceback		generate extra information when a severe error occurs at run time
#	-ftrapuv		initializes stack local variables to an unusual value (error detection)
#	-zero			initializes to zero all local scalar variables (int., real, compl., log.)
### numerical precision: *******************
#	-integer_size <size>	set default integer size
#	-real_size <size>	specifies the default size (in bits) for real and complex variables
#		32	default
#		64	double precision
#		128	extended precision
#	-double_size <size>	defines the default KIND for DOUBLE PRECISION and DOUBLE COMPLEX (64, 128)
#		-real_size 64 is equal to DOUBLE PRECISION without flags
#	-pc<n>			enables control of floating-point significand precision (32, 64, 80)
#	-prec-sqrt		improves precision of square root implementations
#	-fltconsistency		enables improved floating-point consistency and reduce execution speed
### other: *********************************
#	-fpp			Runs the Fortran preprocessor on source files before compilation
#	-ssp			enables Software-based Speculative Pre-computation (SSP) optimization
#	-save			saves all variables in static allocation
#	-recursive		tells the compiler to compile for possible recursive execution (-automatic)
##################################################################################################

fc          = gfortran
optflags    = -O3 # -static -traceback -check all
gc          = cc

#LIBS    = -L$(HOME)/libs/fortran -lfortran

modules_patchy= numeric_kinds.o basic_constants.o global_variables.o spherical_harmonics.o post_proc.o crystal_analysis.o

patchy:  $(modules_patchy)
	$(fc) $(optflags) $(modules_patchy) -o crystal_analysis

%.o: %.c
	$(gc) $< -c -o $@

%.o: %.f90
	$(fc) $(optflags) $< -c -o $@

%.o: %.f
	$(fc) $(optflags) $< -c -o $@

clean:
	rm -f *.o *.mod *~ 
