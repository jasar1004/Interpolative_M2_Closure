OTHERDIR_MN_OPT_ALGO=./
OTHERDIR_MN_INTERP=../M2_Model_Interp
OTHERDIR_PACKAGES=../../Packages

LFLAGS = -Wall -Wextra -g -liomp5 -lpthread

OBJS_M2_OPT = $(OTHERDIR_MN_OPT_ALGO)/NG_MN_Model_3D_OPTIM.o $(OTHERDIR_MN_OPT_ALGO)/Third_Order_Closing_Fluxes.o $(OTHERDIR_MN_OPT_ALGO)/Triangle_Computations.o $(OTHERDIR_MN_OPT_ALGO)/MPI_Utilities.o $(OTHERDIR_MN_OPT_ALGO)/Derivatives_Angular_Moments.o  $(OTHERDIR_MN_OPT_ALGO)/M2_Model_3D_Utilities.o $(OTHERDIR_MN_OPT_ALGO)/circle_rule.o

OBJS_M2_INTERP = $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby.o $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby_Least_Squares_L_I0_Star_N3_ijk.o $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby_Least_Squares_L_norm_f_N3_ijk.o

OBJS_OTHERDIR_PACKAGES = $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.o $(OTHERDIR_PACKAGES)/Finite_Difference/differ.o

Adaptive_Quadrature_Test: Adaptive_Quadrature_Test.o $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES)
		mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Adaptive_Quadrature_Test.cc -L $(CFFC_Path)/nlopt/lib -o Adaptive_Quadrature_Test $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm

Adaptive_Quadrature_Test.o: Adaptive_Quadrature_Test.cc
		mpicc $(LFLAGS) -c -o Adaptive_Quadrature_Test.o Adaptive_Quadrature_Test.cc
		
NG_MN_Model_3D_OPTIM.o: NG_MN_Model_3D_OPTIM.cc
	mpicc $(LFLAGS) -c NG_MN_Model_3D_OPTIM.cc
	
Derivatives_Angular_Moments.o: Derivatives_Angular_Moments.cc
	mpicc $(LFLAGS) -c Derivatives_Angular_Moments.cc
	
Derivatives_Angular_Moments_Finite_Difference.o: Derivatives_Angular_Moments_Finite_Difference.cc
	mpicc $(LFLAGS) -c Derivatives_Angular_Moments_Finite_Difference.cc
	
NG_MN_Model_3D_OPTIM_ClOSING_FLUX.o: NG_MN_Model_3D_OPTIM_ClOSING_FLUX.cc
	mpicc -c NG_MN_Model_3D_OPTIM_ClOSING_FLUX.cc
	
M2_Model_3D_Utilities.o: M2_Model_3D_Utilities.cc
	mpicc -c M2_Model_3D_Utilities.cc
	
chebyshev.o: $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.cpp
	
sphere_lebedev_rule.o: $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.cpp
		
circle_rule.o: circle_rule.cpp
	mpicc $(LFLAGS) -c circle_rule.cpp
		
quadrule.o: $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.cpp
	
differ.o: $(OTHERDIR_PACKAGES)/Finite_Difference/differ.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Finite_Difference/differ.cpp

clean:
	rm -rf *.o
	rm -rf *.o
	rm -rf $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Chebyshev/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Permutations/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Finite_Difference/*.o
