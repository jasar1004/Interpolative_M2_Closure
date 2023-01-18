OTHERDIR_MN_OPT_ALGO=../Non_Gray_Optimization
OTHERDIR_MN_INTERP=../M2_Model_Interp
OTHERDIR_PACKAGES=../../Packages

LFLAGS = -Wall -Wextra -g -liomp5 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lifcore -lpthread

OBJS_M2_OPT = $(OTHERDIR_MN_OPT_ALGO)/NG_MN_Model_3D_OPTIM.o $(OTHERDIR_MN_OPT_ALGO)/Third_Order_Closing_Fluxes.o $(OTHERDIR_MN_OPT_ALGO)/Triangle_Computations.o $(OTHERDIR_MN_OPT_ALGO)/MPI_Utilities.o $(OTHERDIR_MN_OPT_ALGO)/Derivatives_Angular_Moments.o  $(OTHERDIR_MN_OPT_ALGO)/M2_Model_3D_Utilities.o $(OTHERDIR_MN_OPT_ALGO)/circle_rule.o

OBJS_M2_INTERP = $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby.o $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby_Least_Squares_L_I0_Star_N3_ijk.o $(OTHERDIR_MN_INTERP)/Derivatives_M2_Closing_Fluxes_Fit.o $(OTHERDIR_MN_INTERP)/M2_Closure_Fit_Eigenstructure.o $(OTHERDIR_MN_INTERP)/M2_Closure_Rotations_Utilities.o

OBJS_OTHERDIR_PACKAGES = $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.o $(OTHERDIR_PACKAGES)/Finite_Difference/differ.o

Uniform_Non_Gray_M2_3D_L_N3_ijk: Uniform_Non_Gray_M2_3D_L_N3_ijk.o $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES)
		mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Uniform_Non_Gray_M2_3D_L_N3_ijk.cc -L $(CFFC_Path)/nlopt/lib -o Uniform_Non_Gray_M2_3D_L_N3_ijk $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm

Uniform_Non_Gray_M2_3D_HL_LL: Uniform_Non_Gray_M2_3D_HL_LL.o $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES)
		mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Uniform_Non_Gray_M2_3D_HL_LL.cc -L $(CFFC_Path)/nlopt/lib -o Uniform_Non_Gray_M2_3D_HL_LL $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm

Uniform_Non_Gray_M2_N1_1D_L_N3_ijk: Uniform_Non_Gray_M2_N1_1D_L_N3_ijk.o $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES)
		mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Uniform_Non_Gray_M2_N1_1D_L_N3_ijk.cc -L $(CFFC_Path)/nlopt/lib -o Uniform_Non_Gray_M2_N1_1D_L_N3_ijk $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm

Illustration_Energy_Density: Illustration_Energy_Density.o $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES)
		mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Illustration_Energy_Density.cc -L $(CFFC_Path)/nlopt/lib -o Illustration_Energy_Density $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm

Illustration_Length_Scale: Illustration_Length_Scale.o $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES)
		mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Illustration_Length_Scale.cc -L $(CFFC_Path)/nlopt/lib -o Illustration_Length_Scale $(OBJS_M2_OPT) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
			
Polynomial_Interpolation_Non_Gray_M2_Closure: Polynomial_Interpolation_Non_Gray_M2_Closure.o $(OBJS_M2_OPT) $(OBJS_M2_INTERP) $(OBJS_OTHERDIR_PACKAGES)
		mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Polynomial_Interpolation_Non_Gray_M2_Closure.cc -L $(CFFC_Path)/nlopt/lib -o Polynomial_Interpolation_Non_Gray_M2_Closure $(OBJS_M2_OPT) $(OBJS_M2_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
			
Tests_M2_Interpolation: Tests_M2_Interpolation.o $(OBJS_M2_OPT) $(OBJS_M2_INTERP) $(OBJS_OTHERDIR_PACKAGES)
		mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Tests_M2_Interpolation.cc -L $(CFFC_Path)/nlopt/lib -o Tests_M2_Interpolation $(OBJS_M2_OPT) $(OBJS_M2_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
		
Polynomial_Interpolation_Non_Gray_M2_Closure.o: Polynomial_Interpolation_Non_Gray_M2_Closure.cc
		mpicc $(LFLAGS) -c -o Polynomial_Interpolation_Non_Gray_M2_Closure.o Polynomial_Interpolation_Non_Gray_M2_Closure.cc
		
Tests_M2_Interpolation.o: Tests_M2_Interpolation.cc
		mpicc $(LFLAGS) -c -o Tests_M2_Interpolation.o Tests_M2_Interpolation.cc
	       
Uniform_Non_Gray_M2_3D_L_N3_ijk.o: Uniform_Non_Gray_M2_3D_L_N3_ijk.cc
		mpicc $(LFLAGS) -c -o Uniform_Non_Gray_M2_3D_L_N3_ijk.o Uniform_Non_Gray_M2_3D_L_N3_ijk.cc
	       
Uniform_Non_Gray_M2_3D_HL_LL.o: Uniform_Non_Gray_M2_3D_HL_LL.cc
		mpicc $(LFLAGS) -c -o Uniform_Non_Gray_M2_3D_HL_LL.o Uniform_Non_Gray_M2_3D_HL_LL.cc
		
Uniform_Non_Gray_M2_N1_1D_L_N3_ijk.o: Uniform_Non_Gray_M2_N1_1D_L_N3_ijk.cc
		mpicc $(LFLAGS) -c -o Uniform_Non_Gray_M2_N1_1D_L_N3_ijk.o Uniform_Non_Gray_M2_N1_1D_L_N3_ijk.cc
	       
Illustration_Energy_Density.o: Illustration_Energy_Density.cc
		mpicc $(LFLAGS) -c -o Illustration_Energy_Density.o Illustration_Energy_Density.cc
	       
Illustration_Length_Scale.o: Illustration_Length_Scale.cc
		mpicc $(LFLAGS) -c -o Illustration_Length_Scale.o Illustration_Length_Scale.cc
		
Chebyshev_Grid_M2_RT.o: Chebyshev_Grid_M2_RT.cc
		mpicc $(LFLAGS) -c -o Chebyshev_Grid_M2_RT.o Chebyshev_Grid_M2_RT.cc
		
# Objects for OBJS_M2_INTERP
M2_Closure_Rotations_Utilities.o: $(OTHERDIR_MN_INTERP)/M2_Closure_Fit_Eigenstructure.cc
		mpicc $(LFLAGS) -c -o $(OTHERDIR_MN_INTERP)/M2_Closure_Rotations_Utilities.o $(OTHERDIR_MN_INTERP)/M2_Closure_Rotations_Utilities.cc

M2_Closure_Fit_Eigenstructure.o: $(OTHERDIR_MN_INTERP)/M2_Closure_Fit_Eigenstructure.cc
		mpicc $(LFLAGS) -c -o $(OTHERDIR_MN_INTERP)/M2_Closure_Fit_Eigenstructure.o $(OTHERDIR_MN_INTERP)/M2_Closure_Fit_Eigenstructure.cc

Derivatives_M2_Closing_Fluxes_Fit.o: $(OTHERDIR_MN_INTERP)/Derivatives_M2_Closing_Fluxes_Fit.cc
		mpicc $(LFLAGS) -c -o $(OTHERDIR_MN_INTERP)/Derivatives_M2_Closing_Fluxes_Fit.o $(OTHERDIR_MN_INTERP)/Derivatives_M2_Closing_Fluxes_Fit.cc

N3_Non_Gray_M2_3D_RT_Cheby_Least_Squares_L_I0_Star_N3_ijk.o: $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby_Least_Squares_L_I0_Star_N3_ijk.cc
		mpicc $(LFLAGS) -c -o $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby_Least_Squares_L_I0_Star_N3_ijk.o $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby_Least_Squares_L_I0_Star_N3_ijk.cc
		
N3_Non_Gray_M2_3D_RT_Cheby.o: $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby.cc
		mpicc $(LFLAGS) -c $(OTHERDIR_MN_INTERP)/N3_Non_Gray_M2_3D_RT_Cheby.cc
		
# Objects for OBJS_M2_OPT
NG_MN_Model_3D_OPTIM.o: NG_MN_Model_3D_OPTIM.cc
	g++ $(LFLAGS) -c $(OTHERDIR_MN_OPT_ALGO)/NG_MN_Model_3D_OPTIM.cc
	
Third_Order_Closing_Fluxes.o: Third_Order_Closing_Fluxes.cc
	g++ -c $(OTHERDIR_MN_OPT_ALGO)/Third_Order_Closing_Fluxes.cc
	
Triangle_Computations.o: $(OTHERDIR_MN_OPT_ALGO)/Triangle_Computations.cc
	g++ $(LFLAGS) -c $(OTHERDIR_MN_OPT_ALGO)/Triangle_Computations.cc
	
MPI_Utilities.o: $(OTHERDIR_MN_OPT_ALGO)/MPI_Utilities.cc
	g++ $(LFLAGS) -c $(OTHERDIR_MN_OPT_ALGO)/MPI_Utilities.cc
	
Derivatives_Angular_Moments.o: $(OTHERDIR_MN_OPT_ALGO)/Derivatives_Angular_Moments.cc
	g++ $(LFLAGS) -c $(OTHERDIR_MN_OPT_ALGO)/Derivatives_Angular_Moments.cc
	
M2_Model_3D_Utilities.o: M2_Model_3D_Utilities.cc
	g++ -c $(OTHERDIR_MN_OPT_ALGO)/M2_Model_3D_Utilities.cc
		
circle_rule.o: $(OTHERDIR_MN_OPT_ALGO)/circle_rule.cpp
	g++ $(LFLAGS) -c $(OTHERDIR_MN_OPT_ALGO)/circle_rule.cpp
	
# Objects for OBJS_OTHERDIR_PACKAGES
Permutations_With_Order.o: $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.cc
	       mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.cc
	       
Permutations_With_Order_Tests.o: $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order_Tests.cc
	       mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order_Tests.cc
	       
chebyshev.o: $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.cpp
	g++ $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.cpp
	
sphere_lebedev_rule.o: $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.cpp
	g++ $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.cpp
		
quadrule.o: $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.cpp
	
differ.o: $(OTHERDIR_PACKAGES)/Finite_Difference/differ.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Finite_Difference/differ.cpp

clean:
	rm -rf *.o
	rm -rf $(OTHERDIR_MN_OPT_ALGO)/*.o
	rm -rf $(OTHERDIR_MN_INTERP)/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Chebyshev/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/*.o
	rm -rf $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Permutations/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Finite_Difference/*.o
