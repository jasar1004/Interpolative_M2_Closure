#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <limits>

#include "NG_MN_Model_3D_OPTIM.h"

using namespace std;

long double binomialCoefficients(const long double &n, const long double &k);
long double Jacobi_Polynomials(const long double &x, const int &l, const int &alpha, const int &beta);
long double Proriol_Polynomial_Basis(const int &k, const int &l, const long double &gam1, const long double &gam2);
int factorial(const int &n);
int factorial(const int &n, const int &a);
long double factorial_inv(const int &n, const int &a);
long double factorial_ratios(const int &num, const int &den);
int double_factorial(const int &n);
long double SH_Normalization_Constant(const int &l, const int &m);
long double P_l_m_Polynomials(const long double &x, const int &l, const int &m);
long double Associated_Legendre_Polynomials(const long double &x, const int &l, const int &m);

long double Precompute_Legendre_Basis_Coeffs(const int &n, const int &k);
long double Legendre_Polynomial(const long double &x, const int &n);

long double Precompute_First_Kind_Chebyshev_Basis_Coeffs(const int &n, const int &k);
void Chebyshev_First_Kind_to_Monomial_Basis_ratio_E(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &N_Coeffs_Triangle_gam1_gam2);
void Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &N_Coeffs_Triangle_gam1_gam2);
long double Precompute_SH_Basis_Coeffs(const int &l, const int &m, const int &i, const int &j, const int &k);
long double Cartesian_Spherical_harmonics(const long double &x, const long double &y, const long double &z, const int &l, const int &m);
long double Precompute_Proriol_Basis_Coeffs(const int &k, const int &l, const int &i, const int &j, const int &m, const int &n);
long double Proriol_Polynomials(const long double &zeta, const long double &eta, const int &k, const int &l);
void Proriol_to_Monomial_Basis(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &N_Coeffs_Triangle_gam1_gam2, const int &N_Coeffs_gam1);
void Spherical_Harmonics_to_Monomial_Basis(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &N_Coeffs_Triangle_gam1_gam2, const int &Order_SH, const int *Array_l_SH, const int *Array_m_SH);

int SH_Linear_Index(const int &Order_SH, const int &degree_SH);
int SH_Linear_Index_Symmetric(const int &Order_SH, const int &degree_SH);
void Triangle_to_Square_Mapping(long double &zeta, long double &eta, const long double &v_x, const long double &v_y);
long double VanderMonde_Matrix_1_var(const long double &var_1, const int &Order_poly_1);
long double VanderMonde_Matrix_2_vars(const long double &var_1, const long double &var_2, const int &Order_poly_1, const int &Order_poly_2);
long double VanderMonde_Matrix_3_vars(const long double &var_1, const long double &var_2, const long double &var_3, const int &Order_poly_1, const int &Order_poly_2, const int &Order_poly_3);
long double VanderMonde_Matrix_4_vars(const long double &var_1, const long double &var_2, const long double &var_3, const long double &var_4, const int &Order_poly_1, const int &Order_poly_2, const int &Order_poly_3, const int &Order_poly_4);
long double VanderMonde_Vector_N_vars(const int &Index_Entry, const int &Index_Point);
void Solve_A_x_b(const long double *VanderMonde_Matrix, long double *Coeff_Vand_Matrix, const long double *VanderMonde_Vector, const int &n, const int &index_Cheby);
void Check_A_x_b(const long double *VanderMonde_Matrix, const long double *Coeff_Vand_Matrix, const int &n);
void LUdecomposition(const long double *A, long double *L, long double *U, const int &n);
void LU_Solve_Back(const int &n, long double *x, const long double *A, const long double *b);
void LU_Solve_Forward(const int &n, long double *x, const long double *A, const long double *b);
long double Lebedev_Quadrature_Matrix_SH_Temp(const long double *Matrix, const int &Order_SH, const int &degree_SH, const int &Npts_Phi, const int &Npts_Theta, const long double *x_SH, const long double *y_SH, const long double *z_SH);
long double Lebedev_Quadrature_Matrix_SH(const long double *Matrix, const int &Order_SH, const int &degree_SH, const int &Npts_Phi, const int &Npts_Theta);
long double Lebedev_Quadrature_Orthog_SH(const int &order_SH1, const int &order_SH2, const int &m1, const int &m2, const int &Npts_Phi, const int &Npts_Theta);
long double Test_Spherical_Harmonics(const long double *Vals, const int &index_Cheby, const int &N_points_Cheby, const long double &x, const long double &y, const long double &z, const int &Order_SH);
long double Test_Spherical_Harmonics_Even_Odd(const long double *Vals, const int &index_Cheby, const int &N_points_Cheby, const long double &x, const long double &y, const long double &z, const int &Order_SH);
void Orthogonality_Test(const int &Quad_Rule);

long double Triangle_Quadrature_Orthog_Proriol(const int &i_gam1, const int &i_gam2, const int &i_prime_gam1, const int &i_prime_gam2, const int &N_points_Leg);
void Test_Proriol_Polynomials_Orthogonality(const int &Quad_Rule);

long double Cartesian_to_spherical_Coordinates_Conversion_Jacobian(const int &i, const int &j, const long double &r, const long double &mu, const long double &cos_phi, const long double &sin_phi);

long double Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(const int &i, const int &j, const long double &r, const long double &mu, const long double &cos_phi, const long double &sin_phi);

long double Cartesian_to_spherical_Coordinates_Jacobian(const long double &norm_f, const long double &x_val, const long double &y_val, const long double &z_val, const long double &dIn_Plus_dN1_1, const long double &dIn_Plus_dN1_2, const long double &dIn_Plus_dN1_3, const int &VAR_NUM);

long double Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(const long double &norm_f, const long double &x_val, const long double &y_val, const long double &z_val, const long double &dIn_Plus_dN1_1, const long double &dIn_Plus_dN1_2, const long double &dIn_Plus_dN1_3, const int &VAR_NUM);
