#ifndef _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED
#define _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED

/*******************************************************************
  File: N3_Non_Gray_M2_3D_RT_Cheby.h

  Description:  ...  

  Author:  Joachim A.R. Sarr

  Date:    May 05th, 2020
*******************************************************************/

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
#include <new>
#include <mpi.h>

#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#include "../Non_Gray_Optimization/NG_MN_Model_3D_OPTIM.h"
#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED

#ifndef _M2_STATE_PARAMETERS_H_INCLUDED
#include "../Non_Gray_Optimization/M2_State_Parameters.h"
#endif // _M2_STATE_PARAMETERS_H_INCLUDED

#ifndef _M2_CLOSURE_FIT_EIGENSTRUCTURE_H_INCLUDED
#include "./M2_Closure_Fit_Eigenstructure.h"
#endif // _M2_CLOSURE_FIT_EIGENSTRUCTURE_H_INCLUDED

#include "../Non_Gray_Optimization/M2_Model_3D_Utilities.h"

#include "N3_Non_Gray_M2_3D_RT_Uniform.h"

using namespace nlopt;
using namespace std;

#define DISCRETE_MOBIUS_SCALE                   10
#define LEAST_SQUARES                           11

#define EVALUATION_POINTS                       1000
#define INTERPOLATION_POINTS                    1001

#define LAST_ITERATION                          10000

#define MONOMIAL_BASIS                          100000
#define CHEBYSHEV_BASIS                         100001

class N3_Non_Gray_M2_3D_RT_Cheby;

struct N3_Non_Gray_M2_3D_RT_Data_Pointer {                                                                                                   
    N3_Non_Gray_M2_3D_RT_Cheby *M2_3D_Data_Uniform;
    N3_Non_Gray_M2_3D_RT_Cheby *M2_3D_Data_BE;
    N3_Non_Gray_M2_3D_RT_Cheby *M2_3D_Data_HL;
    N3_Non_Gray_M2_3D_RT_Cheby *M2_3D_Data_LL;
};

extern int N3_ijk_Entry_Index[3];

long double myfunc_Least_Squares_L_I0_star_N3_ijk(unsigned n, const long double *x, long double *grad, void *my_func_data);
void add_constraints_Least_Squares_L_N3_ijk(nlopt_opt &opt, my_constraint_data *data, N3_Non_Gray_M2_3D_RT_Data_Pointer *data_realiz);
long double myconstraint_Least_Squares_L_N3_ijk_Hyperbolicity_Lambda_s(unsigned n, const long double *x, long double *grad, void *data);

// One can observe that, in the frame where N1_3 = 0, we also have: N2_13 = N2_23 = 0
// and: N3_113 = N3_223 = N3_333 = N3_123 = 0
// In fact, the basis for polynomials up to second order can be written as follows:
// M = {1, S1, S2, S3, S1^{2}, S1 S2, S1 S3, S2^{2}, S2 S3}
// where: S = [S1, S2, S3] = [sqrt(1 - m^{2}) cos(phi), sqrt(1 - m^{2}) sin(phi), mu]
// and any angular moment can be computed as follows: 
// I^{n} = \int_{0}^{2 PI} \int_{-1}^{1} M(n) I dmu dphi
// In the case where N1_3 = 0, by definition of angular moments, the distribution must
// be an even function of mu.

// Based on these observations, it follows that, for any given set of angular moments up to
// second-order, there exists a frame where the covariance matrix is diagonal and such that,
// N1_3 = 0, i.e., one of the eigenvectors is perpendicular to the flux vector and aligned
// with the z-axis. It follows that such a frame can be obtained by first rotating about either
// the x- or the y-axis, then followed by a rotation around the z-axis, in the intermediate frame.

class N3_Non_Gray_M2_3D_RT_Cheby {
public:
    char path_out[256], prefix[256], extension[256];
    fstream in_out;
    
    MPI_proc_params MPI_proc_parameters;
    
     // Create the rec_N3 datatype for MPI
    static MPI_Datatype rec_N3_type;
    
    static int Problem_Type;
    
    static int closure_type;
    
    static int flag_basis_type;
    
    mutable Record_Moments_Tests Rec_Moms_Test;
    
    int id_proc = -10;
    int num_proc = -10;
    int num_proc_used = -10;
    
    int num_proc_E = -10;
    int num_proc_f = -10;
    int num_proc_Phi = -10;
    int num_proc_Theta = -10;
    int num_proc_Triangle_gam1_gam2 = -10;
    
    int iteration = -10;
    
    int Node_Distribution_E = -10;
    int Node_Distribution_f = -10;
    int Node_Distribution_Theta_Phi = -10;
    int Node_Distribution_gam1 = -10;
    
    int index_f = -10;
    int index_Phi = -10;
    int index_Theta = -10;
    int index_gam1 = -10;
    int index_gam2 = -10;
    int index_triangle = -10;
    
    int N_pts_Mob_Scale = -10;
    int N_Points_E = -10;
    int N_Points_f = -10;
    int N_Points_Phi = -10;
    int N_Points_Theta = -10;
    int N_Points_gam1 = -10;
    int N_Points_Triangle_gam1_gam2 = -10;
    
    int Order_SH = -10;
    int N_Coeffs_SH = -10;
    
    int index_Ncoeffs = 0;
    int index_rec_N3 = 0;
    
    record_Ncoeffs rec_Ncoeffs;
    record_N3 *rec_N3 = NULL;
    
    int Max_Ent_Data_Type;
    
    long double *E_NON_GRAY = NULL;
    long double *N1_1_NON_GRAY = NULL;
    long double *N1_2_NON_GRAY = NULL;
    long double *N1_3_NON_GRAY = NULL;
    
    long double *gam1_NON_GRAY = NULL;
    long double *gam2_NON_GRAY = NULL;
    
    long double *x_SH = NULL;
    long double *y_SH = NULL;
    long double *z_SH = NULL;
    
    long double *N3_111_NON_GRAY = NULL;
    long double *N3_112_NON_GRAY = NULL;
    long double *N3_122_NON_GRAY = NULL;
    long double *N3_123_NON_GRAY = NULL;
    long double *N3_222_NON_GRAY = NULL;
    
    // Derivatives for N3_111
    long double *dN3_111_NON_GRAY_dN1_1 = NULL;
    long double *dN3_111_NON_GRAY_dN1_2 = NULL;
    long double *dN3_111_NON_GRAY_dN1_3 = NULL;
    long double *dN3_111_NON_GRAY_dmu = NULL;
    long double *dN3_111_NON_GRAY_dnorm_f = NULL;
    long double *dN3_111_NON_GRAY_dgam1 = NULL;
    long double *d2_N3_111_NON_GRAY_dnorm_f_dN1_1 = NULL;
    long double *d2_N3_111_NON_GRAY_dnorm_f_dgam1 = NULL;
    long double *d2_N3_111_NON_GRAY_dgam1_dN1_1 = NULL;
    long double *d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1 = NULL;
    
    // Derivatives for N3_122
    long double *dN3_122_NON_GRAY_dN1_1 = NULL;
    long double *dN3_122_NON_GRAY_dN1_2 = NULL;
    long double *dN3_122_NON_GRAY_dN1_3 = NULL;
    long double *dN3_122_NON_GRAY_dmu = NULL;
    long double *dN3_122_NON_GRAY_dnorm_f = NULL;
    long double *dN3_122_NON_GRAY_dgam1 = NULL;
    long double *dN3_122_NON_GRAY_dgam2 = NULL;
    long double *d2_N3_122_NON_GRAY_dgam1_dgam2 = NULL;
    long double *d2_N3_122_NON_GRAY_dnorm_f_dN1_1 = NULL;
    long double *d2_N3_122_NON_GRAY_dnorm_f_dgam1 = NULL;
    long double *d2_N3_122_NON_GRAY_dnorm_f_dgam2 = NULL;
    long double *d2_N3_122_NON_GRAY_dgam1_dN1_1 = NULL;
    long double *d2_N3_122_NON_GRAY_dgam2_dN1_1 = NULL;
    long double *d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1 = NULL;
    long double *d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1 = NULL;
    long double *d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2 = NULL;
    long double *d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1 = NULL;
    long double *d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1 = NULL;
    
    // Derivatives for N3_123
    long double *dN3_123_NON_GRAY_dN1_1 = NULL;
    long double *dN3_123_NON_GRAY_dN1_2 = NULL;
    long double *dN3_123_NON_GRAY_dN1_3 = NULL;
    long double *dN3_123_NON_GRAY_dmu = NULL;
    long double *dN3_123_NON_GRAY_dgam1 = NULL;
    long double *dN3_123_NON_GRAY_dgam2 = NULL;
    long double *dN3_123_NON_GRAY_dgam3 = NULL;
    
    long double *d2_N3_123_NON_GRAY_dN1_1_dN1_2 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_1_dN1_3 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_2_dN1_3 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3 = NULL;
    
    long double *d2_N3_123_NON_GRAY_dgam1_dN1_1 = NULL;
    long double *d2_N3_123_NON_GRAY_dgam2_dN1_1 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_1_dgam3 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_2_dgam1 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_2_dgam2 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_2_dgam3 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_3_dgam1 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_3_dgam2 = NULL;
    long double *d2_N3_123_NON_GRAY_dN1_3_dgam3 = NULL;
    
    long double *d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3 = NULL;
    
    long double *d2_N3_123_NON_GRAY_dgam1_dgam2 = NULL;
    long double *d2_N3_123_NON_GRAY_dgam1_dgam3 = NULL;
    long double *d2_N3_123_NON_GRAY_dgam2_dgam3 = NULL;
    long double *d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1 = NULL;
    long double *d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3 = NULL;
    long double *d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3 = NULL;
    long double *d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3 = NULL;
    
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3 = NULL;
    long double *d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3 = NULL;
    
    long double *d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2 = NULL;
    long double *d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3 = NULL;
    long double *d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3 = NULL;
    
    long double *f_N3_111_NON_GRAY = NULL;
    long double *f_N3_122_NON_GRAY = NULL;
    long double *f_N3_123_NON_GRAY = NULL;
    
    long double *Coeff_Lebedev_Quadrature_N3_111 = NULL;
    long double *Coeff_Lebedev_Quadrature_N3_122 = NULL;
    long double *Coeff_Lebedev_Quadrature_N3_123 = NULL;
    
    long double *Coefficient_Matrix_Fit_N3_111 = NULL;
    long double *Coefficient_Matrix_Fit_N3_122 = NULL;
    long double *Coefficient_Matrix_Fit_N3_123 = NULL;
    
    long double *Coefficient_Matrix_Fit_N3_111_Cheby_Basis = NULL;
    long double *Coefficient_Matrix_Fit_N3_122_Cheby_Basis = NULL;
    long double *Coefficient_Matrix_Fit_N3_123_Cheby_Basis = NULL;
    
    long double *Coefficient_Matrix_Fit_N3_111_HL = NULL;
    long double *Coefficient_Matrix_Fit_N3_111_LL = NULL;
    
    long double *Coefficient_Matrix_Fit_N3_122_HL = NULL;
    long double *Coefficient_Matrix_Fit_N3_122_LL = NULL;
    
    long double *Coefficient_Matrix_Fit_N3_123_HL = NULL;
    long double *Coefficient_Matrix_Fit_N3_123_LL = NULL;
    
    long double *Coefficient_Mobius_Scale_N3_ijk = NULL;
    
    long double *Coefficients_Vander_Matrix_SH = NULL;
    long double *Coefficients_Vander_Matrix_Least_Squares_L_I0_star = NULL;
    long double *Coefficients_Vander_Matrix_L_I0_star_N3_ijk = NULL;
    long double *Coefficients_Vander_Matrix_NG_N3_ijk = NULL;
    
    int *Array_l_SH = NULL;
    int *Array_m_SH = NULL;
    
    long double *Optim_Length_Scale_N3_ijk = NULL;
    
    long double L2_Norm_N3;
    
    long double L_inf_Norm_N3;
    
    long double L_inf_Norm_N3_111, L_inf_Norm_N3_122, L_inf_Norm_N3_123;
    
    // Object member of Closure_RT for Jacobian Computations
    mutable Closure_RT Jacobian_M2;
    
    // Constructor
    N3_Non_Gray_M2_3D_RT_Cheby(const int &Num_Coeffs_E, const int &Num_Coeffs_f, const int &Num_Points_Phi, const int &Num_Points_Theta, const int &Num_Coeffs_gam1, const int &Num_pts_Mobius_Scale, const int Order_SH_Val = 1) {
        N_Points_E = Num_Coeffs_E;
        N_Points_f = Num_Coeffs_f;
        N_Points_Phi = Num_Points_Phi;
        N_Points_Theta = Num_Points_Theta;
        N_Points_gam1 = Num_Coeffs_gam1;
        N_pts_Mob_Scale = Num_pts_Mobius_Scale;
        Order_SH = Order_SH_Val;
        
        allocate();
    }
    
    // Destructor
    ~N3_Non_Gray_M2_3D_RT_Cheby() {
        deallocate();
    }
    
    void allocate();
    
    void deallocate();
    
    void Read_Coefficient_Datas_M2_Fit();
    
    void Setup_SH_Indexes();
    
    void Setup_MPI_Processes();
    
    void Compute_MPI_Processes_Max_Min_Indexes(int &id_E_min, int &id_E_max, 
                                               int &id_f_min, int &id_f_max, 
                                               int &id_phi_min, int &id_phi_max, 
                                               int &id_theta_min, int &id_theta_max, 
                                               int &id_Triangle_gam1_gam2_min, int &id_Triangle_gam1_gam2_max);
    
    int Compute_Num_Pts_Full();
    
    int Compute_Num_Coeffs_Full_Without_SH();
    
    int Compute_Num_Coeffs_Full();
    
    int Compute_Full_Id_Triangle(const int &i_gam1, const int &i_gam2);
    
    int Compute_Id_Triangle_Single_Block(const int &i_gam1, const int &i_gam2);
    
    int Compute_N_pts_Triangle_Single_Block(const int &id_gam1_min, const int &id_gam1_max);

    int Compute_Full_Id(const int &id_Mobius, const int &i_Cheby_E, const int &i_Cheby_f, const int &i_Phi, const int &i_Theta, const int &i_Cheby_gam1_gam2);
    
    int Compute_Full_Id_Single_proc(const int &id_Mobius, const int &i_Cheby_E, const int &i_Cheby_f, const int &i_Phi, const int &i_Theta, const int &i_Cheby_gam1_gam2);
    
    
    static void Create_MPI_Data_Type_rec_N3();
    
    void OpenInputFile(char *filename);
    
    void ReadInputData();
    
    void CloseInputFile();
    
    void Print_Data(const int &index);
    
    // void SetupInterpolant_Values_HL_LL();
    
    void SetupInterpolant_Values_BE();
    
    void Lebedev_Quadrature_Interpolation();
    
    void Setup_Vandermonde_Matrix_Spherical_Harmonic();
    
    void Vandermonde_Interpolation_Spherical_Harmonic();
    
    void Setup_Vandermonde_Matrix_NG_N3_ijk();
    
    void Vandermonde_Interpolation_NG_N3_ijk();
    
    void Polynomial_Interpolation_HL_LL(const int &Maximum_Entropy_Solution_Regime);
    
    void Polynomial_Interpolation_BE(N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_HL, N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_LL);
    
    void Setup_Coefficients_HL_LL(N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_HL, N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_LL);
    
    ////////////////////////////////////////////////////////////////
    // Routines for computation of Length scale
    ////////////////////////////////////////////////////////////////
    long double Recompute_I0_Mapping(const long double &ratio_E, const long double &L_N3_ijk_Actual, const long double &L_N3_ijk_Fit);
    long double Evaluate_Length_Scale_N3(const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double Evaluate_Length_Scale_N3_ijk(const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double Evaluate_diff_Length_Scale_N3_ijk(const int &index_var_fit, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double Exponential_Mapping(const long double &I0_star, const long double &L_I0_star_N3);
    long double Inverse_Exponential_Mapping(const long double &ratio_I0_star, const long double &L_I0_star_N3);
    long double Diff_Exponential_Mapping(const long double &I0_star, const long double &L_I0_star_N3, const int &VAR_NUM);
    
    ////////////////////////////////////////////////////////////////
    // Routines for computation of third-order closing fluxes
    ////////////////////////////////////////////////////////////////
    long double Evaluate_N3_111(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double Evaluate_N3_112(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double Evaluate_N3_122(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double Evaluate_N3_123(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double Evaluate_N3_222(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double Evaluate_N3_ijk(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &VAR_NUM);
    long double Evaluate_f_N3_ijk(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &VAR_NUM);
    long double Evaluate_g_N3_ijk(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &VAR_NUM);
    long double Evaluate_g_N3_ijk_Cheby_Basis(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &VAR_NUM);
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Routines for computation of derivatives of third-order closing fluxes
    /////////////////////////////////////////////////////////////////////////////////////
    void dN3_111_2D_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_111_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    void dN3_112_2D_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_112_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    void dN3_122_2D_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_122_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    void dN3_222_2D_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_222_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    void dN3_ijk_2D_dU(record_d_N3_ijk_RT &dN3_ijk_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX);
    
    void dN3_ijk_2D_dU_Finite_Difference(record_d_N3_ijk_RT &dN3_ijk_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX);
    
    long double df_N3_ijk_interp_dU(const int &index_var_fit, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX);
    long double dg_N3_ijk_interp_dU(const int &index_var_fit, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX);
    
    long double dnorm_f2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3);
    long double dx2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3);
    long double dy2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3);
    long double dz2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3);
    long double dgam1_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double dgam2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double dB_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    
    void set_Closure_RT();
    void Compute_dN3_ijk_dN2_12_Finite_Difference(const long double &I0_star_val, const long double &N1_1, const long double &N1_2, const long double &gam1, const long double &gam2);
    void set_Closure_RT_Derivatives();
    void d_q_RT_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_ijk_RT, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX);
    void Setup_Flux_Jacobian_Matrix();
    
    long double qxxx_RT(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double qxxy_RT(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double qxyy_RT(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    long double qyyy_RT(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    
    long double qxxx();
    long double qxxy();
    long double qxyy();
    long double qyyy();
    
    void Test_Derivatives_N3_ijk_RT(record_d_N3_ijk_RT &dN3_ijk_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX);
    
    void Test_Derivatives_N3_ijk(record_d_N3_ijk_RT &dN3_ijk_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &I0_star_val, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &N2_11, const long double &N2_12, const long double &N2_22);
    
    //////////////////////////////////////////////////////////////
    // Others
    //////////////////////////////////////////////////////////////
    void Least_Square_Optimization();
    
    void Compute_L_ONE_L_TWO_Errors_N3_ijk(const int &VAR_NUM);
    
    void Compute_L_ONE_L_TWO_Errors_Derivatives_N3_ijk(const int &VAR_NUM, const int &var_index);
    
    void Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, const int &flag_Write_Output, ofstream &out_L_inf_Norm, const int &VAR_NUM);
    
    void Copy_to(N3_Non_Gray_M2_3D_RT_Cheby *New_N3_M2_3D_RT_Cheby);
    
    // void Write_Optimal_Length_Scale_Uniform(const int &VAR_NUM);
    
    void Write_Data_Matlab(const int &Max_Ent_Regime);
    
    void Polynomial_Interpolation_Gray_M2_Closure(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, ofstream &output_Opt_Coefficients);
    
    void Write_Maximum_Entropy_Data(const record_N3 *rec_N3_global);
    
    void Read_Maximum_Entropy_Data(const record_N3 *rec_N3_global);
    
    void Reorder_Maximum_Entropy_Data_MPI(record_N3 *rec_N3_global);
    
    void Reorder_Maximum_Entropy_Data_f_N3_ijk_MPI(long double *f_N3_111_NON_GRAY, long double *f_N3_122_NON_GRAY, long double *f_N3_123_NON_GRAY);
    
    //*******************************************
    // Routines Needed for Least Squares Module
    //********************************************
    void Vandermonde_Interpolation_NG_N3_ijk_Least_Squares();
    // long double Evaluate_dratio_E_dLength_Scale(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);
    void Precompute_Max_Ent_Solution_Least_Squares_L_I0_star(N3_Non_Gray_M2_3D_RT_Cheby *N3_3D_RT_Unif, const int &N_Points_f_Least_Squares, const int &N_Points_Phi_Least_Squares, const int &N_Points_Theta_Least_Squares, const int &N_Points_gam1_Least_Squares);
    void Precompute_Final_Max_Ent_Solution(const int &Maximum_Entropy_Solution_Regime = 0);
    void Compute_L_ONE_L_TWO_Errors_Least_Squares_L_N3_ijk(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif);
    
    void Least_Squares_Optimization_Coefficients_L_I0_star(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, N3_Non_Gray_M2_3D_RT_Cheby &N3_M2_3D_HL, N3_Non_Gray_M2_3D_RT_Cheby &N3_M2_3D_LL);
    
    // void Polynomial_Interpolation_Least_Squares(N3_Non_Gray_M2_3D_RT_Cheby &N3_M2_3D_RT_Unif_E_Chebyshev, N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, const int &id_proc, const int &num_proc);
    // void Compute_Optimal_L_N3_ijk_Uniform_N1_gam1_gam2(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, const int &VAR_NUM_VAL, const int &id_proc, const int &num_proc);
    void Polynomial_Interpolation_Non_Gray_M2_N3_ijk(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, ofstream &output_Opt_Coefficients);
    void Write_Coefficients_M2_Closure_Interp(ofstream &output_Opt_Coefficients);
};

inline void N3_Non_Gray_M2_3D_RT_Cheby :: Setup_Coefficients_HL_LL(N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_HL, N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_LL) {
    int N_pts_total;
    N_pts_total = N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
    for (int i = 0; i < N_Points_f; i++) {
        Coefficient_Matrix_Fit_N3_111_HL[i] = M2_3D_Data_N3_HL.Coefficient_Matrix_Fit_N3_111_HL[i];
        Coefficient_Matrix_Fit_N3_111_LL[i] = M2_3D_Data_N3_LL.Coefficient_Matrix_Fit_N3_111_LL[i];
        
        Coefficient_Matrix_Fit_N3_122_HL[i] = M2_3D_Data_N3_HL.Coefficient_Matrix_Fit_N3_122_HL[i];
        Coefficient_Matrix_Fit_N3_122_LL[i] = M2_3D_Data_N3_LL.Coefficient_Matrix_Fit_N3_122_LL[i];
        
        Coefficient_Matrix_Fit_N3_123_HL[i] = M2_3D_Data_N3_HL.Coefficient_Matrix_Fit_N3_123_HL[i];
        Coefficient_Matrix_Fit_N3_123_LL[i] = M2_3D_Data_N3_LL.Coefficient_Matrix_Fit_N3_123_LL[i];
    }   
}

// inline void N3_Non_Gray_M2_3D_RT_Cheby :: Write_Optimal_Length_Scale_Uniform(const int &VAR_NUM) {
//     int index_triangle_temp;
//     int index;
//     char path_N3_ijk[256];
//     strcpy(path_N3_ijk, getenv(PATHVAR));
//     if (VAR_NUM == N3_111_ENTRY) {
//         strcat(path_N3_ijk, "/M2_Model/Non_Gray_Model/3D_gam1_gam2_gam3/Optimal_Length_Scale_N3_111");
//     } else if (VAR_NUM == N3_122_ENTRY) {
//         strcat(path_N3_ijk, "/M2_Model/Non_Gray_Model/3D_gam1_gam2_gam3/Optimal_Length_Scale_N3_122");
//     } else if (VAR_NUM == N3_123_ENTRY) {
//         strcat(path_N3_ijk, "/M2_Model/Non_Gray_Model/3D_gam1_gam2_gam3/Optimal_Length_Scale_N3_123");
//     } else {
//         cout << "VAR_NUM not specified !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
//         exit(0);
//     }
//     strcat(path_N3_ijk, ".dat");
//     
//     fstream out_L_inf_Norm;
// //     out_L_inf_Norm.open(path_N3_ijk, ios::out|ios::binary);
//     out_L_inf_Norm.open(path_N3_ijk, ios::out);
//     long double norm_f_temp;
//     for (int i_f = 0; i_f < N_Points_f; i_f++) {
//         for (int i_Phi = 0; i_Phi < N_Points_Phi; i_Phi++) {
//             for (int i_Theta = 0; i_Theta < N_Points_Theta; i_Theta++) {
//                 index_triangle_temp = 0;
//                 for (int i_gam1 = 0; i_gam1 < N_Points_gam1; i_gam1++) {
//                     for (int i_gam2 = 0; i_gam2 < N_Points_gam1 - i_gam1; i_gam2++) {
//                         index = (i_f*N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
//                         index = index*N_Points_Triangle_gam1_gam2 + index_triangle_temp;
//                         index_triangle_temp++;
//                         
//                         norm_f_temp = pow(N1_1_NON_GRAY[index],2) + pow(N1_2_NON_GRAY[index],2) + pow(N1_3_NON_GRAY[index],2);
//                         
//                         norm_f_temp = sqrt(norm_f_temp);
//                         
//                         out_L_inf_Norm << setw(16) << norm_f_temp << setw(16) << gam1_NON_GRAY[index]<< setw(16) << gam2_NON_GRAY[index] << setw(16) << Optim_Length_Scale_N3_ijk[index] << endl;    
//                     }
//                 }
//             }
//         }
//     }
//     out_L_inf_Norm.close();
// }

inline void N3_Non_Gray_M2_3D_RT_Cheby :: Write_Data_Matlab(const int &Max_Ent_Regime) {
    int index;
    char path_out[256];
    strcpy(path_out, getenv(PATHVAR));
    
    switch (Max_Ent_Regime) {
        case GRAY:
            strcat(path_out, "/M2_Model/Gray_M2_Model/N3_3D_RT_Gray_Datas_For_Matlab.dat");
            break;
        case HYPERBOLIC_LIMIT:
            strcat(path_out, "/M2_Model/Non_Gray_M2_Model/N3_3D_RT_HL_Datas_For_Matlab.dat");
            break;
        case BOSE_EINSTEIN:
            strcat(path_out, "/M2_Model/Non_Gray_M2_Model/N3_3D_RT_BE_Datas_For_Matlab.dat");
            break;
        case LOGARITHMIC_LIMIT:
            strcat(path_out, "/M2_Model/Non_Gray_M2_Model/N3_3D_RT_LL_Datas_For_Matlab.dat");
            break;
        default:
            cout << "Invalid value of Max_Ent_Regime for writing data for visualization" << endl;
            exit(0);
            break;
    }
        
    ofstream in_out;
    in_out.open(path_out);
        
    if (!in_out) {
        cout << "N3_3D_RT_Datas_For_Matlab.dat could not be accessed!" << endl;    
    }
    
    if (in_out.good()) {
        for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
            for (int index_e = 0 ; index_e < N_Points_E; index_e++){
                for (int i_f = 0; i_f < N_Points_f; i_f++) {
                    for (int i_Phi = 0 ; i_Phi < N_Points_Phi; i_Phi++) {
                        for (int i_Theta = 0 ; i_Theta < N_Points_Theta; i_Theta++) {
                            for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_gam1_gam2++) {
                                index = id_Mobius*N_Points_E + index_e;
                                index = index*N_Points_f + i_f;
                                index = (index * N_Points_Phi+ i_Phi) * N_Points_Theta + i_Theta;
                                index = index*N_Points_Triangle_gam1_gam2 + i_gam1_gam2;
                                
                                in_out << 1.0 << setw(18) << N1_1_NON_GRAY[index] << setw(18) << N1_2_NON_GRAY[index] << setw(18) << N1_3_NON_GRAY[index] << setw(18) << gam1_NON_GRAY[index] << setw(18) << gam2_NON_GRAY[index] << setw(18) << N3_111_NON_GRAY[index] << setw(18) << N3_122_NON_GRAY[index] << setw(18) << N3_123_NON_GRAY[index] << setw(18) << f_N3_111_NON_GRAY[index] << setw(18) << f_N3_122_NON_GRAY[index] << setw(18) << f_N3_123_NON_GRAY[index] << setw(18) << dN3_111_NON_GRAY_dgam1[index] << endl;
                                
                                // cout << "N1_1 = " << N1_1_NON_GRAY[index] << "  " << "N1_2 = " << N1_2_NON_GRAY[index] << "  " << "N1_3 = " << N1_3_NON_GRAY[index] << "  " << "gam1 = " << gam1_NON_GRAY[index] << "  " << "gam2 = " << gam2_NON_GRAY[index] << "  " << "dN3_111_NON_GRAY_dgam1 = " << dN3_111_NON_GRAY_dgam1[index] << endl;
                            } // end for i_gam1_gam2
                        } // end for i_Theta
                    } // end for i_Phi
                } // end for i_f
            } // end for index_e
        } // end for id_Mobius
        in_out.close();
    }
}

inline void N3_Non_Gray_M2_3D_RT_Cheby :: Read_Coefficient_Datas_M2_Fit() {
    int N_coeffs_fit_N3_ijk, N_coeffs_fit_Length_Scale_N3_ijk;
    switch(Problem_Type) {
        case GRAY:
            if (N_Points_E != 1) {
                cout << "Issue with N_Points_E for gray radiation !!!!!!!!!!" << endl;
                exit(0);
            }
            break;
        case NON_GRAY:
            if (N_Points_E == 1) {
                cout << "Issue with N_Points_E for non-gray radiation !!!!!!!!!!" << endl;
                exit(0);
            }
            break;
    }
    
    // Set the appropriate path for the data corresponding to the 
    // gray M2 closure
    N_coeffs_fit_N3_ijk = N_Points_E*N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
    
    Coefficient_Matrix_Fit_N3_111 = new long double [N_coeffs_fit_N3_ijk];
    Coefficient_Matrix_Fit_N3_122 = new long double [N_coeffs_fit_N3_ijk];
    Coefficient_Matrix_Fit_N3_123 = new long double [N_coeffs_fit_N3_ijk];
        
    char M2_Gray_N3_ijk_BE_Coefficients_File[256];
    ifstream in_M2_Gray_N3_ijk_BE;
    
    char M2_Non_Gray_N3_ijk_BE_Coefficients_File[256];
    ifstream in_M2_Non_Gray_N3_ijk_BE;
    
    switch (Problem_Type) {
        case GRAY:
            strcpy(M2_Gray_N3_ijk_BE_Coefficients_File,getenv(PATHVAR));
            strcat(M2_Gray_N3_ijk_BE_Coefficients_File,"/M2_Model/Gray_M2_Model/Coefficients_Fit_N3_ijk_Gray_BE.dat");
            
            in_M2_Gray_N3_ijk_BE.open(M2_Gray_N3_ijk_BE_Coefficients_File);
            
            if (!in_M2_Gray_N3_ijk_BE) {
                cout << "Coefficients_Fit_N3_ijk_Non_Gray_BE.dat could not be accessed!" << endl;
                exit(0);
            }
            
            if (in_M2_Gray_N3_ijk_BE.good()) {
                for (int i = 0; i < N_coeffs_fit_N3_ijk ; i++) {
                    in_M2_Gray_N3_ijk_BE >> setprecision(12) >> Coefficient_Matrix_Fit_N3_111[i];
                }
                
                for (int i = 0; i < N_coeffs_fit_N3_ijk ; i++) {
                    in_M2_Gray_N3_ijk_BE >> setprecision(12) >> Coefficient_Matrix_Fit_N3_122[i];
                }
                
                for (int i = 0; i < N_coeffs_fit_N3_ijk ; i++) {
                    in_M2_Gray_N3_ijk_BE >> setprecision(12) >> Coefficient_Matrix_Fit_N3_123[i];
                }
                
                for (int i = 0; i < N_coeffs_fit_N3_ijk ; i++) {
                    cout << "i = " << i << "  " << "Coefficient_Matrix_Fit_N3_111 = " << Coefficient_Matrix_Fit_N3_111[i] << "  " << "Coefficient_Matrix_Fit_N3_122 = " << Coefficient_Matrix_Fit_N3_122[i] << "  " << "Coefficient_Matrix_Fit_N3_123 = " << Coefficient_Matrix_Fit_N3_123[i] << endl;
                }
                
//                 cout << "Coefficient_Matrix_Fit_N3_111[0] = " << Coefficient_Matrix_Fit_N3_111[0] << "  " << "Coefficient_Matrix_Fit_N3_122[i] = " << Coefficient_Matrix_Fit_N3_122[N_coeffs_fit_N3_ijk - 1] << endl;
//                 exit(0);
                in_M2_Gray_N3_ijk_BE.close();
            }
            break;
        case NON_GRAY:
            N_coeffs_fit_Length_Scale_N3_ijk = N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
            Coefficient_Mobius_Scale_N3_ijk = new long double[N_coeffs_fit_Length_Scale_N3_ijk];
    
            strcpy(M2_Non_Gray_N3_ijk_BE_Coefficients_File,getenv(PATHVAR));
            strcat(M2_Non_Gray_N3_ijk_BE_Coefficients_File,"/M2_Model/Non_Gray_M2_Model/Coefficients_Fit_N3_ijk_Non_Gray_BE.dat");
            
            in_M2_Non_Gray_N3_ijk_BE.open(M2_Non_Gray_N3_ijk_BE_Coefficients_File);
            
            if (!in_M2_Non_Gray_N3_ijk_BE) {
                cout << "Coefficients_Fit_N3_ijk_Non_Gray_BE.dat could not be accessed!" << endl;
                exit(0);
            }
            
            if (in_M2_Non_Gray_N3_ijk_BE.good()) {
                in_M2_Non_Gray_N3_ijk_BE >> Coefficient_Mobius_Scale_N3_ijk[0];
                // for (int i = 0; i < N_coeffs_fit_Length_Scale_N3_ijk ; i++) {
                //    in_M2_Non_Gray_N3_ijk_BE >> Coefficient_Mobius_Scale_N3_ijk[i];
                //    if (fabs(Coefficient_Mobius_Scale_N3_ijk[i]) < 1.0e-6) {
                //        Coefficient_Mobius_Scale_N3_ijk[i] = 0.0;
                //    }
                // }
                
                for (int i = 0; i < N_coeffs_fit_N3_ijk ; i++) {
                    in_M2_Non_Gray_N3_ijk_BE >> Coefficient_Matrix_Fit_N3_111[i];
                }
                
                for (int i = 0; i < N_coeffs_fit_N3_ijk ; i++) {
                    in_M2_Non_Gray_N3_ijk_BE >> Coefficient_Matrix_Fit_N3_122[i];
                }
                
                for (int i = 0; i < N_coeffs_fit_N3_ijk ; i++) {
                    in_M2_Non_Gray_N3_ijk_BE >> Coefficient_Matrix_Fit_N3_123[i];
                }
                
//                 cout << "Coefficient_Matrix_Fit_N3_111[0] = " << Coefficient_Matrix_Fit_N3_111[0] << "  " << "Coefficient_Matrix_Fit_N3_122[i] = " << Coefficient_Matrix_Fit_N3_122[N_coeffs_fit_N3_ijk - 1] << endl;
//                 exit(0);
                in_M2_Non_Gray_N3_ijk_BE.close();
            }
            break;
    }
}

#endif // _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED
