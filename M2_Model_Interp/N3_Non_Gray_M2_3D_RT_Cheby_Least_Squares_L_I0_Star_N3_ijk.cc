/*******************************************************************
  File: N3_Non_Gray_M2_3D_RT_Cheby.cc

  Description:  ...  

  Author:  Joachim A.R. Sarr

  Date:    November 15th, 2020
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

#ifndef _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED
#include "N3_Non_Gray_M2_3D_RT_Cheby.h"
#endif // _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED

// ******************************************************************************************
// This routine sets up and call the optimization algorithm for computing the optimizal length
// scale involved in the exponential mapping of the radiative energy density
// ******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Least_Squares_Optimization_Coefficients_L_I0_star(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, N3_Non_Gray_M2_3D_RT_Cheby &N3_M2_3D_HL, N3_Non_Gray_M2_3D_RT_Cheby &N3_M2_3D_LL) {
    int N_Coeffs_SH_L_N3_ijk = 1;
    int NVARS = N_Points_f*N_Coeffs_SH_L_N3_ijk*N_Points_Triangle_gam1_gam2;
    long double x[NVARS];
    long double tol_x[NVARS];
    N3_Non_Gray_M2_3D_RT_Data_Pointer N3_M2_3D_RT_Data_Pointer;
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_Cheby(N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    
    long double minf = 1.0e6, min_grad = 1.0e6; /* the minimum objective value, upon return */
    nlopt_opt opt = NULL;
    
    for (int i = 0; i < NVARS; i++) {
        x[i] = 0.0;
        tol_x[i] = 1.0e-12;
    }
    
    x[0] = 2.0;
    
    opt = nlopt_create(NLOPT_LN_COBYLA, NVARS); /* algorithm and dimensionality */
    
    Copy_to(&N3_M2_3D_RT_Cheby);
    
    N3_M2_3D_RT_Data_Pointer.M2_3D_Data_BE = &N3_M2_3D_RT_Cheby;
    N3_M2_3D_RT_Data_Pointer.M2_3D_Data_HL = &N3_M2_3D_HL;
    N3_M2_3D_RT_Data_Pointer.M2_3D_Data_LL = &N3_M2_3D_LL;
    N3_M2_3D_RT_Data_Pointer.M2_3D_Data_Uniform = &N3_3D_RT_Unif;
    
    nlopt_set_maxeval(opt, 100);
    
    nlopt_set_min_objective(opt, myfunc_Least_Squares_L_I0_star_N3_ijk, &N3_M2_3D_RT_Data_Pointer);
    
    nlopt_set_xtol_abs(opt, tol_x);
    nlopt_set_ftol_abs(opt, 1.0e-12);
    
    // my_constraint_data data[1];
    // add_constraints_Least_Squares_L_N3_ijk(opt, data, &N3_M2_3D_RT_Data_Pointer);
                                
    if (nlopt_optimize(opt, x, &minf, &min_grad) < 0) {
        printf("....................nlopt failed!.....................\n");
        exit(0);
    } else {
        if (id_proc == 0) {
            printf("***********found minimum at f(");
            for (int index_vars = 0; index_vars < NVARS; index_vars++) {
                if (index_vars < NVARS - 1) {
                    printf("%Lf,", x[index_vars]);
                } else {
                    printf("%Lf", x[index_vars]);
                }
            }
            printf(") = %0.10Le***********\n ", minf);
        }
    }
    
    // Update coefficients for the approximation of optimal length scale L_N3
    for (int i = 0; i < NVARS; i++) {
        Coefficient_Mobius_Scale_N3_ijk[i] = N3_M2_3D_RT_Cheby.Coefficient_Mobius_Scale_N3_ijk[i];
    }
}
 
// ******************************************************************************************
// This routine computes the objective function for the optimization algorithm used for the
// computation of the optimal Length scale for the exponential mapping of the radiative
// energy density
// ******************************************************************************************
long double myfunc_Least_Squares_L_I0_star_N3_ijk(unsigned n, const long double *x, long double *grad, void *my_func_data) {
    ofstream out_L_inf_Norm_temp;
    N3_Non_Gray_M2_3D_RT_Data_Pointer *N3_M2_3D_RT_Data_Pointer = (N3_Non_Gray_M2_3D_RT_Data_Pointer *) my_func_data;
    
    N3_Non_Gray_M2_3D_RT_Cheby *N3_M2_3D_RT_Cheby = N3_M2_3D_RT_Data_Pointer->M2_3D_Data_BE;
    N3_Non_Gray_M2_3D_RT_Cheby *N3_M2_3D_HL = N3_M2_3D_RT_Data_Pointer->M2_3D_Data_HL;
    N3_Non_Gray_M2_3D_RT_Cheby *N3_M2_3D_LL = N3_M2_3D_RT_Data_Pointer->M2_3D_Data_LL;
    N3_Non_Gray_M2_3D_RT_Cheby *N3_3D_RT_Unif = N3_M2_3D_RT_Data_Pointer->M2_3D_Data_Uniform;
    
    // Setup the current iterate for the length scale of the exponential mapping of the radiative
    // energy density
    for (int i = 0; i < n; i++) {
        N3_M2_3D_RT_Cheby->Coefficient_Mobius_Scale_N3_ijk[i] = x[i];
        if (N3_M2_3D_RT_Cheby->id_proc == PRIMARY_ID) {
            cout << "i = " << i << "  " << "Coefficient_Matrix_Fit_Mob_Scale = " << x[i] << endl;
        }
    }
    
    long double objective_function;
    int NF_obj = 1;
    
    //N3_M2_3D_RT_Cheby->f_L_N3_ijk = N3_M2_3D_RT_Cheby->OPTIM_NON_GRAY_M2_Least_Squares_f_L_N3_ijk(*N3_M2_3D_RT_Cheby, *N3_M2_3D_HL, *N3_M2_3D_LL);
    
    // Now perform polynomial interpolation
    N3_M2_3D_RT_Cheby->Polynomial_Interpolation_BE(*N3_M2_3D_HL, *N3_M2_3D_LL);
    
    // Now assess the error of our interpolative-based approximation of the Eddington factor based on the curremt
    // iterate for the coefficients of the interpolative-based approximation of the optimal length scale for the 
    // exponential mapping of the radiative energy density, as well as the associated optimal value for f_L_N3
    
    N3_M2_3D_RT_Cheby->Compute_L_ONE_L_TWO_Errors_Least_Squares_L_N3_ijk(*N3_3D_RT_Unif);
    
    objective_function = N3_3D_RT_Unif->L2_Norm_N3;
    objective_function = objective_function*objective_function;
    
    return objective_function;   
}

void add_constraints_Least_Squares_L_N3_ijk(nlopt_opt &opt, my_constraint_data *data, N3_Non_Gray_M2_3D_RT_Data_Pointer *data_realiz) {
    // nlopt_add_inequality_constraint(opt, myconstraint_Least_Squares_L_N3_ijk_Hyperbolicity_Lambda_s, data_realiz, 0.0);
}
 
// ******************************************************************************************
// This routine enforces constraints for the eigenvalues of the flux Jacobian of the system
// of closed moment equations based on the M2 closure, for hyperbolicity
// ******************************************************************************************
long double myconstraint_Least_Squares_L_N3_ijk_Hyperbolicity_Lambda_s(unsigned n, const long double *x, long double *grad, void *data) {
    long double constr_val;
    long double lam_Imag_min;
    long double ratio_E, ratio_E_orig, I0_star_val;
    long double norm_f, theta, phi;
    long double N1_1, N1_2, N1_3;
    long double gam1, gam2;
    long double max_lambda;
    int num_pts_temp;
    long double lam_Real[6], lam_Imag[6];
    
    N3_Non_Gray_M2_3D_RT_Data_Pointer *N3_M2_3D_RT_Data_Pointer = (N3_Non_Gray_M2_3D_RT_Data_Pointer *) data;
    
    N3_Non_Gray_M2_3D_RT_Cheby *N3_M2_3D_RT_Cheby = N3_M2_3D_RT_Data_Pointer->M2_3D_Data_BE;
    N3_Non_Gray_M2_3D_RT_Cheby *N3_M2_3D_HL = N3_M2_3D_RT_Data_Pointer->M2_3D_Data_HL;
    N3_Non_Gray_M2_3D_RT_Cheby *N3_M2_3D_LL = N3_M2_3D_RT_Data_Pointer->M2_3D_Data_LL;
    N3_Non_Gray_M2_3D_RT_Cheby *N3_3D_RT_Unif = N3_M2_3D_RT_Data_Pointer->M2_3D_Data_Uniform;
    
    max_lambda = 0.0;
    num_pts_temp = 100;
    
    
    int N_pts_Mob_Scale = 100;
    long double Mobius_Scale_Actual, L_I0_star_N3;
    long double alpha_Lag;
    alpha_Lag = 0;
    long double x_L_N3[N_pts_Mob_Scale], w_L_N3[N_pts_Mob_Scale];
    gen_laguerre_ek_compute ( N_pts_Mob_Scale, alpha_Lag, x_L_N3, w_L_N3 );
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_N3[id_Mobius];
        for (int i_E = 0; i_E < num_pts_temp; i_E++){
            ratio_E = zeros_shifted(i_E, num_pts_temp, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            for (int i_f = 0; i_f < num_pts_temp; i_f++) {
                norm_f = zeros_shifted(i_f, num_pts_temp - 1, 0.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
                for (int i_Phi = 0; i_Phi < num_pts_temp; i_Phi++) {
                    phi = zeros_shifted(i_Phi, num_pts_temp - 1, 0.0, 2.0*PI, UNIFORM_DISTRIBUTION);
                    for (int i_Theta = 0; i_Theta < num_pts_temp; i_Theta++) {
                        theta = zeros_shifted(i_Theta, num_pts_temp - 1, 0.0, PI, UNIFORM_DISTRIBUTION);
                        theta = PI/2.0;
                        for (int i_gam1 = 0; i_gam1 < num_pts_temp; i_gam1++) {
                            for (int i_gam2 = 0; i_gam2 < num_pts_temp - i_gam1; i_gam2++) {
                                Triangle_to_Square_Mapping_Nodes(gam1, gam2, i_gam1, i_gam2, num_pts_temp, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
                                
                                N1_1 = norm_f*sin(theta)*cos(phi);
                                N1_2 = norm_f*sin(theta)*sin(phi);
                                N1_3 = norm_f*cos(theta);
                                
                                if (fabs(N1_3) > 1.0e-12) {
                                    cout << "Problem in 2D !!!!!!!!!!!!!!!!!!!!!!" << endl;
                                    cout << "N1_1 = " << N1_1 << endl;
                                    exit(0);
                                }
                                
                                L_I0_star_N3 = N3_M2_3D_RT_Cheby->Evaluate_Length_Scale_N3(N1_1, N1_2, N1_3, gam1, gam2);
                                ratio_E = N3_M2_3D_RT_Cheby->Recompute_I0_Mapping(ratio_E_orig, Mobius_Scale_Actual, L_I0_star_N3);
                                
                                cout << "Compute I0_star_val with current lenght scale distribution" << endl;
                                I0_star_val = 0.0;
                                
                                N3_M2_3D_RT_Cheby->Jacobian_M2.Compute_W_array(I0_star_val, N1_1, N1_2, gam1, gam2);
                                N3_M2_3D_RT_Cheby->Setup_Flux_Jacobian_Matrix();
                                
                                Eigenvalues_overwrite( N3_M2_3D_RT_Cheby->Jacobian_M2.dFdU, lam_Real, lam_Imag);
                                
                                for (int i = 0; i < 6; i++) {
                                    lam_Imag_min = min(lam_Imag_min, lam_Imag[i]);
                                }
                                constr_val = min(constr_val, lam_Imag_min);
                                // max_lambda = max(max_lambda, fabs(Lambdas));
                            }
                        }
                    }
                }
            }
        }
    }
    
    constr_val = 1.0e6*constr_val;
    
    if (N3_M2_3D_RT_Cheby->id_proc == 0) {
        cout << "f_L_N3 = " << x[0] << "  " << "diff_lambdas_max = " << constr_val << endl;
    }
    
    return constr_val;
}

// ******************************************************************************************
// This routine computes the L1 and L2 errors of our proposed interpolative-based approximation
// of the closing fluxes for the non-gray M2 closure for values of the radiative energy density
// ranging from the hyperbolic limit to the logarithmic limit, while the values of N1_1, N1_2, 
// N1_3, \gamma1, and \gamma2 are fixed.
// ******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Compute_L_ONE_L_TWO_Errors_Least_Squares_L_N3_ijk(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif) {
    long double max_err_N3_N1_1, max_err_N3_N1_2, max_err_N3_N1_3, max_err_N3_gam1, max_err_N3_gam2;
    long double L2_Norm_N3_111, L_inf_Norm_N3_111;
    long double L2_Norm_N3_122, L_inf_Norm_N3_122;
    long double L2_Norm_N3_123, L_inf_Norm_N3_123;
    long double error_N3_111, error_N3_122, error_N3_123;
    long double ratio_E, N1_1, N1_2, N1_3, gam1, gam2;
    long double N3_111_Fit, N3_111_Numerical;
    long double N3_122_Fit, N3_122_Numerical;
    long double N3_123_Fit, N3_123_Numerical;
    long double Mobius_Scale_Actual, Mobius_Scale_Fit;
    int index;
    int id_E_min, id_E_max;
    int id_f_min, id_f_max;
    int id_phi_min, id_phi_max;
    int id_theta_min, id_theta_max;
    int id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max;
    long double *L2_Norm_N3_array = NULL, *L_inf_Norm_N3_111_array = NULL, *L_inf_Norm_N3_122_array = NULL, *L_inf_Norm_N3_123_array = NULL;
    long double w_N3_111, w_N3_122, w_N3_123;
    long double L2_Norm_N3;
    long double weight_total;
    long double weight_L_N3_val, weight_r_I0_val, weight_N1_val;
    w_N3_111 = 0.4*10.0;
    w_N3_122 = 0.3*10.0;
    w_N3_123 = 0.3*10.0;
    
    if (id_proc == PRIMARY_ID) {
        L2_Norm_N3_array = new long double[num_proc];
        L_inf_Norm_N3_111_array = new long double[num_proc];
        L_inf_Norm_N3_122_array = new long double[num_proc];
        L_inf_Norm_N3_123_array = new long double[num_proc];
    }
    
    error_N3_111 = 0.0;
    error_N3_122 = 0.0;
    error_N3_123 = 0.0;
    
    L_inf_Norm_N3_111 = 0.0;
    L_inf_Norm_N3_122 = 0.0;
    L_inf_Norm_N3_123 = 0.0;
    
    long double x_L_N3[N3_3D_RT_Unif.N_pts_Mob_Scale], weight_L_N3[N3_3D_RT_Unif.N_pts_Mob_Scale];
    long double x_r_I0[N3_3D_RT_Unif.N_Points_E], weight_r_I0[N3_3D_RT_Unif.N_Points_E];
    long double x_N1[2*(N3_3D_RT_Unif.N_Points_f - 1) + 1], weight_N1[2*(N3_3D_RT_Unif.N_Points_f - 1) + 1];
    
    chebyshev_quadrature ( weight_r_I0, x_r_I0, N3_3D_RT_Unif.N_Points_E, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
    chebyshev_quadrature ( weight_N1, x_N1, 2*(N3_3D_RT_Unif.N_Points_f - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
    gen_laguerre_ek_compute ( N3_3D_RT_Unif.N_pts_Mob_Scale, 0.0, x_L_N3, weight_L_N3 );
    
    N3_3D_RT_Unif.Compute_MPI_Processes_Max_Min_Indexes(id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, id_theta_min, id_theta_max, id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max);
    
    for (int id_Mobius = 0; id_Mobius < N3_3D_RT_Unif.N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_N3[id_Mobius];
        for (int i_E = id_E_min; i_E < id_E_max; i_E++) {
            for (int i_f = id_f_min; i_f < id_f_max; i_f++) {
                for (int i_Phi = id_phi_min; i_Phi < id_phi_max; i_Phi++) {
                    for (int i_Theta = id_theta_min; i_Theta < id_theta_max; i_Theta++) {
                        for (int i_gam1_gam2 = id_Triangle_gam1_gam2_min; i_gam1_gam2 < id_Triangle_gam1_gam2_max; i_gam1_gam2++) {
                            index = N3_3D_RT_Unif.Compute_Full_Id(id_Mobius, i_E, i_f, i_Phi, i_Theta, i_gam1_gam2);
                                
                            ratio_E = N3_3D_RT_Unif.E_NON_GRAY[index];
                            N1_1 = N3_3D_RT_Unif.N1_1_NON_GRAY[index];
                            N1_2 = N3_3D_RT_Unif.N1_2_NON_GRAY[index];
                            N1_3 = N3_3D_RT_Unif.N1_3_NON_GRAY[index];
                            gam1 = N3_3D_RT_Unif.gam1_NON_GRAY[index];
                            gam2 = N3_3D_RT_Unif.gam2_NON_GRAY[index];
                            
                            Mobius_Scale_Fit = Evaluate_Length_Scale_N3_ijk(N1_1, N1_2, N1_3, gam1, gam2);
                            
                            if (fabs(ratio_E) == 1.0) {
                                ratio_E = ratio_E;
                            } else {
                                //                     cout << "ratio_E Before = " << ratio_E << endl;
                                ratio_E = Inverse_Mobius_Transformation(ratio_E, Mobius_Scale_Actual);
                                //                     cout << "I0 Before = " << ratio_E << endl;
                                ratio_E = Mobius_Transformation(ratio_E, Mobius_Scale_Fit);
                                //                     cout << "ratio_E after = " << ratio_E << endl;
                            }
                            
                            // N3_111_Numerical = N3_3D_RT_Unif.N3_111_NON_GRAY[index];
                            // N3_111_Fit = Evaluate_N3_111(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            N3_111_Numerical = N3_3D_RT_Unif.f_N3_111_NON_GRAY[index];
                            N3_111_Fit = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
                            
                            // N3_122_Numerical = N3_3D_RT_Unif.N3_122_NON_GRAY[index];
                            // N3_122_Fit = Evaluate_N3_122(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            N3_122_Numerical = N3_3D_RT_Unif.f_N3_122_NON_GRAY[index];
                            N3_122_Fit = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
                            
                            N3_123_Numerical = 0.0; // N3_3D_RT_Unif.N3_123_NON_GRAY[index];
                            N3_123_Fit = 0.0; // Evaluate_N3_123(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            
                            weight_L_N3_val = exp(x_L_N3[id_Mobius]) * weight_L_N3[id_Mobius];
                            if (isinf(weight_L_N3_val) || isnan(weight_L_N3_val)) {
                                weight_L_N3_val = 0.0;
                                // cout << "weight_L_N3_val = " << weight_L_N3_val << endl;
                                // cout << "id_Mobius = " << id_Mobius << "  " << "x_L_N3 = " << x_L_N3[id_Mobius] << "  " << "weight_L_N3 = " << weight_L_N3[id_Mobius] << endl;
                            }
                            
                            weight_r_I0_val = weight_r_I0[i_E];
                            weight_N1_val = weight_N1[i_f + (N3_3D_RT_Unif.N_Points_f - 1)];
                            weight_total = weight_L_N3_val * weight_r_I0_val * weight_N1_val;
                            
                            error_N3_111 += /*weight_total **/ pow(N3_111_Fit - N3_111_Numerical, 2);
                            error_N3_122 += /*weight_total **/ pow(N3_122_Fit - N3_122_Numerical, 2);
                            error_N3_123 += /*weight_total **/ pow(N3_123_Fit - N3_123_Numerical, 2);
                            
//                             if (id_proc == PRIMARY_ID) {
//                                 cout << "  " << "ratio_E = " << ratio_E << "  " << "N1_1 = " << N1_1 << "  " << "N1_2 = " << N1_2 << "  " << "N1_3 = " << N1_3 << "  " << "gam1 = " << gam1 << "  " << "gam2 = " << gam2 << "  " << "N3_111_Numerical = " << N3_111_Numerical << "  " << "N3_122_Numerical = " << N3_122_Numerical << "  " << "N3_111_Fit = " << N3_111_Fit << "  " << "N3_122_Fit = " << N3_122_Fit << endl;
//                             }
                            
                            L_inf_Norm_N3_111 = max(L_inf_Norm_N3_111, fabs(N3_111_Fit - N3_111_Numerical));
                            L_inf_Norm_N3_122 = max(L_inf_Norm_N3_122, fabs(N3_122_Fit - N3_122_Numerical));
                            L_inf_Norm_N3_123 = max(L_inf_Norm_N3_123, fabs(N3_123_Fit - N3_123_Numerical));
//                                 
//                             if (L_inf_Norm_N3 == fabs(N3_Fit - N3_Numerical)) {
//                                 max_err_N3_N1_1 = N1_1;
//                                 max_err_N3_N1_2 = N1_2;
//                                 max_err_N3_N1_3 = N1_3;
//                                 max_err_N3_gam1 = gam1;
//                                 max_err_N3_gam2 = gam2;
//                             }
                        }
                    }
                }
            }
        }
    }
    
    L2_Norm_N3 = w_N3_111*error_N3_111 + w_N3_122*error_N3_122 + w_N3_123*error_N3_123;
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Gather all maximum entropy solutions to primary processor 
    MPI_Gather(&L2_Norm_N3, 1, MPI_LONG_DOUBLE, L2_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Gather(&L_inf_Norm_N3_111, 1, MPI_LONG_DOUBLE, L_inf_Norm_N3_111_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Gather(&L_inf_Norm_N3_122, 1, MPI_LONG_DOUBLE, L_inf_Norm_N3_122_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Gather(&L_inf_Norm_N3_123, 1, MPI_LONG_DOUBLE, L_inf_Norm_N3_123_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Compute L2 and L_inf norm among processors
    L2_Norm_N3 = 0.0;
    L_inf_Norm_N3_111 = 0.0;
    L_inf_Norm_N3_122 = 0.0;
    L_inf_Norm_N3_123 = 0.0;
    if (id_proc == PRIMARY_ID) {
        for (int id_prc = 0; id_prc < N3_3D_RT_Unif.num_proc_used; id_prc++) {
            L2_Norm_N3 += L2_Norm_N3_array[id_prc];
            L_inf_Norm_N3_111 = max(L_inf_Norm_N3_111, fabs(L_inf_Norm_N3_111_array[id_prc]));
            L_inf_Norm_N3_122 = max(L_inf_Norm_N3_122, fabs(L_inf_Norm_N3_122_array[id_prc]));
            L_inf_Norm_N3_123 = max(L_inf_Norm_N3_123, fabs(L_inf_Norm_N3_123_array[id_prc]));
        }
        
        L2_Norm_N3 /= N3_3D_RT_Unif.N_pts_Mob_Scale*N3_3D_RT_Unif.N_Points_E*N3_3D_RT_Unif.N_Points_f*N3_3D_RT_Unif.N_Points_Phi*N3_3D_RT_Unif.N_Points_Theta*N3_3D_RT_Unif.N_Points_Triangle_gam1_gam2;
        L2_Norm_N3 = sqrt(L2_Norm_N3);
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Broadcast errors on all the other processors
    MPI_Bcast(&L_inf_Norm_N3_111, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(&L_inf_Norm_N3_122, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(&L_inf_Norm_N3_123, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(&L2_Norm_N3, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    N3_3D_RT_Unif.L_inf_Norm_N3_111 = L_inf_Norm_N3_111;
    N3_3D_RT_Unif.L_inf_Norm_N3_122 = L_inf_Norm_N3_122;
    N3_3D_RT_Unif.L_inf_Norm_N3_123 = L_inf_Norm_N3_123;
    N3_3D_RT_Unif.L2_Norm_N3 = L2_Norm_N3;
        
    if (id_proc == PRIMARY_ID) {
        cout << "Convergence Stats ........." << "L2_Norm_N3 = " << L2_Norm_N3 << "     "  << "L_inf_Norm_N3_111 = " << L_inf_Norm_N3_111 << "     "  << "L_inf_Norm_N3_122 = " << L_inf_Norm_N3_122 << "     "  << "L_inf_Norm_N3_123 = " << L_inf_Norm_N3_123 << "    " << "max_err_N3_N1_1 = " << max_err_N3_N1_1 << "    " << "max_err_N3_N1_2 = " << max_err_N3_N1_2 << "    " << "max_err_N3_N1_3 = " << max_err_N3_N1_3 << "    " << "max_err_N3_gam1 = " << max_err_N3_gam1 << "    " << "max_err_N3_gam2 = " << max_err_N3_gam2 << endl;
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        delete[] L2_Norm_N3_array; L2_Norm_N3_array = NULL;
        delete[] L_inf_Norm_N3_111_array; L_inf_Norm_N3_111_array = NULL;
        delete[] L_inf_Norm_N3_122_array; L_inf_Norm_N3_122_array = NULL;
        delete[] L_inf_Norm_N3_123_array; L_inf_Norm_N3_123_array = NULL;
    }
}
