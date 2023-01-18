#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#include "NG_MN_Model_3D_OPTIM.h"
#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED

// ******************************************************************************
// This routine computes entries of the Hessian matrix of second-derivatives
// of the objective function (or first-derivatives of the angular moments)
// ===> H^{(1)}_{ij} = {\partial I^{(i)}} {\partial \lam_j}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int H_ONE_Matrix_Entries(long double *H_ONE, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     int index_p, index_a;
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     index_p = M2_State->index_p;
     index_a = M2_State->index_a;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double poly_Sk_cur, poly_ak;
     poly = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     poly_Sk = new long double[M2_State->NVARS];
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
     }
     
     poly_Sk_cur = 0.0;
     poly_ak = 0.0;
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk_cur += M2_State->Sk_cur[index_p*M2_State->NVARS+i]*poly[i];
         poly_ak += M2_State->ak[index_a*M2_State->NVARS+i]*poly[i];
     }
     
     long double summation;
     
     switch (M2_State->Regime) {
         case GRAY:
//              if (sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS) > 0.0) {
//                  cout << "sum_Lagrange = " << sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS) << endl; 
//                  exit(0);
//              }
             H_ONE[0] = -4.0*poly_Sk_cur*poly_ak/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 5));
             break;
         case BOSE_EINSTEIN:
             // use this to avoid overflow for the exponential function
             summation = sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS);
             H_ONE[0] = coeff_1*poly_Sk_cur*poly_ak*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))/pow((1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))), 2.0);
             break;
         case HYPERBOLIC_LIMIT:
             H_ONE[0] = coeff_1*poly_Sk_cur*poly_ak*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             H_ONE[0] = coeff_1*poly_Sk_cur*poly_ak/pow(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 2.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }

// ******************************************************************************
// This routine computes entries of the second-derivatives of the angular moments
// with respect to the Lagrange multipliers
// ===> H^{(2)}_{ijk} = {\partial^{2} I^{(i)}} {\partial \lam_j\partial \lam_k}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int H_TWO_Matrix_Entries(long double *H_TWO, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1, exp_x;
     int index_Lag1, index_Lag2, index_Moment;
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     index_Moment = M2_State->index_Moment;
     index_Lag1 = M2_State->index_Lag_i;
     index_Lag2 = M2_State->index_Lag_j;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double poly_Sk_cur, poly_Lag_1, poly_Lag_2;
     poly = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     poly_Sk = new long double[M2_State->NVARS];
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
     }
     
     poly_Sk_cur = 0.0;
     poly_Lag_1 = 0.0;
     poly_Lag_2 = 0.0;
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk_cur += M2_State->Sk_cur[index_Moment*M2_State->NVARS+i]*poly[i];
         poly_Lag_1 += M2_State->ak[index_Lag1*M2_State->NVARS+i]*poly[i];
         poly_Lag_2 += M2_State->ak[index_Lag2*M2_State->NVARS+i]*poly[i];
     }
     
     switch (M2_State->Regime) {
         case GRAY:
             H_TWO[0] = 20.0*poly_Sk_cur*poly_Lag_1*poly_Lag_2/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 6));
             break;
         case BOSE_EINSTEIN:
             exp_x = exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             H_TWO[0] = coeff_1*poly_Sk_cur*poly_Lag_1*poly_Lag_2*exp_x*(1.0 + exp_x)/pow((1.0 - exp_x), 3.0);
             break;
         case HYPERBOLIC_LIMIT:
             H_TWO[0] = coeff_1*poly_Sk_cur*poly_Lag_1*poly_Lag_2*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             H_TWO[0] = 2.0*coeff_1*poly_Sk_cur*poly_Lag_1*poly_Lag_2/pow(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 3.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
// ******************************************************************************
// This routine computes entries of the third-derivatives of the angular moments
// with respect to the Lagrange multipliers
// ===> H^{(3)}_{ijkl} = {\partial^{3} I^{(i)}} {\partial \lam_j \partial \lam_k \partial \lam_l}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int H_THREE_Matrix_Entries(long double *H_THREE, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1, exp_x;
     int index_Lag1, index_Lag2, index_Lag3, index_Moment;
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     index_Moment = M2_State->index_Moment;
     index_Lag1 = M2_State->index_Lag_i;
     index_Lag2 = M2_State->index_Lag_j;
     index_Lag3 = M2_State->index_Lag_k;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double poly_Sk_cur, poly_Lag_1, poly_Lag_2, poly_Lag_3;
     poly = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     poly_Sk = new long double[M2_State->NVARS];
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
     }
     
     poly_Sk_cur = 0.0;
     poly_Lag_1 = 0.0;
     poly_Lag_2 = 0.0;
     poly_Lag_3 = 0.0;
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk_cur += M2_State->Sk_cur[index_Moment*M2_State->NVARS+i]*poly[i];
         poly_Lag_1 += M2_State->ak[index_Lag1*M2_State->NVARS+i]*poly[i];
         poly_Lag_2 += M2_State->ak[index_Lag2*M2_State->NVARS+i]*poly[i];
         poly_Lag_3 += M2_State->ak[index_Lag3*M2_State->NVARS+i]*poly[i];
     }
     
     switch (M2_State->Regime) {
         case GRAY:
             H_THREE[0] = -120.0*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 7));
             break;
         case BOSE_EINSTEIN:
             exp_x = exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             H_THREE[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur*exp_x*(1.0 + 4.0*exp_x + exp_x*exp_x)/pow((1.0 - exp_x), 4.0);
             break;
         case HYPERBOLIC_LIMIT:
             H_THREE[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             H_THREE[0] = 6.0*coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur/pow(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
// ******************************************************************************
// This routine computes entries of the first-derivatives of the highest order
// moments with respect to the Lagrange multipliers
// ===> {\partial I^{(N)}} {\partial \lam_j}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int Higher_Order_Moments_Derivatives(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double *poly_Sk_cur, poly_ak;
     poly = new long double[M2_State->NVARS];
     poly_Sk = new long double[M2_State->NVARS];
     poly_Sk_cur = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
     }
     
     poly_ak = 0.0;
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk_cur[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk_cur[i] += M2_State->Sk_cur[i*M2_State->NVARS+j]*poly[j];
         }
     }
     poly_ak = generate_higher_order_monomials_basis(M2_State->Omega1[Id_angle],M2_State->Omega2[Id_angle],M2_State->Omega3[Id_angle], M2_State->Index_Higher_Order_Moments);
     
     switch (M2_State->Regime) {
         case GRAY:
             for (int i = 0; i < M2_State->NVARS; i++) {
                 F_Hess[i] = -4.0*poly_Sk_cur[i]*poly_ak/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 5));
             }
             break;
         case BOSE_EINSTEIN:
             for (int i = 0; i < M2_State->NVARS; i++) {
                 F_Hess[i] = coeff_1*poly_Sk_cur[i]*poly_ak*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))/pow((1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))), 2.0);
             }
             break;
         case HYPERBOLIC_LIMIT:
             for (int i = 0; i < M2_State->NVARS; i++) {
                 F_Hess[i] = coeff_1*poly_Sk_cur[i]*poly_ak*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             }
             break;
         case LOGARITHMIC_LIMIT:
             for (int i = 0; i < M2_State->NVARS; i++) {
                 F_Hess[i] = coeff_1*poly_Sk_cur[i]*poly_ak/pow(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 2.0);
            }
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk_cur;
     delete[] poly_Sk;
     return 0;
 }
 
// ******************************************************************************
// This routine computes entries of the second-derivatives of the highest order
// moments with respect to the Lagrange multipliers
// ===> {\partial^{2} I^{(N)}} {\partial \lam_j \partial \lam_k}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int Higher_Order_Moments_Second_Derivatives(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1, exp_x;
     int index_Lag1, index_Lag2;
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     index_Lag1 = M2_State->index_Lag_i;
     index_Lag2 = M2_State->index_Lag_j;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double poly_Sk_cur, poly_Lag_1, poly_Lag_2;
     poly = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     poly_Sk = new long double[M2_State->NVARS];
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
     }
     
     poly_Lag_1 = 0.0;
     poly_Lag_2 = 0.0;
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Lag_1 += M2_State->ak[index_Lag1*M2_State->NVARS+i]*poly[i];
         poly_Lag_2 += M2_State->ak[index_Lag2*M2_State->NVARS+i]*poly[i];
     }
     
     poly_Sk_cur = generate_higher_order_monomials_basis(M2_State->Omega1[Id_angle],M2_State->Omega2[Id_angle],M2_State->Omega3[Id_angle], M2_State->Index_Higher_Order_Moments);
     
     switch (M2_State->Regime) {
         case GRAY:
             F_Hess[0] = 20.0*poly_Lag_1*poly_Lag_2*poly_Sk_cur/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 6));
             break;
         case BOSE_EINSTEIN:
             exp_x = exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             F_Hess[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Sk_cur*exp_x*(1.0 + exp_x)/pow((1.0 - exp_x), 3.0);
             break;
         case HYPERBOLIC_LIMIT:
             F_Hess[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Sk_cur*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             F_Hess[0] = 2.0*coeff_1*poly_Lag_1*poly_Lag_2*poly_Sk_cur/pow(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 3.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
// ******************************************************************************
// This routine computes entries of the third-derivatives of the highest order
// moments with respect to the Lagrange multipliers
// ===> {\partial^{3} I^{(N)}} {\partial \lam_j \partial \lam_k \partial \lam_l}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int Higher_Order_Moments_Third_Derivatives(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1, exp_x;
     int index_Lag1, index_Lag2, index_Lag3;
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     index_Lag1 = M2_State->index_Lag_i;
     index_Lag2 = M2_State->index_Lag_j;
     index_Lag3 = M2_State->index_Lag_k;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double poly_Sk_cur, poly_Lag_1, poly_Lag_2, poly_Lag_3;
     poly = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     poly_Sk = new long double[M2_State->NVARS];
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
     }
     
     poly_Lag_1 = 0.0;
     poly_Lag_2 = 0.0;
     poly_Lag_3 = 0.0;
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Lag_1 += M2_State->ak[index_Lag1*M2_State->NVARS+i]*poly[i];
         poly_Lag_2 += M2_State->ak[index_Lag2*M2_State->NVARS+i]*poly[i];
         poly_Lag_3 += M2_State->ak[index_Lag3*M2_State->NVARS+i]*poly[i];
     }
     
     poly_Sk_cur = generate_higher_order_monomials_basis(M2_State->Omega1[Id_angle],M2_State->Omega2[Id_angle],M2_State->Omega3[Id_angle], M2_State->Index_Higher_Order_Moments);
     
     switch (M2_State->Regime) {
         case GRAY:
             F_Hess[0] = -120.0*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 7));
             break;
         case BOSE_EINSTEIN:
             exp_x = exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             F_Hess[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur*exp_x*(1.0 + 4.0*exp_x + exp_x*exp_x)/pow((1.0 - exp_x), 4.0);
             break;
         case HYPERBOLIC_LIMIT:
             F_Hess[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             F_Hess[0] = 6.0*coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur/pow(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
// ******************************************************************************
// This routine computes the Hessian matrix of first-derivatives of the 
// angular moments with respect to the Lagrange multipliers
// ===> H^{(1)}_{ij} = {\partial^{2} I^{(i)}} {\partial \lam_j}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 void Compute_H_ONE_Matrix(long double *H_ONE, const long double *x, const long double *Sk, M2_State_Param &M2_State) {
     long double H_ONE_ij;
     M2_State.set_x(x);
     M2_State.set_Sk(Sk);
     M2_State.set_Sk_cur(Sk);
     M2_State.set_ak(Sk);
     
     for (int i = 0; i < M2_State.NVARS; i++) {
         for (int j = 0; j < M2_State.NVARS; j++) {
             M2_State.index_p = i;
             M2_State.index_a = j;
             Lebedev_Quadrature_Func(&H_ONE_ij, 1, H_ONE_Matrix_Entries, &M2_State);
             H_ONE[i*M2_State.NVARS + j] = H_ONE_ij;
         }
     }
 }
 
// ******************************************************************************
// This routine computes the matrix of second-derivatives of the angular moments
// with respect to the Lagrange multipliers
// ===> H^{(2)}_{ijk} = {\partial^{2} I^{(i)}} {\partial \lam_j\partial \lam_k}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 void Compute_H_TWO_Matrix(long double *H_TWO, const long double *x, const long double *Sk, M2_State_Param &M2_State) {
     long double H_TWO_ijk;
     M2_State.set_x(x);
     M2_State.set_Sk(Sk);
     M2_State.set_Sk_cur(Sk);
     M2_State.set_ak(Sk);
     
     for (int i = 0; i < M2_State.NVARS; i++) {
         for (int j = 0; j < M2_State.NVARS; j++) {
             for (int k = 0; k < M2_State.NVARS; k++) {
                 M2_State.index_Moment = i;
                 M2_State.index_Lag_i = j;
                 M2_State.index_Lag_j = k;
                 Lebedev_Quadrature_Func(&H_TWO_ijk, 1, H_TWO_Matrix_Entries, &M2_State);
                 H_TWO[(i*M2_State.NVARS + j)*M2_State.NVARS + k] = H_TWO_ijk;
             }
         }
     } 
 }

// ******************************************************************************
// This routine computes the matrix of first-derivatives of the Lagrange 
// multipliers with respect to the angular moments (which is also the inverse of
// the Hessian matrix)
// ===> A^{(1)}_{ij} = {\partial^{2} \lam_i} {\partial I^{(j)}}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 void Compute_A_ONE_Matrix(long double *A_ONE, const long double *x, const long double *Sk, M2_State_Param &M2_State) {
     long double *Q, *R, *H_ONE;
     H_ONE = new long double[M2_State.NVARS*M2_State.NVARS];
     Q = new long double[M2_State.NVARS*M2_State.NVARS];
     R = new long double[M2_State.NVARS*M2_State.NVARS];
     
     // Compute Hessian Matrix H_ONE
     Compute_H_ONE_Matrix(H_ONE, x, Sk, M2_State);
     
     
//      cout << "************************** Original !!!!!!!! ***********************" << endl;
// 
//      for (int i = 0; i < M2_State.NVARS; i++) {
//          for (int j = 0; j < M2_State.NVARS; j++) {
//              cout << "i = " << i << "   " << "j = " << j << "   " << "H_ONE = " << H_ONE[i*M2_State.NVARS+j] << endl;
//         }
//     }
     
     // ******************************************************************************
     // Perform QR decomposition on H_ONE using Modified Gram Schmidt Factorization
     // ==> H_ONE = Q R
     // where R is upper triangular
     // ******************************************************************************
     Modified_Gram_Schmidt_Factorization(M2_State.NVARS, H_ONE, Q, R);
     
//      for (int i = 0; i < M2_State.NVARS; i++) {
//          for (int j = 0; j < M2_State.NVARS; j++) {
//              cout << "i = " << i << "   " << "j = " << j << "   " << "Q = " << Q[i*M2_State.NVARS+j] << "   " << "R = " << R[i*M2_State.NVARS+j] << endl;
//         }
//     }
     
     // ******************************************************************************
     // The inverse of H_ONE satisfies: 
     // ==> H_ONE inv(H_ONE) = \delta_{ij}
     // Using the QR decomposition of H_ONE, we then have
     // Q R inv(H_ONE) = \delta_{ij}
     // This system is solved in two steps here:
     // First : Q A = \delta_{ij} where A = R inv(H_ONE)
     // And then : R inv(H_ONE) = A
     // ******************************************************************************
     
     
     // ******************************************************************************
     // Here we solve Q A = \delta_{ij} ==> A = inv(Q) \delta_{ij} = Q^{T} \delta_{ij}
     // In fact inv(Q) = Q^{T} since Q is an orthogonal matrix
     // ******************************************************************************
     transpose_matrix(M2_State.NVARS, M2_State.NVARS, H_ONE, Q); // now A ==> H_ONE
     
     // ******************************************************************************
     // Compute the product H_ONE * dIn_dInm1
     // ******************************************************************************
     for (int i = 0; i < M2_State.NVARS; i++) {
         for (int j = 0; j < M2_State.NVARS; j++) {
             Q[i*M2_State.NVARS + j] = H_ONE[i*M2_State.NVARS + j];
         }
     }
     
     for (int i = 0; i < M2_State.NVARS; i++) {
         for (int j = 0; j < M2_State.NVARS; j++) {
             H_ONE[i*M2_State.NVARS + j] = 0.0;
             for (int k = 0; k < M2_State.NVARS; k++) {
                 H_ONE[i*M2_State.NVARS + j] += Q[i*M2_State.NVARS + k]*dIn_dInm1(k, j, M2_State.I0, M2_State.N1_1, M2_State.N1_2, M2_State.N1_3, M2_State.gamma_1, M2_State.gamma_2);
             }
         }
     }
     
     // ******************************************************************************
     // Now we solve : R inv(H_ONE) = A using backward substitution since
     // R is an upper triangular matrix
     // ******************************************************************************
     Backward_Substitution(M2_State.NVARS, M2_State.NVARS, H_ONE, A_ONE, R);
     
//      for (int i = 0; i < M2_State.NVARS; i++) {
//          for (int j = 0; j < M2_State.NVARS; j++) {
//              cout << "i = " << i << "   " << "j = " << j << "   " << "A_ONE = " << A_ONE[i*M2_State.NVARS+j] << endl;
//         }
//     }
     
     delete[] H_ONE;
     delete[] Q;
     delete[] R;
 }
 
// ******************************************************************************
// This routine computes the matrix of second-derivatives of the Lagrange 
// multipliers with respect to the angular moments
// ===> A^{(2)}_{ijk} = {\partial^{2} \lam_i} {\partial I^{(j)} \partial I^{(k)}}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 void Compute_A_TWO_Matrix(long double *A_TWO, long double *A_ONE, const int &j, const int &l, M2_State_Param &M2_State) {
     long double H_TWO_imk;
     long double *temp_mat_M_ijl;
     temp_mat_M_ijl = new long double[M2_State.NVARS];
     
     // A^{(1)}_{ij} = Inv_Hessian_{ij}
     
     for (int i = 0; i < M2_State.NVARS; i++) {
         M2_State.index_Moment = i;
         temp_mat_M_ijl[i] = 0.0;
         for (int m = 0; m < M2_State.NVARS; m++) {
             for (int k = 0; k < M2_State.NVARS; k++) {
                 M2_State.index_Lag_i = m;
                 M2_State.index_Lag_j = k;
                 // Compute H^{(2)}_{imk}
                 Lebedev_Quadrature_Func(&H_TWO_imk, 1, H_TWO_Matrix_Entries, &M2_State);  
                 // Compute M_{ijl} =  H^{(1)}_{ip} A^{(2)}_{pjl} = - H^{(2)}_{imk} A^{(1)}_{mj} A^{(1)}_{kl}
                 temp_mat_M_ijl[i] += -H_TWO_imk*A_ONE[m*M2_State.NVARS+j]*A_ONE[k*M2_State.NVARS+l];
             }
         }
     }
     
     for (int k = 0; k < M2_State.NVARS; k++) {
         A_TWO[k] = 0.0;
         // Now compute A^{(2)}_{kjl} = inv(H^{(1)})_{ki} M_{ijl}
         for (int i = 0; i < M2_State.NVARS; i++) {
             A_TWO[k] += temp_mat_M_ijl[i]*A_ONE[k*M2_State.NVARS+i];
         }
     }
     
     delete[] temp_mat_M_ijl;
 }
  
void Compute_A_TWO_Matrix(long double *A_TWO, long double *A_ONE, M2_State_Param &M2_State) {
     long double H_TWO_ipq;
     long double *temp_mat_M_ijl;
     int index_temp;
     temp_mat_M_ijl = new long double[M2_State.NVARS*M2_State.NVARS*M2_State.NVARS];
     
     // A^{(1)}_{ij} = Inv_Hessian_{ij}
     
     for (int i = 0; i < M2_State.NVARS; i++) {
         M2_State.index_Moment = i;
         for (int j = 0; j < M2_State.NVARS; j++) {
             for (int k = 0; k < M2_State.NVARS; k++) {
                 index_temp = (i*M2_State.NVARS + j)*M2_State.NVARS + k;
                 temp_mat_M_ijl[index_temp] = 0.0;
                 for (int p = 0; p < M2_State.NVARS; p++) {
                     for (int q = 0; q < M2_State.NVARS; q++) {
                         M2_State.index_Lag_i = p;
                         M2_State.index_Lag_j = q;
                         // Compute H^{(2)}_{ipq}
                         Lebedev_Quadrature_Func(&H_TWO_ipq, 1, H_TWO_Matrix_Entries, &M2_State);  
                         // Compute M_{ijk} =  H^{(1)}_{ip} A^{(2)}_{pjk} = - H^{(2)}_{ipq} A^{(1)}_{pj} A^{(1)}_{qk}
                         temp_mat_M_ijl[index_temp] += -H_TWO_ipq*A_ONE[p*M2_State.NVARS+j]*A_ONE[q*M2_State.NVARS+k];
                    }
                 }
             }
         }
     }
     
     int index_temp_2;
     for (int i = 0; i < M2_State.NVARS; i++) {
         for (int j = 0; j < M2_State.NVARS; j++) {
             for (int k = 0; k < M2_State.NVARS; k++) {
                 index_temp = (i*M2_State.NVARS + j)*M2_State.NVARS + k;
                 A_TWO[index_temp] = 0.0;
                 // Now compute A^{(2)}_{ijk} = inv(H^{(1)})_{ip} M_{pjk}
                 for (int p = 0; p < M2_State.NVARS; p++) {
                     index_temp_2 = (p*M2_State.NVARS + j)*M2_State.NVARS + k;
                     A_TWO[index_temp] += A_ONE[i*M2_State.NVARS+p]*temp_mat_M_ijl[index_temp_2];
                }
            }
         }
     }
     
     delete[] temp_mat_M_ijl;
 }
 
// *************************************************************************************************
// This routine computes the matrix of third-derivatives of the Lagrange 
// multipliers with respect to the angular moments
// ===> A^{(2)}_{ijkl} = {\partial^{2} \lam_i} {\partial I^{(j)} \partial I^{(k)} \partial I^{(l)}}
// where \lam_j are the Lagrange multipliers
// *************************************************************************************************
 void Compute_A_THREE_Matrix(long double *A_THREE, const long double *A_ONE, const long double *A_TWO, const int &j, const int &k, const int &l, const int &p, M2_State_Param &M2_State) {
     long double d2_In_lower_dlambda_2, d3_In_lower_dlambda_3;
     long double *temp_mat_V_ijkl;
     temp_mat_V_ijkl = new long double[M2_State.NVARS];
     int index_temp_1, index_temp_2, index_temp_3;
     
     // A^{(1)}_{ij} = Inv_Hessian_{ij}
     
     for (int i = 0; i < M2_State.NVARS; i++) {
         M2_State.index_Moment = i;
         temp_mat_V_ijkl[i] = 0.0;
         for (int m = 0; m < M2_State.NVARS; m++) {
             for (int p = 0; p < M2_State.NVARS; p++) {
                 for (int q = 0; q < M2_State.NVARS; q++) {
                     M2_State.index_Lag_i = m;
                     M2_State.index_Lag_j = p;
                     M2_State.index_Lag_k = q; 
                     
                     // Compute H^{(3)}_{impq}
                     Lebedev_Quadrature_Func(&d3_In_lower_dlambda_3, 1, H_THREE_Matrix_Entries, &M2_State);  
                     // Compute M_{ijkl} = - H^{(3)}_{impq} A^{(1)}_{mj} A^{(1)}_{pk} A^{(1)}_{ql}
                     
                     temp_mat_V_ijkl[i] += -d3_In_lower_dlambda_3*A_ONE[m*M2_State.NVARS+j]*A_ONE[p*M2_State.NVARS+k]*A_ONE[q*M2_State.NVARS+l];
                 }
             }
         }
     }
     
     // Compute the matrix A^{(2)}: a third-order tensor
     // Compute_A_TWO_Matrix(A_TWO, A_ONE, M2_State);
         
     for (int i = 0; i < M2_State.NVARS; i++) {
         M2_State.index_Moment = i;
         
         for (int p = 0; p < M2_State.NVARS; p++) {
             for (int q = 0; q < M2_State.NVARS; q++) {
                 M2_State.index_Lag_i = p;
                 M2_State.index_Lag_j = q;
                 
                 // Compute H^{(2)}_{imp}
                 Lebedev_Quadrature_Func(&d2_In_lower_dlambda_2, 1, H_TWO_Matrix_Entries, &M2_State);  
                      
                 // Compute V_{ijkl} = M_{ijkl} - H^{(2)}_{ipq} (A^{(2)}_{qkl} A^{(1)}_{pj} + A^{(2)}_{pjl} A^{(1)}_{qk} + A^{(2)}_{pjk} A^{(1)}_{ql})
                 // where V_{ijkl} = H^{(1)}_{ip} A^{(3)}_{pjkl}
                 
                 index_temp_1 = (q*M2_State.NVARS + k)*M2_State.NVARS + l;
                 index_temp_2 = (p*M2_State.NVARS + j)*M2_State.NVARS + l;
                 index_temp_3 = (p*M2_State.NVARS + j)*M2_State.NVARS + k;
                 
                 temp_mat_V_ijkl[i] += -d2_In_lower_dlambda_2*(A_TWO[index_temp_1]*A_ONE[p*M2_State.NVARS+j] + A_TWO[index_temp_2]*A_ONE[q*M2_State.NVARS+k] + A_TWO[index_temp_3]*A_ONE[q*M2_State.NVARS+l]);
             }
         }
     }
     
     for (int i = 0; i < M2_State.NVARS; i++) {
         A_THREE[i] = 0.0;
         // Now compute A^{(3)}_{ijkl} = inv(H^{(1)})_{ip} V_{pjkl}
         for (int p = 0; p < M2_State.NVARS; p++) {
             A_THREE[i] += A_ONE[i*M2_State.NVARS+p]*temp_mat_V_ijkl[p];
         }
     }
     
     delete[] temp_mat_V_ijkl;
 }

// ******************************************************************************
// This routine computes the matrix of second-derivatives of the Third-order 
// angular moments with respect to the lower-order angular moments
// ===> {\partial^{2} I^{(3)}} {\partial I^{(j)} \partial I^{(k)}}
// ******************************************************************************
 long double Calculate_N3_Second_Derivatives(long double *dN3dlambda, long double *A_ONE, const int &j, const int &k, M2_State_Param &M2_State) {
     long double d2_N3_d_lambda_2;
     int index_temp;
     long double *A_TWO;
     long double d2N3dlIn = 0.0;
     A_TWO = new long double[M2_State.NVARS*M2_State.NVARS*M2_State.NVARS];
     
     Compute_A_TWO_Matrix(A_TWO, A_ONE, j, k, M2_State);
     
     for (int p = 0; p < M2_State.NVARS; p++) {
         for (int q = 0; q < M2_State.NVARS; q++) {
             M2_State.index_Lag_i = p;
             M2_State.index_Lag_j = q;
             
             // Compute {\partial^{2} I^{(3)}} {\partial \lam_{q} \partial \lam_{p}}
             Lebedev_Quadrature_Func(&d2_N3_d_lambda_2, 1, Higher_Order_Moments_Second_Derivatives, &M2_State);
             
             // Now compute {\partial^{2} I^{(3)}} {\partial I^{(k)} \partial I^{(j)}} = 
             // A^{(1)}_{pj} A^{(1)}_{qk} {\partial^{2} I^{(3)}} {\partial \lam_{q} \partial \lam_{p}}
             d2N3dlIn += d2_N3_d_lambda_2*A_ONE[p*M2_State.NVARS+j]*A_ONE[q*M2_State.NVARS+k];
        }
    }
    
    for (int p = 0; p < M2_State.NVARS; p++) {
        // Now compute {\partial^{2} I^{(3)}} {\partial I^{(k)} \partial I^{(j)}} = 
        // A^{(2)}_{pjk} A^{(1)}_{qk} {\partial I^{(3)}} {\partial \lam_{p}}
        index_temp = (p*M2_State.NVARS + j)*M2_State.NVARS + k;
        d2N3dlIn += dN3dlambda[k]*A_TWO[index_temp];
    }
     
     delete[] A_TWO;
     
     return d2N3dlIn;
 }
 
// *********************************************************************************
// This routine computes the matrix of third-derivatives of the Third-order 
// angular moments with respect to the lower-order angular moments
// ===> {\partial^{3} I^{(3)}} {\partial I^{(j)} \partial I^{(k)} \partial I^{(l)}}
// *********************************************************************************
 long double Calculate_N3_Third_Derivatives(long double *dN3dlambda, long double *A_ONE, const int &i, const int &j, const int &k, M2_State_Param &M2_State) {
     long double d_N3_d_lambda, d2_N3_d_lambda_2, d3_N3_d_lambda_3;
     long double *A_TWO, *A_THREE;
     long double d3N3dlIn = 0.0;
     int index_temp;
     A_TWO = new long double[M2_State.NVARS];
     A_THREE = new long double[M2_State.NVARS];
     
     Compute_A_TWO_Matrix(A_TWO, A_ONE, j, k, M2_State);
     
//      Compute_A_THREE_Matrix(A_THREE, A_ONE, j, k, M2_State);
     
     for (int p = 0; p < M2_State.NVARS; p++) {
         for (int q = 0; q < M2_State.NVARS; q++) {
             for (int r = 0; r < M2_State.NVARS; r++) {
                 M2_State.index_Lag_i = p;
                 M2_State.index_Lag_j = q;
                 M2_State.index_Lag_k = r;
                 
                 Lebedev_Quadrature_Func(&d3_N3_d_lambda_3, 1, Higher_Order_Moments_Third_Derivatives, &M2_State);
                 
                 d3N3dlIn += d3_N3_d_lambda_3*A_ONE[p*M2_State.NVARS+i]*A_ONE[q*M2_State.NVARS+j]*A_ONE[r*M2_State.NVARS+k];
             }
        }
    }
    
    for (int p = 0; p < M2_State.NVARS; p++) {
         for (int q = 0; q < M2_State.NVARS; q++) {
             M2_State.index_Lag_i = p;
             M2_State.index_Lag_j = q;
             
             // Compute {\partial^{2} I^{(3)}} {\partial \lam_{q} \partial \lam_{p}}
             Lebedev_Quadrature_Func(&d2_N3_d_lambda_2, 1, Higher_Order_Moments_Second_Derivatives, &M2_State);
             
             // Now compute {\partial^{2} I^{(3)}} {\partial I^{(k)} \partial I^{(j)}} = 
             // A^{(1)}_{pj} A^{(1)}_{qk} {\partial^{2} I^{(3)}} {\partial \lam_{q} \partial \lam_{p}}
             d3N3dlIn += d2_N3_d_lambda_2*A_ONE[p*M2_State.NVARS+j]*A_ONE[q*M2_State.NVARS+k];
        }
    }
    
    
    
    for (int p = 0; p < M2_State.NVARS; p++) {
        index_temp = ((p*M2_State.NVARS + i)*M2_State.NVARS + j)*M2_State.NVARS + k;
        d3N3dlIn += dN3dlambda[k]*A_THREE[index_temp];
    }
     
     delete[] A_TWO;
     delete[] A_THREE;
     
     return d3N3dlIn;
 }
 
 void transpose_matrix(const int &nl, const int &nc, long double *A_transpose, const long double *A) {
     // This routine aims at storing elements of the transpose of 
     // matrix A in the matrix A_transpose
     
     if (nl != nc) {
         cout << "Cannot take transpose of matrix as it is not square" << endl; 
         exit(0);
     }
     
     for (int i = 0; i < nl; i++) {
         for (int j = 0; j < nc; j++) {
             A_transpose[i*nc + j] = A[j*nc + i];
         }
     }
 }
 
 void Backward_Substitution(const int &nl, const int &nc, long double *A, long double *Q, const long double *R) {
     // This routine aims at solving the system of linear equations: R Q = A for Q
     // where R is an upper triangular matrix
     // This can be achieved using backward substitution
     
     // nl: number of lines
     // nc: number of columns
     
     /* Backward substitution for discovering values of unknowns */
     for(int k = 0; k < nc; k++) {                   
         for(int i = nl-1; i >= 0; i--) {                     
             Q[i*nc + k] = A[i*nc + k];
             for(int j = i+1; j < nl;j++) {
                 if(i != j) {
                     Q[i*nc + k] = Q[i*nc + k] - R[i*nl + j]*Q[j*nc + k];
                }          
            }
            Q[i*nc + k] = Q[i*nc + k]/R[i*nl + i];
        }
     }
 }
 
// ******************************************************************************
// This routine computes derivatives, up to third-order, of the highest-order 
// angular moments with respect to the lower-order angular moments
// ******************************************************************************
 void Calculate_Higher_Order_Moments_Derivatives(record_N3 &rec_N3, const long double *x, const long double *Sk, M2_State_Param &M2_State) {
     long double norm_f_2;
     long double *Q_data, *Q_data_test, *Higher_order_Q_datas;
     long double dN3_111_dN2_11, dN3_122_dN2_11, dN3_123_dN2_11;
     long double dN3_111_dN2_12, dN3_122_dN2_12, dN3_123_dN2_12;
     long double dN3_111_dN2_22, dN3_122_dN2_22, dN3_123_dN2_22;
     long double dI3_111_dI0, dI3_122_dI0;
     Q_data = new long double[M2_State.NVARS*M2_State.NVARS];
     Q_data_test = new long double[M2_State.NVARS*M2_State.NVARS];
     Higher_order_Q_datas = new long double[M2_State.NVARS];
     long double d2_N3_122_dN2_11_dN1_1, d2_N3_122_dN2_22_dN1_1;
     long double d3_N3_122_dN2_11_dN2_22_dN1_1;
     
     M2_State.set_x(x);
     M2_State.set_Sk(Sk);
     M2_State.set_Sk_cur(Sk);
     M2_State.set_ak(Sk);
     
     Compute_A_ONE_Matrix(Q_data, x, Sk, M2_State);
     Compute_H_ONE_Matrix(Q_data_test, x, Sk, M2_State);
     
//      long double temp_val_Moments;
//      long double diff_temp_val_Moments;
//      for (int i = 0; i < M2_State.NVARS; i++) {
//          for (int j = 0; j < M2_State.NVARS; j++) {
//              temp_val_Moments = 0.0;
//              for (int k = 0; k < M2_State.NVARS; k++) {
//                  temp_val_Moments += Q_data_test[j*M2_State.NVARS+k]*Q_data[k*M2_State.NVARS+i];
//             }
//             
// //             if (i == j) {
// //                 if (fabs(1.0 - temp_val_Moments) > 1.0e-8) {
// //                     cout << "i = " << i << "   " << "j = " << j << "   " << "temp_val_Moments = " << diff_temp_val_Moments << endl;
// //                     exit(0);
// //                 }
// //             } else {
// //                 if (fabs(temp_val_Moments) > 1.0e-4) {
// //                     cout << "i = " << i << "   " << "j = " << j << "   " << "temp_val_Moments = " << diff_temp_val_Moments << endl;
// //                     exit(0);
// //                 }
// //             }
//             
//             diff_temp_val_Moments = fabs(temp_val_Moments - dIn_dInm1(M2_State.Index[j], M2_State.Index[i], M2_State.I0, M2_State.N1_1, M2_State.N1_2, M2_State.N1_3, M2_State.gamma_1, M2_State.gamma_2));
//             
//             if (diff_temp_val_Moments > 1.0e-6) {
//                 cout << "i = " << i << "   " << "j = " << j << "   " << "temp_val_Moments = " << diff_temp_val_Moments << endl;
//                 exit(0);
//             }
//         }
//      }
     
     norm_f_2 = pow(M2_State.N1_1, 2) + pow(M2_State.N1_2, 2) + pow(M2_State.N1_3, 2);
        
     M2_State.Index_Higher_Order_Moments = 0;
     
     dI3_111_dI0 = 0.0;
     rec_N3.dN3_111_dN1_1 = 0.0;
     rec_N3.dN3_111_dN1_2 = 0.0;
     rec_N3.dN3_111_dN1_3 = 0.0;
     
     dN3_111_dN2_11 = 0.0;
     dN3_111_dN2_12 = 0.0;
     dN3_111_dN2_22 = 0.0;
     
     switch (M2_State.Domain_Type) {
         case GAM1_GAM2_GAM3:
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 dI3_111_dI0 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 0];
                 rec_N3.dN3_111_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
                 rec_N3.dN3_111_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
                 rec_N3.dN3_111_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_111_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_111_dN2_12 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 5];
//                  dN3_111_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
            }
            break;
         case BOUNDARY_GAM1:
             cout << "Update computation of the derivatives here !!!!!!!!!!!!!!!!!!!" << endl;
             // N3_111 = N1_1^3
             rec_N3.dN3_111_dN1_1 = 3.0*pow(M2_State.N1_1, 2);
             rec_N3.dN3_111_dN1_2 = 0.0;
             rec_N3.dN3_111_dN1_3 = 0.0;
             
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 dN3_111_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
             }
             break;
         case BOUNDARY_GAM2:
             cout << "Update computation of the derivatives here !!!!!!!!!!!!!!!!!!!" << endl;
             // N3_111 ==> maximum entropy
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 rec_N3.dN3_111_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
                 rec_N3.dN3_111_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  dN3_111_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
            }
            break;
         case BOUNDARY_GAM3:
             cout << "Update computation of the derivatives here !!!!!!!!!!!!!!!!!!!" << endl;
             // gamma_3 = 0 ==> gamma_1 + gamma_2 = 1 and N2_33 = (N1_3)^2;
             //              ==> N2_11 + N2_22 = 1 - (N1_3)^2
             //              ==> N2_22 = 1 - N2_11 - (N1_3)^2
             // N3_111 ==> maximum entropy
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 rec_N3.dN3_111_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
//                  rec_N3.dN3_111_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  dN3_111_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
            }
            break;
         default:
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 rec_N3.dN3_111_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
                 rec_N3.dN3_111_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
                 rec_N3.dN3_111_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_111_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_111_dN2_12 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 5];
//                  dN3_111_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
            }
             break;
    };
     
    if (M2_State.display) {
        if (M2_State.id_proc == M2_State.proc_display) {
            cout << "dI3_111_dI0 = " << dI3_111_dI0 << "  " << "dN3_111_dN1_1 = " << rec_N3.dN3_111_dN1_1 << "  " << "dN3_111_dN1_2 = " << rec_N3.dN3_111_dN1_2 << "  " << "dN3_111_dN1_3 = " << rec_N3.dN3_111_dN1_3<< "  " << "dN3_111_dN2_11 = " << dN3_111_dN2_11 << "  " << "dN3_111_dN2_12 = " << dN3_111_dN2_12 << "  " << "dN3_111_dN2_22 = " << dN3_111_dN2_22 << endl;
        }
    }
     
     M2_State.Index_Higher_Order_Moments = 1;
     
     dI3_122_dI0 = 0.0;
     rec_N3.dN3_122_dN1_1 = 0.0;
     rec_N3.dN3_122_dN1_2 = 0.0;
     rec_N3.dN3_122_dN1_3 = 0.0;
     
     dN3_122_dN2_11 = 0.0;
     dN3_122_dN2_12 = 0.0;
     dN3_122_dN2_22 = 0.0;
     
     switch (M2_State.Domain_Type) {
         case GAM1_GAM2_GAM3:
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 dI3_122_dI0 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 0];
                 rec_N3.dN3_122_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
                 rec_N3.dN3_122_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
                 rec_N3.dN3_122_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_122_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_122_dN2_12 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 5];
//                  dN3_122_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
            }
            
            
            if (M2_State.finite_diff_Triangle) {
                d2_N3_122_dN2_11_dN1_1 = Calculate_N3_Second_Derivatives(Higher_order_Q_datas, Q_data, 4, 1, M2_State);
                d2_N3_122_dN2_22_dN1_1 = Calculate_N3_Second_Derivatives(Higher_order_Q_datas, Q_data, 7, 1, M2_State);
                
                d3_N3_122_dN2_11_dN2_22_dN1_1 = Calculate_N3_Third_Derivatives(Higher_order_Q_datas, Q_data, 4, 7, 1, M2_State);
            }
            break;
         case BOUNDARY_GAM1:
             // N3_122 = N1_1 N2_22 = N1_1 [(N1_2)^2 + gam2 (1 - norm_f^2)]
             rec_N3.dN3_122_dN1_1 = pow(M2_State.N1_2, 2) + M2_State.gamma_2*(1.0 - norm_f_2 - 2.0*pow(M2_State.N1_1, 2));
             rec_N3.dN3_122_dN1_2 = 2.0*M2_State.N1_1*(M2_State.N1_2 - M2_State.gamma_2*M2_State.N1_2);
             rec_N3.dN3_122_dN1_3 = -2.0*M2_State.gamma_2*M2_State.N1_1*M2_State.N1_3;
             dN3_122_dN2_22 = M2_State.N1_1;
             break;
         case BOUNDARY_GAM2:
             // N3_122 = N1_1 (N1_2)^2
             rec_N3.dN3_122_dN1_1 = pow(M2_State.N1_2, 2);
             rec_N3.dN3_122_dN1_2 = 2.0*M2_State.N1_1*M2_State.N1_2;
             rec_N3.dN3_122_dN1_3 = 0.0;
             // dN3_122_dN2_11 = M2_State.N1_1;
             
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 dN3_122_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
             }
             break;
         case BOUNDARY_GAM3:
             cout << "Update computation of the derivatives here !!!!!!!!!!!!!!!!!!!" << endl;
             // N3_122 ==> maximum entropy
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 rec_N3.dN3_122_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
                 rec_N3.dN3_122_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  dN3_122_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
            }
//             rec_N3.d2N3_122_dN1_2_dN1_2 = Calculate_N3_Second_Derivatives(Higher_order_Q_datas, Q_data, 2, 2, M2_State);
            break;
         default:
             Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
             for (int i = 0; i < M2_State.NVARS; i++) {
                 rec_N3.dN3_122_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
                 rec_N3.dN3_122_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
                 rec_N3.dN3_122_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_122_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_122_dN2_12 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 5];
//                  dN3_122_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
            }
            break;
    };
    
    if (M2_State.display) {
        if (M2_State.id_proc == M2_State.proc_display) {
            cout << "dI3_122_dI0 = " << dI3_122_dI0 << "  " << "dN3_122_dN1_1 = " << rec_N3.dN3_122_dN1_1 << "  " << "dN3_122_dN1_2 = " << rec_N3.dN3_122_dN1_2 << "  " << "dN3_122_dN1_3 = " << rec_N3.dN3_122_dN1_3 << "  " << "dN3_122_dN2_11 = " << dN3_122_dN2_11 << "  " << "dN3_122_dN2_12 = " << dN3_122_dN2_12 << "  " << "dN3_122_dN2_22 = " << dN3_122_dN2_22 << endl;
        }
    }
     
//      M2_State.Index_Higher_Order_Moments = 2;
//      Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
//      
//      rec_N3.dN3_123_dN1_1 = 0.0;
//      rec_N3.dN3_123_dN1_2 = 0.0;
//      rec_N3.dN3_123_dN1_3 = 0.0;
//      
//      dN3_123_dN2_11 = 0.0;
//      dN3_123_dN2_22 = 0.0;
//      
//      switch (M2_State.Domain_Type) {
//          case GAM1_GAM2_GAM3:
//              Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
//              for (int i = 0; i < M2_State.NVARS; i++) {
//                  rec_N3.dN3_123_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
//                  rec_N3.dN3_123_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  rec_N3.dN3_123_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_123_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_123_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
//             }
//             break;
//          case BOUNDARY_GAM1:
//              // N3_123 = N1_1 N1_2 N1_3
//              rec_N3.dN3_123_dN1_1 = M2_State.N1_2*M2_State.N1_3;
//              rec_N3.dN3_123_dN1_2 = M2_State.N1_1*M2_State.N1_3;
//              rec_N3.dN3_123_dN1_3 = M2_State.N1_1*M2_State.N1_2;
//              break;
//          case BOUNDARY_GAM2:
//              // N3_123 = N1_1 N1_2 N1_3
//              rec_N3.dN3_123_dN1_1 = M2_State.N1_2*M2_State.N1_3;
//              rec_N3.dN3_123_dN1_2 = M2_State.N1_1*M2_State.N1_3;
//              rec_N3.dN3_123_dN1_3 = M2_State.N1_1*M2_State.N1_2;
//              break;
//          case BOUNDARY_GAM3:
//              // N3_123 = N1_1 N1_2 N1_3
//              rec_N3.dN3_123_dN1_1 = M2_State.N1_2*M2_State.N1_3;
//              rec_N3.dN3_123_dN1_2 = M2_State.N1_1*M2_State.N1_3;
//              rec_N3.dN3_123_dN1_3 = M2_State.N1_1*M2_State.N1_2;
//             break;
//          default:
//              Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
//              for (int i = 0; i < M2_State.NVARS; i++) {
//                  rec_N3.dN3_123_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
//                  rec_N3.dN3_123_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  rec_N3.dN3_123_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_123_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_123_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
//             }
//             break;
//     };
//     
//     
//     M2_State.Index_Higher_Order_Moments = 3;
//     long double dI3_112_dI0, dN3_112_dN1_1, dN3_112_dN1_2, dN3_112_dN1_3, dN3_112_dN2_11, dN3_112_dN2_12, dN3_112_dN2_22;
//     dI3_112_dI0 = 0.0;
//     dN3_112_dN1_1 = 0.0;
//     dN3_112_dN1_2 = 0.0;
//     dN3_112_dN1_3 = 0.0;
//     
//     dN3_112_dN2_11 = 0.0;
//     dN3_112_dN2_12 = 0.0;
//     dN3_112_dN2_22 = 0.0;
//      
//      switch (M2_State.Domain_Type) {
//          case GAM1_GAM2_GAM3:
//              Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
//              for (int i = 0; i < M2_State.NVARS; i++) {
//                  dI3_112_dI0 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 0];
//                  dN3_112_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
//                  dN3_112_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  dN3_112_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_112_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_112_dN2_12 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 5];
//                  dN3_112_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
//             }
//             break;
//          case BOUNDARY_GAM1:
//              
//              break;
//          case BOUNDARY_GAM2:
//              
//              break;
//          case BOUNDARY_GAM3:
//              
//             break;
//          default:
//              Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
//              for (int i = 0; i < M2_State.NVARS; i++) {
//                  dI3_112_dI0 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 0];
//                  dN3_112_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
//                  dN3_112_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  dN3_112_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_112_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_112_dN2_12 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 5];
//                  dN3_112_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
//             }
//             break;
//     };
//     
//     if (M2_State.display) {
//         if (M2_State.id_proc == M2_State.proc_display) {
//             cout << "dI3_112_dI0 = " << dI3_112_dI0 << "  " << "dN3_112_dN1_1 = " << dN3_112_dN1_1 << "  " << "dN3_112_dN1_2 = " << dN3_112_dN1_2 << "  " << "dN3_112_dN1_3 = " << dN3_112_dN1_3 << "  " << "dN3_112_dN2_11 = " << dN3_112_dN2_11 << "  " << "dN3_112_dN2_12 = " << dN3_112_dN2_12 << "  " << "dN3_112_dN2_22 = " << dN3_112_dN2_22 << endl;
//         }
//     }
//     
//     M2_State.Index_Higher_Order_Moments = 4;
//     long double dI3_222_dI0, dN3_222_dN1_1, dN3_222_dN1_2, dN3_222_dN1_3, dN3_222_dN2_11, dN3_222_dN2_12, dN3_222_dN2_22;
//     dI3_222_dI0 = 0.0;
//     dN3_222_dN1_1 = 0.0;
//     dN3_222_dN1_2 = 0.0;
//     dN3_222_dN1_3 = 0.0;
//     
//     dN3_222_dN2_11 = 0.0;
//     dN3_222_dN2_12 = 0.0;
//     dN3_222_dN2_22 = 0.0;
//      
//      switch (M2_State.Domain_Type) {
//          case GAM1_GAM2_GAM3:
//              Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
//              for (int i = 0; i < M2_State.NVARS; i++) {
//                  dI3_222_dI0 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 0];
//                  dN3_222_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
//                  dN3_222_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  dN3_222_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_222_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_222_dN2_12 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 5];
//                  dN3_222_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
//             }
//             break;
//          case BOUNDARY_GAM1:
//              
//              break;
//          case BOUNDARY_GAM2:
//              
//              break;
//          case BOUNDARY_GAM3:
//              
//             break;
//          default:
//              Lebedev_Quadrature_Func(Higher_order_Q_datas, M2_State.NVARS, Higher_Order_Moments_Derivatives, &M2_State);
//              for (int i = 0; i < M2_State.NVARS; i++) {
//                  dI3_222_dI0 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 0];
//                  dN3_222_dN1_1 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 1];
//                  dN3_222_dN1_2 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 2];
//                  dN3_222_dN1_3 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 3];
//                  dN3_222_dN2_11 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 4];
//                  dN3_222_dN2_12 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 5];
//                  dN3_222_dN2_22 += Higher_order_Q_datas[i]*Q_data[i*M2_State.NVARS + 7];
//             }
//             break;
//     };
//     
//     if (M2_State.display) {
//         if (M2_State.id_proc == M2_State.proc_display) {
//             cout << "dI3_222_dI0 = " << dI3_222_dI0 << "  " << "dN3_222_dN1_1 = " << dN3_222_dN1_1 << "  " << "dN3_222_dN1_2 = " << dN3_222_dN1_2 << "  " << "dN3_222_dN1_3 = " << dN3_222_dN1_3 << "  " << "dN3_222_dN2_11 = " << dN3_222_dN2_11 << "  " << "dN3_222_dN2_12 = " << dN3_222_dN2_12 << "  " << "dN3_222_dN2_22 = " << dN3_222_dN2_22 << endl;
//         }
//     }
    
    // For any given set of values for N1_1, N1_2, and N1_3, we have
    // (we can only assume independency of independent variables)
    // N2_11 = N1_1^2 + gam1*(1 - norm_f^2)
    // N2_22 = N1_2^2 + gam2*(1 - norm_f^2)
    // Therefore, considering the set (N1_1, N1_2, N1_3, gam1, gam2)
    // df/dN2_11 = (df/d_N1_1) (d_N1_1/dN2_11) + (df/d_N1_2) (d_N1_2/dN2_11) + ...
    //             (df/d_N1_3) (d_N1_3/dN2_11) + (df/dgam1) (dgam1/dN2_11)
    //
    // df/dN2_22 = (df/d_N1_1) (d_N1_1/dN2_22) + (df/d_N1_2) (d_N1_2/dN2_22) + ...
    //             (df/d_N1_3) (d_N1_3/dN2_22) + (df/dgam2) (dgam2/dN2_22)
    // So
    // (df/dgam1) = (df/dN2_11) * (dN2_11/dgam1) ==> (df/dgam1) = (df/dN2_11) * (1 - norm_f^2)
    // (df/dgam2) = (df/dN2_22) * (dN2_22/dgam2) ==> (df/dgam2) = (df/dN2_22) * (1 - norm_f^2))
    
    if (M2_State.Domain_Type == GAM1_GAM2_GAM3) {
        rec_N3.dN3_111_dgam1 = dN3_111_dN2_11 * (1.0 - norm_f_2);
        rec_N3.dN3_122_dgam1 = dN3_122_dN2_11 * (1.0 - norm_f_2);
        rec_N3.dN3_123_dgam1 = dN3_123_dN2_11 * (1.0 - norm_f_2);
        
        rec_N3.dN3_111_dgam2 = dN3_111_dN2_22 * (1.0 - norm_f_2);
        rec_N3.dN3_122_dgam2 = dN3_122_dN2_22 * (1.0 - norm_f_2);
        rec_N3.dN3_123_dgam2 = dN3_123_dN2_22 * (1.0 - norm_f_2);
        
        if (M2_State.finite_diff_Triangle) {
            rec_N3.d2_N3_122_dgam1_dN1_1 = d2_N3_122_dN2_11_dN1_1 * (1.0 - norm_f_2);
            rec_N3.d2_N3_122_dgam2_dN1_1 = d2_N3_122_dN2_22_dN1_1 * (1.0 - norm_f_2);
        }
    }
    
    if (M2_State.display) {
        if (M2_State.id_proc == M2_State.proc_display) {
            cout << "********************************************************************" << endl;
            cout << endl;
        }
    }
     
     delete[] Q_data;
     delete[] Q_data_test;
     delete[] Higher_order_Q_datas;
 }

 long double dIn_dInm1(const int &i, const int &j, const long double &I0, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
     long double dIn;
     dIn = 0.0;
     
     if (i == j) {
         dIn = 1.0;
     }
     
     switch (i) {
         case 0:
             switch (j) {
                 case 0:
                     dIn = 1.0;
                     break;
            };
            break;
         case 1:
             switch (j) {
                 case 1:
                     dIn = 1.0;
                     break;
            };
            break;
         case 2:
             switch (j) {
                 case 2:
                     dIn = 1.0;
                     break;
            };
            break;
         case 3:
             switch (j) {
                 case 3:
                     dIn = 1.0;
                     break;
            };
            break;
         case 4:
             switch (j) {
                 case 0:
                     dIn = -pow(N1_1,2) + gam1*(1.0 + (pow(N1_1,2) + pow(N1_2,2) + pow(N1_3,2)));
                     break;
                 case 1:
                     dIn = 2.0*N1_1*(1.0 - gam1);
                     break;
                 case 2:
                     dIn = -2.0*N1_2*gam1;
                     break;
                 case 3:
                     dIn = -2.0*N1_3*gam1;
                     break;
                 case 4:
                     dIn = 1.0;
                     break;
            };
            break;
         case 5:
             switch (j) {
                 case 0:
                     dIn = -N1_1*N1_2;
                     break;
                 case 1:
                     dIn = N1_2;
                     break;
                 case 2:
                     dIn = N1_1;
                     break;
                 case 5:
                     dIn = 1.0;
                     break;
            };
            break;
         case 6:
             switch (j) {
                 case 0:
                     dIn = -N1_1*N1_3;
                     break;
                 case 1:
                     dIn = N1_3;
                     break;
                 case 3:
                     dIn = N1_1;
                     break;
                 case 6:
                     dIn = 1.0;
                     break;
            };
            break;
         case 7:
             switch (j) {
                 case 0:
                     dIn = -pow(N1_2,2) + gam2*(1.0 + (pow(N1_1,2) + pow(N1_2,2) + pow(N1_3,2)));
                     break;
                 case 1:
                     dIn = -2.0*N1_1*gam2;
                     break;
                 case 2:
                     dIn = 2.0*N1_2*(1.0 - gam2);
                     break;
                 case 3:
                     dIn = -2.0*N1_3*gam2;
                     break;
                 case 7:
                     dIn = 1.0;
                     break;
            };
            break;
         case 8:
             switch (j) {
                 case 0:
                     dIn = -N1_2*N1_3;
                     break;
                 case 2:
                     dIn = N1_3;
                     break;
                 case 3:
                     dIn = N1_2;
                     break;
                 case 8:
                     dIn = 1.0;
                     break;
            };
            break;
     };
     return dIn;
 }
