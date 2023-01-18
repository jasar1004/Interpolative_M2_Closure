/*******************************************************************
  File: M1_Model_1D_Utilities.cc

  Description:  ...  

  Author:  Joachim A.R. Sarr

  Date:    May 05th, 2020
*******************************************************************/

#include "M2_Model_3D_Utilities.h"

// int binomialCoefficients(const long double &n, const long double &k) {
//    if (k == 0 || k == n)
//    return 1;
//    return binomialCoefficients(n - 1, k - 1) + binomialCoefficients(n - 1, k);
// }

long double binomialCoefficients(const long double &n, const long double &k) {
    long double c;
    if (k < 0) {
        return 0;
    }
    
    if (k == 0 || k == n) {
        return 1;
    }
    
    c = 1.0;
    for (int i = 0; i < k; i++) {
        c = c * (n - i) / (i + 1);
    }
    return c;
}

// long double Jacobi_Polynomials(const long double &x, const int &l, const int &alpha, const int &beta) {
//     long double B_a_b = 0.0;
//     long double Coeff_1, Coeff_2;
//     for (int i = 0; i < l; i++) {
//         Coeff_1 = binomialCoefficients(l+alpha, l - i);
//         Coeff_2 = binomialCoefficients(l+beta, i);
//         B_a_b += (long double)(Coeff_1)*(long double)(Coeff_2)*pow(x, i)*pow(x - 1.0, l - i);
//     }
//     return B_a_b;
// }
// 
// long double Precompute_Legendre_Basis_Coeffs(const int &n, const int &k) {
//     long double Binomial_Coeffs;
//     long double factorial_ratio;
//     long double Leg_Basis_Coeffs;
//     
//     Leg_Basis_Coeffs = pow(2, n);
//     
//     Binomial_Coeffs = binomialCoefficients(n, k);
//     Leg_Basis_Coeffs *= Binomial_Coeffs;
//     
//     Binomial_Coeffs = binomialCoefficients((long double)((n+k-1)/2.0), n);
//     Leg_Basis_Coeffs *= Binomial_Coeffs;
//     
//     return Leg_Basis_Coeffs;
// }
// 
// long double Legendre_Polynomial(const long double &x, const int &n) {
//     long double Coeff_Leg_Poly;
//     long double Leg_Poly = 0.0;
//     for (int k = n-1; k >= 0; k++) {
//         Coeff_Leg_Poly = Precompute_Legendre_Basis_Coeffs(n, k);
//         if (k == n - 1) {
//             Leg_Poly = Coeff_Leg_Poly;
//         } else {
//             Leg_Poly = Coeff_Leg_Poly + x*Leg_Poly;
//         }
//     }
//     return Leg_Poly;
// }
// 
// long double Proriol_Polynomial_Basis(const int &k, const int &l, const long double &gam1, const long double &gam2) {
//     long double gam1_square, gam2_square;
//     long double Leg_gam1_p, jacob_gam2_p;
//     long double poly;
//     gam1_square = 2.0*gam1/(1.0 - gam2) - 1.0;
//     gam2_square = 2.0*gam2 - 1.0;
//     
//     Leg_gam1_p = Legendre_Polynomial(gam1, k);
//     jacob_gam2_p = Jacobi_Polynomials(gam2, l, 2*k + 1, 0);
//     
//     poly = pow(1.0 - gam2, k)*Leg_gam1_p*jacob_gam2_p;
//     return poly;
// }

int factorial(const int &n) { 
    if (n == 0) 
        return 1; 
    return n * factorial(n - 1); 
} 

int factorial(const int &n, const int &a) { 
    if (n == a) 
        return 1; 
    return n * factorial(n - 1, a); 
} 

long double factorial_inv(const int &n, const int &a) { 
    if (n == a) 
        return 1.0; 
    return (long double)(1.0/n) * factorial_inv(n - 1, a); 
} 


long double factorial_ratios(const int &num, const int &den) { 
    long double ratio;
//         cout << "num = " << num << "   " << "den = " << den << endl;
    if (num > den) { 
        ratio = (long double)(factorial(num, den));
    } else if (num < den) {
        ratio = (long double)(factorial_inv(den, num));
//         cout << "ratio = " << ratio << endl;
    } else {
        ratio = 1.0;
    }
    
    return ratio;
}

int double_factorial(const int &n) { 
    if (n == 0 || n == 1) 
        return 1; 
    return n * double_factorial(n - 2); 
} 

long double SH_Normalization_Constant(const int &l, const int &m) {
    long double K_l_m;
    K_l_m = (2.0*l+1.0)/(4.0*PI);
//     K_l_m *= (long double)(factorial(l - m))/(long double)(factorial(l + m));
    K_l_m *= factorial_ratios(l - m, l + m);
    K_l_m = sqrt(K_l_m);
    
//     cout << "factorial_ratios(l - m, l + m) = " << factorial_ratios(l - m, l + m) << endl;
    
    return K_l_m;
}

long double P_l_m_Polynomials(const long double &x, const int &l, const int &m) {
    long double P_l_m;
    switch (l) {
        case 0:
            switch (m) {
                case 0:
                    P_l_m = 1.0;
                    break;
            };
            break;
        default:
            if (l == m) {
                P_l_m = pow(-1.0, m)*double_factorial(2.0*m-1.0)*pow((1.0 - pow(x, 2)), m/2.0);
            } else if (l == m+1) {
                P_l_m = x*(2.0*m+1.0)*P_l_m_Polynomials(x, l-1, m);
            }
            break;
    };
    return P_l_m;
}

long double Associated_Legendre_Polynomials(const long double &x, const int &l, const int &m) {
    long double P_l_m;
    long double P_lminus1_m, P_lminus2_m;
    long double coeff1, coeff2;
    
    if (l == fabs(m) || l == fabs(m)+1) {
        P_l_m = P_l_m_Polynomials(x, l, fabs(m));
    } else {
        P_lminus2_m = Associated_Legendre_Polynomials(x, l - 2, fabs(m));
        P_lminus1_m = Associated_Legendre_Polynomials(x, l - 1, fabs(m));
        
        coeff1 = x*(long double)(2.0*l - 1.0)/(long double)(l - fabs(m));
        coeff2 = (long double)(l + fabs(m) - 1.0)/(long double)(l - fabs(m));
            
        P_l_m = coeff1*P_lminus1_m - coeff2*P_lminus2_m;
    }
    
    return P_l_m;
}

long double Precompute_SH_Basis_Coeffs(const int &l, const int &m, const int &i, const int &j, const int &k) {
    long double Binomial_Coeffs;
    long double factorial_ratio;
    long double SH_Basis_Coeffs;
    
    if (m == 0) {
        SH_Basis_Coeffs = pow(2,l);
        SH_Basis_Coeffs *= sqrt((2.0*l+1.0)/(4.0*PI));
        
        Binomial_Coeffs = binomialCoefficients(l, k);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)((l+k-1)/2.0), l);
        SH_Basis_Coeffs *= Binomial_Coeffs;
    } else if (m > 0) {
        SH_Basis_Coeffs = m*pow(-1.0,k+j)*pow(2.0, l+m-2*j-1);
        SH_Basis_Coeffs *= sqrt(2.0)*sqrt((2.0*l+1.0)/(4.0*PI));
        
        factorial_ratio = factorial_ratios(l - m, l + m);
        SH_Basis_Coeffs *= sqrt(factorial_ratio);
        
        factorial_ratio = factorial_ratios(i, i - m);
        SH_Basis_Coeffs *= factorial_ratio;
        
        factorial_ratio = factorial_ratios(m-j-1, m-2*j);
        SH_Basis_Coeffs *= factorial_ratio/factorial(j);
        
        Binomial_Coeffs = binomialCoefficients(l, i);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)((l+i-1)/2.0), l);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients(j, k);
        SH_Basis_Coeffs *= Binomial_Coeffs;
    } else {
        SH_Basis_Coeffs = pow(-1.0,k+j)*pow(2.0, l+fabs(m)-2*j-1);
        SH_Basis_Coeffs *= sqrt(2.0)*sqrt((2.0*l+1.0)/(4.0*PI));
        
        factorial_ratio = factorial_ratios(l - fabs(m), l + fabs(m));
        SH_Basis_Coeffs *= sqrt(factorial_ratio);
        
        factorial_ratio = factorial_ratios(i, i - fabs(m));
        SH_Basis_Coeffs *= factorial_ratio;
        
        Binomial_Coeffs = binomialCoefficients((long double)(l), i);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)((l+i-1)/2.0), l);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)(fabs(m)-j-1), j);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)(j), k);
        SH_Basis_Coeffs *= Binomial_Coeffs;
    }
    
    return SH_Basis_Coeffs;
}

long double Precompute_First_Kind_Chebyshev_Basis_Coeffs(const int &n, const int &k) {
    long double factorial_ratio;
    long double Chebyshev_Basis_Coeffs;
    
    if (n == 0) {
        Chebyshev_Basis_Coeffs = 1.0;
    } else {
        Chebyshev_Basis_Coeffs = (n/2.0)*pow(-1.0,k)*pow(2.0, n-2*k);
        
        factorial_ratio = factorial_ratios(n-k-1, n-2*k);
        Chebyshev_Basis_Coeffs *= factorial_ratio/(long double)(factorial(k));
    }
    
    return Chebyshev_Basis_Coeffs;
}

void Chebyshev_First_Kind_to_Monomial_Basis_ratio_E(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &N_Coeffs_Triangle_gam1_gam2) {
    long double Coeffs_Cheby;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_SH*N_Coeffs_Triangle_gam1_gam2;
    int N_Coeffs_Total_Without_SH = N_Coeffs_E*N_Coeffs_f*N_Coeffs_Triangle_gam1_gam2;
    int index_Monomial_Basis, index_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
        for (int i_fit_SH = 0; i_fit_SH < N_Coeffs_SH; i_fit_SH++) {
            for (int i_fit_gam1_gam2 = 0; i_fit_gam1_gam2 < N_Coeffs_Triangle_gam1_gam2; i_fit_gam1_gam2++) {
                
                for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
                    index_Monomial_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_Triangle_gam1_gam2 + i_fit_gam1_gam2;
                    
                    for (int i_fit_E_Orthog = 0; i_fit_E_Orthog < N_Coeffs_E; i_fit_E_Orthog++) {
                        index_Orthog_Basis = (i_fit_E_Orthog*N_Coeffs_f + i_fit_f)*N_Coeffs_Triangle_gam1_gam2 + i_fit_gam1_gam2;
                        for (int k = 0; k <= floor(i_fit_E_Orthog/2.0); k++) {
                            if (i_fit_E_Orthog - 2*k == i_fit_E) {
                                Coeffs_Cheby = Precompute_First_Kind_Chebyshev_Basis_Coeffs(i_fit_E_Orthog, k);
                                Coefficients_Fit_Monomials_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Monomial_Basis] += Coeffs_Cheby*Coefficients_Fit_Orthog_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Orthog_Basis];
                            }
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

void Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &N_Coeffs_Triangle_gam1_gam2) {
    long double Coeffs_Cheby;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_SH*N_Coeffs_Triangle_gam1_gam2;
    int N_Coeffs_Total_Without_SH = N_Coeffs_E*N_Coeffs_f*N_Coeffs_Triangle_gam1_gam2;
    int index_Monomial_Basis, index_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
        for (int i_fit_SH = 0; i_fit_SH < N_Coeffs_SH; i_fit_SH++) {
            for (int i_fit_gam1_gam2 = 0; i_fit_gam1_gam2 < N_Coeffs_Triangle_gam1_gam2; i_fit_gam1_gam2++) {
                
                for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
                    index_Monomial_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_Triangle_gam1_gam2 + i_fit_gam1_gam2;
                    
                    for (int i_fit_f_Orthog = 0; i_fit_f_Orthog < N_Coeffs_f; i_fit_f_Orthog++) {
                        index_Orthog_Basis = (i_fit_E*N_Coeffs_f + i_fit_f_Orthog)*N_Coeffs_Triangle_gam1_gam2 + i_fit_gam1_gam2;
                        for (int k = 0; k <= i_fit_f_Orthog; k++) {
                            if (2*i_fit_f_Orthog - 2*k == 2*i_fit_f) {
                                Coeffs_Cheby = Precompute_First_Kind_Chebyshev_Basis_Coeffs(2*i_fit_f_Orthog, k);
                                Coefficients_Fit_Monomials_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Monomial_Basis] += Coeffs_Cheby*Coefficients_Fit_Orthog_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Orthog_Basis];
                            }
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

long double Cartesian_Spherical_harmonics(const long double &x, const long double &y, const long double &z, const int &l, const int &m) {
    long double Coeffs_SH, Cartes_SH;
    long double Cartes_SH_j, Cartes_SH_k;
    long double x2, y2;
    
    // x ---> y
    // y ---> z
    // z ---> x
    
    if (m == 0) {
        for (int k = l; k >= 0; k--) {
            Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, 0, 0, k);
//             cout << "k = " << k << "   " << "l = " << l << "   " << "Coeffs_SH = " << Coeffs_SH << endl;
            if (k == l) {
                Cartes_SH = Coeffs_SH;
            } else {
                Cartes_SH = Coeffs_SH + x*Cartes_SH;
            }
        }
    } else if (m > 0) {
        x2 = x*x;
        y2 = y*y;
        for (int i = l; i >= m; i--) {
            for (int j = 0; j <= floor(m/2); j++) {
                for (int k = j; k >= 0; k--) {
                    Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
                    if (k == j) {
                        Cartes_SH_k = Coeffs_SH;
                    } else {
                        Cartes_SH_k = Coeffs_SH + x2*Cartes_SH_k;
                    }
                }
                if (j == 0) {
                    Cartes_SH_j = Cartes_SH_k;
                } else {
                    Cartes_SH_j = Cartes_SH_k + y2*Cartes_SH_j;
                }
            }
            
            if ((m % 2) != 0) { //then m is odd
                Cartes_SH_j *= y;
            }
            
            if (i == l) {
                Cartes_SH = Cartes_SH_j;
            } else {
                Cartes_SH = Cartes_SH_j + x*Cartes_SH;
            }
        }
    } else {
        x2 = x*x;
        y2 = y*y;
        for (int i = l; i >= fabs(m); i--) {
            for (int j = 0; j <= floor((fabs(m)-1)/2); j++) {
                // j should take values up to (m-1)/2
                for (int k = j; k >= 0; k--) {
                    Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
                    if (k == j) {
                        Cartes_SH_k = Coeffs_SH;
                    } else {
                        Cartes_SH_k = Coeffs_SH + x2*Cartes_SH_k;
                    }
                }
                
                if (j == 0) {
                    Cartes_SH_j = Cartes_SH_k;
                } else {
                    Cartes_SH_j = Cartes_SH_k + y2*Cartes_SH_j;
                }
            }
            if ((-m % 2) == 0) { //then m is even
                Cartes_SH_j *= y;
            }
            if (i == l) {
                Cartes_SH = Cartes_SH_j;
            } else {
                Cartes_SH = Cartes_SH_j + x*Cartes_SH;
            }
        }
        Cartes_SH *= z;
    }
    
//     cout << "l = " << l << "   " << "m = " << m << "   " << "x = " << x << "   " << "y = " << y << "   " << "z = " << z << "   " << "Cartes_SH = " << Cartes_SH << endl;
    
    return Cartes_SH;
}

// long double Cartesian_Spherical_harmonics(const long double &x, const long double &y, const long double &z, const int &l, const int &m) {
//     long double Coeffs_SH, Cartes_SH;
//     long double Cartes_SH_j, Cartes_SH_k;
//     long double x2, z2;
//     
//     if (m == 0) {
//         for (int k = l; k >= 0; k--) {
//             Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, 0, 0, k);
// //             cout << "k = " << k << "   " << "l = " << l << "   " << "Coeffs_SH = " << Coeffs_SH << endl;
//             if (k == l) {
//                 Cartes_SH = Coeffs_SH;
//             } else {
//                 Cartes_SH = Coeffs_SH + z*Cartes_SH;
//             }
//         }
//     } else if (m > 0) {
//         x2 = x*x;
//         z2 = z*z;
//         for (int i = l; i >= m; i--) {
//             for (int j = 0; j <= floor(m/2); j++) {
//                 for (int k = j; k >= 0; k--) {
//                     Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
//                     if (k == j) {
//                         Cartes_SH_k = Coeffs_SH;
//                     } else {
//                         Cartes_SH_k = Coeffs_SH + z2*Cartes_SH_k;
//                     }
//                 }
//                 if (j == 0) {
//                     Cartes_SH_j = Cartes_SH_k;
//                 } else {
//                     Cartes_SH_j = Cartes_SH_k + x2*Cartes_SH_j;
//                 }
//             }
//             
//             if ((m % 2) != 0) { //then m is odd
//                 Cartes_SH_j *= x;
//             }
//             
//             if (i == l) {
//                 Cartes_SH = Cartes_SH_j;
//             } else {
//                 Cartes_SH = Cartes_SH_j + z*Cartes_SH;
//             }
//         }
//     } else {
//         x2 = x*x;
//         z2 = z*z;
//         for (int i = l; i >= fabs(m); i--) {
//             for (int j = 0; j <= floor((fabs(m)-1)/2); j++) {
//                 // j should take values up to (m-1)/2
//                 for (int k = j; k >= 0; k--) {
//                     Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
//                     if (k == j) {
//                         Cartes_SH_k = Coeffs_SH;
//                     } else {
//                         Cartes_SH_k = Coeffs_SH + z2*Cartes_SH_k;
//                     }
//                 }
//                 
//                 if (j == 0) {
//                     Cartes_SH_j = Cartes_SH_k;
//                 } else {
//                     Cartes_SH_j = Cartes_SH_k + x2*Cartes_SH_j;
//                 }
//             }
//             if ((-m % 2) == 0) { //then m is even
//                 Cartes_SH_j *= x;
//             }
//             if (i == l) {
//                 Cartes_SH = Cartes_SH_j;
//             } else {
//                 Cartes_SH = Cartes_SH_j + z*Cartes_SH;
//             }
//         }
//         Cartes_SH *= y;
//     }
//     
// //     cout << "l = " << l << "   " << "m = " << m << "   " << "x = " << x << "   " << "y = " << y << "   " << "z = " << z << "   " << "Cartes_SH = " << Cartes_SH << endl;
//     
//     return Cartes_SH;
// }

long double Precompute_Proriol_Basis_Coeffs(const int &k, const int &l, const int &i, const int &j, const int &m, const int &n) {
    long double Proriol_Basis_Coeffs;
    long double Binomial_Coeffs;
    
    if (k == 0 && l == 0) {
        Proriol_Basis_Coeffs = 1.0;
    } else {
        Proriol_Basis_Coeffs = pow(2.0, k+m)*pow(-1.0, i+j-m+n);
    
        Binomial_Coeffs = binomialCoefficients(k, i);
        Proriol_Basis_Coeffs *= Binomial_Coeffs;
    
        Binomial_Coeffs = binomialCoefficients((long double)((k+i-1)/2.0), k);
        Proriol_Basis_Coeffs *= Binomial_Coeffs;
    
        Binomial_Coeffs = binomialCoefficients(l+2*k+1, l-j);
        Proriol_Basis_Coeffs *= Binomial_Coeffs;
    
        Binomial_Coeffs = binomialCoefficients(l, j);
        Proriol_Basis_Coeffs *= Binomial_Coeffs;
    
        Binomial_Coeffs = binomialCoefficients(i, m);
        Proriol_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients(k+j-m, n);
        Proriol_Basis_Coeffs *= Binomial_Coeffs;
    }
    
    // Add normalization constant
    Proriol_Basis_Coeffs *= sqrt(2.0*(2.0*k+1.0)*(k+l+1.0));
        
    return Proriol_Basis_Coeffs;
}

long double Proriol_Polynomials(const long double &zeta, const long double &eta, const int &k, const int &l) {
    long double Coeffs_Proriol, Proriol_poly;
    long double Proriol_poly_j, Proriol_poly_m, Proriol_poly_n;
    
    for (int i = k; i >= 0; i--) {
        for (int j = 0; j <= l; j++) {
            for (int m = i; m >= 0; m--) {
                for (int n = k+j-m; n >= 0; n--) {
                    Coeffs_Proriol = Precompute_Proriol_Basis_Coeffs(k, l, i, j, m, n);
                    if (n == k+j-m) {
                        Proriol_poly_n = Coeffs_Proriol;
                    } else {
                        Proriol_poly_n = Coeffs_Proriol + eta*Proriol_poly_n;
                    }
                }
                if (m == i) {
                    Proriol_poly_m = Proriol_poly_n;
                } else {
                    Proriol_poly_m = Proriol_poly_n + zeta*Proriol_poly_m;
                }
            }
            
            if (j == 0) {
                Proriol_poly_j = Proriol_poly_m;
            } else {
                Proriol_poly_j = Proriol_poly_m + eta*Proriol_poly_j;
            }
        }
        
        if (i == k) {
            Proriol_poly = Proriol_poly_j;
        } else {
            Proriol_poly = Proriol_poly + Proriol_poly_j;
        }
    }
        
    return Proriol_poly;
}

void Spherical_Harmonics_to_Monomial_Basis(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &N_Coeffs_Triangle_gam1_gam2, const int &Order_SH, const int *Array_l_SH, const int *Array_m_SH) {
    long double Coeffs_SH;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_SH*N_Coeffs_Triangle_gam1_gam2;
    int N_Coeffs_Total_Without_SH = N_Coeffs_E*N_Coeffs_f*N_Coeffs_Triangle_gam1_gam2;
    int index_Basis;
    
    int index_SH_Monomial_Basis, index_SH_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
        for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
            for (int i_fit_gam1_gam2 = 0; i_fit_gam1_gam2 < N_Coeffs_Triangle_gam1_gam2; i_fit_gam1_gam2++) {
                index_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_Triangle_gam1_gam2 + i_fit_gam1_gam2;
                
                index_SH_Monomial_Basis = 0;
                for (int p = 0; p <= Order_SH; p+=2) {
                    for (int q = 0; q <= Order_SH - p; q+=2) {
                        index_SH_Orthog_Basis = 0;
                        for (int l = 0; l <= Order_SH; l+=2) {
                            for (int m = 0; m <= l; m+=2) {
                                if (m == 0) {
                                    for (int k = 0; k <= l; k++) {
                                        if (p == k && q == 0) {
                                            Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, 0, 0, k);
                                            Coefficients_Fit_Monomials_Basis[index_SH_Monomial_Basis*N_Coeffs_Total_Without_SH+index_Basis] += Coeffs_SH*Coefficients_Fit_Orthog_Basis[index_SH_Orthog_Basis*N_Coeffs_Total_Without_SH+index_Basis];  
                                        }
                                    }
                                } else {
                                    for (int i = m; i <= l; i++) {
                                        for (int j = 0; j <= floor(m/2.0); j++) {
                                            for (int k = 0; k <= j; k++) {
                                                if (p == i-m+2*k && q == m - 2*j) {
                                                    Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
                                                    Coefficients_Fit_Monomials_Basis[index_SH_Monomial_Basis*N_Coeffs_Total_Without_SH+index_Basis] += Coeffs_SH*Coefficients_Fit_Orthog_Basis[index_SH_Orthog_Basis*N_Coeffs_Total_Without_SH+index_Basis];  
                                                } 
//                                                 else if ((i-m+2*k) % 2 != 0 || (m-2*j) % 2 != 0) {
//                                                     Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
//                                                     cout << "l = " << l << "   " << "m = " << m << "   " << "(i-m+2*k) = " << (i-m+2*k) << "   " << "(m-2*j) = " << (m-2*j) << "   " << "Coeffs_SH = " << Coeffs_SH << endl;
//                                                 }
                                            }
                                        }
                                    }
                                }
                                index_SH_Orthog_Basis++;
                            }
                        }
                        index_SH_Monomial_Basis++;
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

void Proriol_to_Monomial_Basis(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &N_Coeffs_Triangle_gam1_gam2, const int &N_Coeffs_gam1) {
    long double Coeffs_Proriol;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_SH*N_Coeffs_Triangle_gam1_gam2;
    int N_Coeffs_Total_Without_SH = N_Coeffs_E*N_Coeffs_f*N_Coeffs_Triangle_gam1_gam2;
    int index_Monomial_Basis, index_Orthog_Basis;
    
    int index_Triangle_Monomial_Basis, index_Triangle_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
        for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
            for (int i_fit_SH = 0; i_fit_SH < N_Coeffs_SH; i_fit_SH++) {
                
                index_Triangle_Monomial_Basis = 0;
                for (int p = 0; p < N_Coeffs_gam1; p++) {
                    for (int q = 0; q < N_Coeffs_gam1-p; q++) {
                        index_Monomial_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_Triangle_gam1_gam2 + index_Triangle_Monomial_Basis;
                        index_Triangle_Monomial_Basis++;
                        
                        index_Triangle_Orthog_Basis = 0;
                        for (int k = 0; k < N_Coeffs_gam1; k++) {
                            for (int l = 0; l < N_Coeffs_gam1-k; l++) {
                                index_Orthog_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_Triangle_gam1_gam2 + index_Triangle_Orthog_Basis;
                                
                                index_Triangle_Orthog_Basis++;
                                
                                for (int i = 0; i <= k; i++) {
                                    for (int j = 0; j <= l; j++) {
                                        for (int m = 0; m <= i; m++) {
                                            for (int n = 0; n <= k+j-m; n++) {
                                                
                                                if (m == p && l-j+n == q) {
                                                    Coeffs_Proriol = Precompute_Proriol_Basis_Coeffs(k, l, i, j, m, n);
                                                    Coefficients_Fit_Monomials_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Monomial_Basis] += Coeffs_Proriol*Coefficients_Fit_Orthog_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Orthog_Basis];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

int SH_Linear_Index(const int &Order_SH, const int &degree_SH) {
   int index_SH;
   index_SH = Order_SH*Order_SH + Order_SH + degree_SH;
   return index_SH;
}

int SH_Linear_Index_Symmetric(const int &Order_SH, const int &degree_SH) {
   int index_SH;
   index_SH = Order_SH*(2*Order_SH+1) + degree_SH;
//    index_SH = Order_SH*Order_SH + degree_SH;
   return index_SH;
}

void Triangle_to_Square_Mapping(long double &zeta, long double &eta, const long double &v_x, const long double &v_y) {
    zeta = 1.0 + (v_x - v_y) - sqrt(pow(v_x - v_y, 2) + 4.0*(1.0 - v_x - v_y));
    eta = 1.0 - (v_x - v_y) - sqrt(pow(v_x - v_y, 2) + 4.0*(1.0 - v_x - v_y));
//     zeta = 2.0*v_x/(1.0 - v_y) - 1.0;
//     eta = 2.0*v_y - 1.0;
}

long double VanderMonde_Matrix_1_var(const long double &var_1, const int &Order_poly_1) {
    long double Van_Entry;
    long double poly_1;
    poly_1 = Chebyshev_Polynomial_Basis(var_1, Order_poly_1);
    
    Van_Entry = poly_1;
    return Van_Entry;
}

long double VanderMonde_Matrix_2_vars(const long double &var_1, const long double &var_2, const int &Order_poly_1, const int &Order_poly_2) {
    long double Van_Entry;
    long double poly_1, poly_2;
    poly_1 = Chebyshev_Polynomial_Basis(var_1, Order_poly_1);
    poly_2 = Chebyshev_Polynomial_Basis(var_2, Order_poly_2);
    
    Van_Entry = poly_1*poly_2;
    return Van_Entry;
}

long double VanderMonde_Matrix_3_vars(const long double &var_1, const long double &var_2, const long double &var_3, const int &Order_poly_1, const int &Order_poly_2, const int &Order_poly_3) {
    long double Van_Entry;
    long double poly_1, poly_2, poly_3;
    poly_1 = Chebyshev_Polynomial_Basis(var_1, Order_poly_1);
    poly_2 = Chebyshev_Polynomial_Basis(var_2, Order_poly_2);
    poly_3 = Chebyshev_Polynomial_Basis(var_3, Order_poly_3);
    
    Van_Entry = poly_1*poly_2*poly_3;
    return Van_Entry;
}

long double VanderMonde_Matrix_4_vars(const long double &var_1, const long double &var_2, const long double &var_3, const long double &var_4, const int &Order_poly_1, const int &Order_poly_2, const int &Order_poly_3, const int &Order_poly_4) {
    long double Van_Entry;
    long double poly_1, poly_2, poly_3, poly_4;
    poly_1 = Chebyshev_Polynomial_Basis(var_1, Order_poly_1);
    poly_2 = Chebyshev_Polynomial_Basis(var_2, Order_poly_2);
    poly_3 = Chebyshev_Polynomial_Basis(var_3, Order_poly_3);
    poly_4 = Chebyshev_Polynomial_Basis(var_4, Order_poly_4);
    
    Van_Entry = poly_1*poly_2*poly_3*poly_4;
    return Van_Entry;
}

long double VanderMonde_Vector_N_vars(const int &Index_Entry, const int &Index_Point) {
    long double Van_Vector;
    
    if (Index_Entry == Index_Point) {
        Van_Vector = 1.0;
    } else {
        Van_Vector = 0.0;
    }
    return Van_Vector;
}

void Solve_A_x_b(const long double *VanderMonde_Matrix, long double *Coeff_Vand_Matrix, const long double *VanderMonde_Vector, const int &n, const int &index_Cheby) {
    // First perform Cholesky decomposition of A
//                 cout << "ratio_EBBBBBBBBBBBBBBB = " << ratio_E << endl;
    long double *Temp_Mat_L, *Temp_Mat_U, *Temp_Vec;
    Temp_Mat_L = new long double[n*n];
    Temp_Mat_U = new long double[n*n];
    Temp_Vec = new long double[n];
    
    LUdecomposition(VanderMonde_Matrix, Temp_Mat_L, Temp_Mat_U, n);
    
    LU_Solve_Forward(n, Temp_Vec, Temp_Mat_L, VanderMonde_Vector);
    
    LU_Solve_Back(n, Temp_Vec, Temp_Mat_U, Temp_Vec);
    
    for (int i = 0; i < n; i++) {
        Coeff_Vand_Matrix[i*n + index_Cheby] = Temp_Vec[i];
    }
    
    delete[] Temp_Mat_L;
    delete[] Temp_Mat_U;
    delete[] Temp_Vec;
}

void Check_A_x_b(const long double *VanderMonde_Matrix, const long double *Coeff_Vand_Matrix, const int &n) {
    // First perform Cholesky decomposition of A
    long double temp_val;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            temp_val = 0.0;
            for (int k = 0; k < n; k++) {
                temp_val += VanderMonde_Matrix[i*n + k] * Coeff_Vand_Matrix[k*n + j];
            }
            if (i == j) {
                if (fabs(1.0 - temp_val) > 1.0e-6) {
                    cout << "Issue with Vandermonde System:" << endl;
                    cout << "i = " << i << "  " << "j = " << j << "  " << "temp_val = " << temp_val << endl;
                    exit(0);
                }
            } else {
                if (fabs(temp_val) > 1.0e-6) {
                    cout << "Issue with Vandermonde System:" << endl;
                    cout << "i = " << i << "  " << "j = " << j << "  " << "temp_val = " << temp_val << endl;
                    exit(0);
                }
            }
        }
    }
}

void LUdecomposition(const long double *A, long double *L, long double *U, const int &n) {
   for (int i = 0; i < n; i++) {
       for (int j = 0; j < n; j++) {
           if (j < i) {
               L[j*n+i] = 0;
           } else {
               L[j*n+i] = A[j*n+i];
               for (int k = 0; k < i; k++) {
                   L[j*n+i] = L[j*n+i] - L[j*n+k] * U[k*n+i];
            }
         }
      }
      for (int j = 0; j < n; j++) {
          if (j < i) {
              U[i*n+j] = 0;
          } else if (j == i) {
              U[i*n+j] = 1;
          } else {
              U[i*n+j] = A[i*n+j] / L[i*n+i];
              for (int k = 0; k < i; k++) {
                  U[i*n+j] = U[i*n+j] - ((L[i*n+k] * U[k*n+j]) / L[i*n+i]);
            }
         }
      }
   }
}

void LU_Solve_Back(const int &n, long double *x, const long double *A, const long double *b)
{
    // Back solve A x = b
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < n; j++) {
            x[i] -= A[i*n+j] * x[j];
        }
//         if (fabs(A[i*n+i]) < 1.0e-10) {
//             x[i] = 0.0;
//         } else {
            x[i] /= A[i*n+i];
//         }
    }
}

void LU_Solve_Forward(const int &n, long double *x, const long double *A, const long double *b)
{
    // Forward solve A x = b
    for (int i = 0; i < n; i++) {
        x[i] = b[i];
        for (int j = 0; j < i; j++) {
            x[i] -= A[i*n+j] * x[j];
        }
        
//         if (fabs(A[i*n+i]) < 1.0e-10) {
//             x[i] = 0.0;
//         } else {
            x[i] /= A[i*n+i];
//         }
    }
}

long double Lebedev_Quadrature_Matrix_SH_Temp(const long double *Matrix, const int &Order_SH, const int &degree_SH, const int &Npts_Phi, const int &Npts_Theta, const long double *x_SH, const long double *y_SH, const long double *z_SH) {
    long double *x, *y, *z, *w;
    long double integral_approx, Y_l_m;
    
    x = new long double[Npts_Phi*Npts_Theta];
    y = new long double[Npts_Phi*Npts_Theta];
    z = new long double[Npts_Phi*Npts_Theta];
    w = new long double[Npts_Phi*Npts_Theta];
    
    setup_spherical_harmonics_data (Npts_Phi, Npts_Theta, x, y, z, w);
    
    integral_approx = 0.0;
    
    for ( int m = 0; m < Npts_Phi*Npts_Theta; m++ ) {
        Y_l_m = Cartesian_Spherical_harmonics(x[m], y[m], z[m], Order_SH, degree_SH);
        integral_approx = integral_approx + w[m] * Matrix[m] * Y_l_m;
        
        if (x_SH[m] != x[m] || y_SH[m] != y[m] || z_SH[m] != z[m]) {
            cout << "x = " << x[m] << "   " << "x_SH = " << x_SH[m] << endl;
            cout << "y = " << y[m] << "   " << "y_SH = " << y_SH[m] << endl;
            cout << "z = " << z[m] << "   " << "z_SH = " << z_SH[m] << endl;
            exit(0);
        }
        
//         cout << "x = " << x[m] << "   " << "y = " << y[m] << "   " << "z = " << z[m] << "   " << "l = " << Order_SH << "   " << "m = " << degree_SH << "   " << "Y_l_m = " << Y_l_m << "   " << "Matrix = " << Matrix[m] << "   " << "integral_approx = " << integral_approx << endl;
    }
    
//     integral_approx *= 4.0;
    
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] w;
    
    return integral_approx;
}

long double Lebedev_Quadrature_Matrix_SH(const long double *Matrix, const int &Order_SH, const int &degree_SH, const int &Npts_Phi, const int &Npts_Theta) {
    long double *x, *y, *z, *w;
    long double integral_approx, Y_l_m;
    
    x = new long double[Npts_Phi*Npts_Theta];
    y = new long double[Npts_Phi*Npts_Theta];
    z = new long double[Npts_Phi*Npts_Theta];
    w = new long double[Npts_Phi*Npts_Theta];
    
    setup_spherical_harmonics_data (Npts_Phi, Npts_Theta, x, y, z, w);
    
    integral_approx = 0.0;
    
    for ( int m = 0; m < Npts_Phi*Npts_Theta; m++ ) {
        Y_l_m = Cartesian_Spherical_harmonics(x[m], y[m], z[m], Order_SH, degree_SH);
        integral_approx = integral_approx + w[m] * Matrix[m] * Y_l_m;
    }
    
    
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] w;
    
    return integral_approx;
}

long double Lebedev_Quadrature_Orthog_SH(const int &order_SH1, const int &order_SH2, const int &m1, const int &m2, const int &Npts_Phi, const int &Npts_Theta) {
    long double *x, *y, *z, *w;
    long double integral_approx, Y_l_m_1, Y_l_m_2;
    
    x = new long double[Npts_Phi*Npts_Theta];
    y = new long double[Npts_Phi*Npts_Theta];
    z = new long double[Npts_Phi*Npts_Theta];
    w = new long double[Npts_Phi*Npts_Theta];
    
    setup_spherical_harmonics_data (Npts_Phi, Npts_Theta, x, y, z, w);
    
    integral_approx = 0.0;
    
    for ( int m = 0; m < Npts_Phi*Npts_Theta; m++ ) {
        Y_l_m_1 = Cartesian_Spherical_harmonics(x[m], y[m], z[m], order_SH1, m1);
        Y_l_m_2 = Cartesian_Spherical_harmonics(x[m], y[m], z[m], order_SH2, m2);
//         cout << "Y_l_m_1 = " << Y_l_m_1 << "   " << "Y_l_m_2 = " << Y_l_m_2 << endl;
        integral_approx = integral_approx + w[m] * Y_l_m_1 * Y_l_m_2;
    }
    
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] w;
    
    return integral_approx;
}

long double Triangle_Quadrature_Orthog_Proriol(const int &i_gam1, const int &i_gam2, const int &i_prime_gam1, const int &i_prime_gam2, const int &N_points_Leg) {
    long double *mu_quad, *w_Leg;
    long double integral_approx, Proriol_poly_1, Proriol_poly_2;
    long double zeta, eta;
    
    mu_quad = new long double[N_points_Leg];
    w_Leg = new long double[N_points_Leg];
    
    legendre_set ( N_points_Leg, mu_quad, w_Leg );
    
    integral_approx = 0.0;
    
    // Duffy transform for mapping the square to the standard triangle
    
    for ( int m = 0; m < N_points_Leg; m++ ) {
        for ( int n = 0; n < N_points_Leg; n++ ) {
            zeta = (1.0 + mu_quad[m])*(1.0 - mu_quad[n])/4.0;
            eta = (1.0 + mu_quad[n])/2.0;
            
            // Duffy Transform ==> d zeta = d zeta_p (1.0 - eta_p)/4.0 - d eta_p (1.0 + zeta_p)/4.0;
            //                 and d eta = d eta_p/4.0;
            //                 ==> d zeta d eta = (1.0 - eta_p)/8.0 d zeta_p d eta_p
            
            Proriol_poly_1 = Proriol_Polynomials(zeta, eta, i_gam1, i_gam2);
            Proriol_poly_2 = Proriol_Polynomials(zeta, eta, i_prime_gam1, i_prime_gam2);
            //         cout << "Proriol_poly_1 = " << Proriol_poly_1 << "   " << "Proriol_poly_2 = " << Proriol_poly_2 << endl;
            integral_approx = integral_approx + w_Leg[m] * w_Leg[n] * Proriol_poly_1 * Proriol_poly_2 * (1.0 - mu_quad[n])/8.0;
        }
    }
    
    delete[] mu_quad;
    delete[] w_Leg;
    
    return integral_approx;
}

long double Test_Spherical_Harmonics_Even_Odd(const long double *Vals, const int &index_Cheby, const int &N_points_Cheby, const long double &x, const long double &y, const long double &z, const int &Order_SH) {
    long double f_SH;
    f_SH = 0.0;
    int index = 0;
    
    for (int i_SH = 0; i_SH <= Order_SH; i_SH++) {
        for (int l_SH = -i_SH; l_SH <= i_SH; l_SH++) {
            f_SH += Vals[index*N_points_Cheby + index_Cheby]*Cartesian_Spherical_harmonics(x,y,z,i_SH,l_SH);
            index++;
        }
    }
    return f_SH;
}

long double Test_Spherical_Harmonics(const long double *Vals, const int &index_Cheby, const int &N_points_Cheby, const long double &x, const long double &y, const long double &z, const int &Order_SH) {
    long double f_SH;
    f_SH = 0.0;
    int index = 0;
    
    for (int i_SH = 0; i_SH <= Order_SH; i_SH+=2) {
        for (int l_SH = 0; l_SH <= i_SH; l_SH+=2) {
            f_SH += Vals[index*N_points_Cheby + index_Cheby]*Cartesian_Spherical_harmonics(x,y,z,i_SH,l_SH);
            index++;
        }
    }
    return f_SH;
}

void Orthogonality_Test(const int &Quad_Rule) {
    long double delta_ij;
    for (int i_SH1 = 0; i_SH1 <= 10; i_SH1+=1) {
        for (int i_SH2 = 0; i_SH2 <= 10; i_SH2+=1) {
            for (int m_SH1 = -i_SH1; m_SH1 <= i_SH1; m_SH1+=1) {
                for (int m_SH2 = -i_SH2; m_SH2 <= i_SH2; m_SH2+=1) {
                    
//                     if (m_SH1 >= 0 && m_SH2 >= 0) {
                        delta_ij = Lebedev_Quadrature_Orthog_SH(i_SH1, i_SH2, m_SH1, m_SH2, Quad_Rule, Quad_Rule);
                        
//                         cout << "i_SH1, i_SH2 = " << i_SH1 << "    " << i_SH2 << "      " << "m_SH1, m_SH2 = " << m_SH1 << "    " << m_SH2 << "     " << "delta_ij = " << delta_ij << endl;
                        
                        if ((i_SH1 == i_SH2) && (m_SH1 == m_SH2)) {
                            if (fabs(delta_ij - 1.0) > 1.0e-8) {
                                cout << "i_SH1, i_SH2 = " << i_SH1 << "    " << i_SH2 << "      " << "m_SH1, m_SH2 = " << m_SH1 << "    " << m_SH2 << "     " << "delta_ij = " << delta_ij << endl;
                                exit(0);
                            }
                        } else {
                            if (fabs(delta_ij) > 1.0e-8) {
                                cout << "i_SH1, i_SH2 = " << i_SH1 << "    " << i_SH2 << "      " << "m_SH1, m_SH2 = " << m_SH1 << "    " << m_SH2 << "     " << "delta_ij = " << delta_ij << endl;
                                exit(0);
                            }
                        }
//                     }
                }
            }
        }
    }
}

void Test_Proriol_Polynomials_Orthogonality(const int &Quad_Rule) {
    long double delta_ij;
    for (int i_gam1 = 0; i_gam1 < 10; i_gam1++) {
        for (int i_gam2 = 0; i_gam2 < 10 - i_gam1; i_gam2++) {
            
            for (int i_prime_gam1 = 0; i_prime_gam1 < 10; i_prime_gam1++) {
                for (int i_prime_gam2 = 0; i_prime_gam2 < 10 - i_prime_gam1; i_prime_gam2++) {
                    
                    delta_ij = Triangle_Quadrature_Orthog_Proriol(i_gam1, i_gam2, i_prime_gam1, i_prime_gam2, Quad_Rule);
                            
                    if ((i_gam1 == i_prime_gam1) && (i_gam2 == i_prime_gam2)) {
                        if (fabs(delta_ij - 1.0) > 1.0e-8) {
                            cout << "i_gam1, i_prime_gam1 = " << i_gam1 << "    " << i_prime_gam1 << "      " << "i_gam2, i_prime_gam2 = " << i_gam2 << "    " << i_prime_gam2 << "     " << "delta_ij = " << delta_ij << endl;
                            exit(0);
                        }
                    } else {
                        if (fabs(delta_ij) > 1.0e-8) {
                            cout << "i_gam1, i_prime_gam1 = " << i_gam1 << "    " << i_prime_gam1 << "      " << "i_gam2, i_prime_gam2 = " << i_gam2 << "    " << i_prime_gam2 << "     " << "delta_ij = " << delta_ij << endl;
                            exit(0);
                        }
                    }
                }
            }
        }
    }
}


long double Cartesian_to_spherical_Coordinates_Conversion_Jacobian(const int &i, const int &j, const long double &r, const long double &mu, const long double &cos_phi, const long double &sin_phi) {
    long double mat = 0.0;
    long double sin_theta = sqrt(1.0 - pow(mu, 2));
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = mu;
                    break;
                case 1:
                    mat = r;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 0:
                    mat = sin_theta*cos_phi;
                    break;
                case 1:
                    mat = r*cos_phi*mu/sin_theta;
                    break;
                case 2:
                    mat = -r*sin_theta*sin_phi;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 0:
                    mat = sin_theta*sin_phi;
                    break;
                case 1:
                    mat = r*sin_phi*mu/sin_theta;
                    break;
                case 2:
                    mat = r*sin_theta*cos_phi;
                    break;
            };
            break;
    };
    
    return mat;
}

long double Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(const int &i, const int &j, const long double &r, const long double &cos_theta, const long double &cos_phi, const long double &sin_phi) {
    long double mat = 0.0;
    long double sin_theta = sqrt(1.0 - pow(cos_theta, 2));
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = cos_theta;
                    break;
                case 1:
                    mat = -r*cos_theta;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 0:
                    mat = sin_theta*cos_phi;
                    break;
                case 1:
                    mat = -r*sin_theta*cos_phi;
                    break;
                case 2:
                    mat = -r*sin_theta*sin_phi;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 0:
                    mat = sin_theta*sin_phi;
                    break;
                case 1:
                    mat = -r*sin_theta*sin_phi;
                    break;
                case 2:
                    mat = r*sin_theta*cos_phi;
                    break;
            };
            break;
    };
    
    cout << "Update this !!!!!!!!!!!!" << endl;
    return mat;
}

long double Cartesian_to_spherical_Coordinates_Jacobian(const long double &norm_f, const long double &x_val, const long double &y_val, const long double &z_val, const long double &dIn_Plus_dN1_1, const long double &dIn_Plus_dN1_2, const long double &dIn_Plus_dN1_3, const int &VAR_NUM) {
    long double mu, sin_theta, cos_phi, sin_phi;
    long double dIn_Plus_dmu;
    long double dN1_1_dmu, dN1_2_dmu, dN1_3_dmu;
    
    mu = x_val;
    sin_theta = sqrt(1.0 - pow(mu,2));
    
    cos_phi = y_val/sin_theta;
    sin_phi = z_val/sin_theta;
    
    if (sin_theta < 1.0e-8) {
        cos_phi = 1.0;
        sin_phi = 0.0;
    }
    
    dN1_1_dmu = Cartesian_to_spherical_Coordinates_Conversion_Jacobian(0, VAR_NUM, norm_f, mu, cos_phi, sin_phi);
    dN1_2_dmu = Cartesian_to_spherical_Coordinates_Conversion_Jacobian(1, VAR_NUM, norm_f, mu, cos_phi, sin_phi);
    dN1_3_dmu = Cartesian_to_spherical_Coordinates_Conversion_Jacobian(2, VAR_NUM, norm_f, mu, cos_phi, sin_phi);
    
    dIn_Plus_dmu = dIn_Plus_dN1_1*dN1_1_dmu + dIn_Plus_dN1_2*dN1_2_dmu + dIn_Plus_dN1_3*dN1_3_dmu;
    
    return dIn_Plus_dmu;
}

long double Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(const long double &norm_f, const long double &x_val, const long double &y_val, const long double &z_val, const long double &dIn_Plus_dN1_1, const long double &dIn_Plus_dN1_2, const long double &dIn_Plus_dN1_3, const int &VAR_NUM) {
    long double cos_theta, sin_theta, cos_phi, sin_phi;
    long double dIn_Plus_dmu;
    long double dN1_1_dtheta, dN1_2_dtheta, dN1_3_dtheta;
    
    cos_theta = x_val;
    sin_theta = sqrt(1.0 - pow(cos_theta,2));
    
    cos_phi = y_val/sin_theta;
    sin_phi = z_val/sin_theta;
    
    if (sin_theta < 1.0e-8) {
        cos_phi = 1.0;
        sin_phi = 0.0;
    }
    
    cout << "Update this !!!!!!!!!!!!" << endl;
    
    dN1_1_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(0, VAR_NUM, norm_f, cos_theta, cos_phi, sin_phi);
    dN1_2_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(1, VAR_NUM, norm_f, cos_theta, cos_phi, sin_phi);
    dN1_3_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(2, VAR_NUM, norm_f, cos_theta, cos_phi, sin_phi);
    
    dIn_Plus_dmu = dIn_Plus_dN1_1*dN1_1_dtheta + dIn_Plus_dN1_2*dN1_2_dtheta + dIn_Plus_dN1_3*dN1_3_dtheta;
    
    return dIn_Plus_dmu;
}
