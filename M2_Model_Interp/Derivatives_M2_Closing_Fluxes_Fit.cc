#ifndef _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED
#include "N3_Non_Gray_M2_3D_RT_Cheby.h"
#endif // _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED

//***************************************************************************************************
//
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Recompute_I0_Mapping(const long double &ratio_E, const long double &L_N3_ijk_Actual, const long double &L_N3_ijk_Fit) {
    long double ratio_E_temp;
    
    if (fabs(ratio_E) == 1.0) {
        ratio_E_temp = ratio_E;
    } else {
        // Compute the energy density associated with the length scale associated
        // with the uniform distribution
        ratio_E_temp = Inverse_Mobius_Transformation(ratio_E, L_N3_ijk_Actual);
        
        // Compute the value of the mapping of I0 associated with the approximated
        // optimal length scale
        ratio_E_temp = Mobius_Transformation(ratio_E_temp, L_N3_ijk_Fit);
    }
    
    return ratio_E_temp;
}

//***************************************************************************************************
// This routine computes the interpolative-based approximation of the optimal length scale, L_N3_ijk, 
// for the interpolative-based approximation of the third-order closing flux, N3_ijk
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_Length_Scale_N3(const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double Length_Scale;
    
    switch (closure_type) {
        case MOMENT_CLOSURE_PROJECTION:
            Length_Scale = 1.0;
            break;
        case MOMENT_CLOSURE_M2:
            Length_Scale = Evaluate_Length_Scale_N3_ijk(N1_1, N1_2, N1_3, gam1, gam2);
            break;
        default:
            cout << "Incorrect value for Closure Type" << endl;
            exit(0);
            break;
    };
    
    return Length_Scale;
}

//***************************************************************************************************
// This routine computes the interpolative-based approximation of the optimal length scale, L_N3_ijk, 
// for the interpolative-based approximation of the third-order closing flux, N3_ijk
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_Length_Scale_N3_ijk(const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double L, L_m, L_f, L_gam1, L_gam2;
    long double norm_f_2, norm_f;
    int index_Coeffs;
    int N_Points_Gray_Cheby_without_Lebed = N_Points_f*N_Points_Triangle_gam1_gam2;
    long double x, y, z, x2, y2, z2;
    int Order_SH_L_N3_ijk = 0; //Order_SH; //
    int N_Coeffs_SH_L_N3_ijk = 1; //N_Coeffs_SH; //
    
    index_Coeffs = N_Points_Gray_Cheby_without_Lebed*N_Coeffs_SH_L_N3_ijk - 1;
    
    norm_f_2 = pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2);
    norm_f = sqrt(norm_f_2);
    
    x = N1_1/norm_f;
    y = N1_2/norm_f;
    z = N1_3/norm_f;
    
    if (norm_f < 1.0e-8) {
        x = 1.0;
        y = 0.0;
        z = 0.0;
    }
    
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    for (int l = Order_SH_L_N3_ijk; l >= 0; l-=2) {
        for (int m = Order_SH_L_N3_ijk - l; m >= 0; m-=2) {
            for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
                for (int i_fit_gam1 = N_Points_gam1 - 1; i_fit_gam1 >= 0; i_fit_gam1--) {
                    for (int i_fit_gam2 = N_Points_gam1 - i_fit_gam1 - 1; i_fit_gam2 >= 0; i_fit_gam2--) {
                        if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
                            L_gam2 = Coefficient_Mobius_Scale_N3_ijk[index_Coeffs];
                        } else {
                            L_gam2 = Coefficient_Mobius_Scale_N3_ijk[index_Coeffs] + L_gam2*gam2;
                        }
                        index_Coeffs--;
                    }
                    if (i_fit_gam1 == N_Points_gam1 - 1) {
                        L_gam1 = L_gam2;
                    } else {
                        L_gam1 = L_gam2 + L_gam1*gam1;
                    }
                }
                if (i_fit_f == N_Points_f - 1) {
                    L_f = L_gam1;
                } else {
                    L_f = L_gam1 + L_f*norm_f_2;
                } 
            }
            if (m == Order_SH_L_N3_ijk - l) {
                L_m = L_f;
            } else {
                // L_m = L_f + L_m*x2;
                L_m = L_f + L_m*y2;
            }
        }
        if (l == Order_SH_L_N3_ijk) {
            L = L_m;
        } else {
            // L = L_m + L*z2;
            L = L_m + L*x2;
        }
    }
    
    // cout << "Update this for INDEX_y_SH_2 !!!!!!!!!!!!!" << endl;
    // exit(0);
    
    // L = Coefficient_Mobius_Scale_N3_ijk[0];
    
    L = exp(L);
    
    return L;
}

//***************************************************************************************************
// This routine computes derivatives of the interpolative-based approximation of the optimal length 
// scale, L_N3_111, with respect to the conserved solutions
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_diff_Length_Scale_N3_ijk(const int &index_var_fit, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double L, L_m, L_f, L_gam1, L_gam2;
    long double norm_f_2, norm_f;
    int index_Coeffs;
    int N_Points_Gray_Cheby_without_Lebed = N_Points_f*N_Points_Triangle_gam1_gam2;
    long double x, y, z, x2, y2, z2;
    long double Length_Scale;
    int Order_SH_L_N3_ijk = 0;
    int N_Coeffs_SH_L_N3_ijk = 1;
    
    index_Coeffs = N_Points_Gray_Cheby_without_Lebed*N_Coeffs_SH_L_N3_ijk - 1;
    
    norm_f_2 = pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2);
    norm_f = sqrt(norm_f_2);
    
    x = N1_1/norm_f;
    y = N1_2/norm_f;
    z = N1_3/norm_f;
    
    if (norm_f < 1.0e-8) {
        x = 1.0;
        y = 0.0;
        z = 0.0;
    }
    
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    Length_Scale = Evaluate_Length_Scale_N3_ijk(N1_1, N1_2, N1_3, gam1, gam2);
    
    for (int l = Order_SH_L_N3_ijk; l >= 0; l-=2) {
        for (int m = Order_SH_L_N3_ijk - l; m >= 0; m-=2) {
            for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
                for (int i_fit_gam1 = N_Points_gam1 - 1; i_fit_gam1 >= 0; i_fit_gam1--) {
                    for (int i_fit_gam2 = N_Points_gam1 - i_fit_gam1 - 1; i_fit_gam2 >= 0; i_fit_gam2--) {
                        if (index_var_fit == INDEX_GAM2) {
                            if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
                                L_gam2 = i_fit_gam2*Coefficient_Mobius_Scale_N3_ijk[index_Coeffs];
                            } else {
                                L_gam2 = i_fit_gam2*Coefficient_Mobius_Scale_N3_ijk[index_Coeffs] + L_gam2*gam2;
                            }
                        } else {
                            if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
                                L_gam2 = Coefficient_Mobius_Scale_N3_ijk[index_Coeffs];
                            } else {
                                L_gam2 = Coefficient_Mobius_Scale_N3_ijk[index_Coeffs] + L_gam2*gam2;
                            }
                        }
                        index_Coeffs--;
                    }
                    if (index_var_fit == INDEX_GAM1) {
                        if (i_fit_gam1 == N_Points_gam1 - 1) {
                            L_gam1 = i_fit_gam1*L_gam2;
                        } else {
                            L_gam1 = i_fit_gam1*L_gam2 + L_gam1*gam1;
                        }
                    } else {
                        if (i_fit_gam1 == N_Points_gam1 - 1) {
                            L_gam1 = L_gam2;
                        } else {
                            L_gam1 = L_gam2 + L_gam1*gam1;
                        }
                    }
                }
                if (index_var_fit == INDEX_NORM_f_2) {
                    if (i_fit_f == N_Points_f - 1) {
                        L_f = i_fit_f*L_gam1;
                    } else {
                        L_f = i_fit_f*L_gam1 + L_f*norm_f_2;
                    }
                } else {
                    if (i_fit_f == N_Points_f - 1) {
                        L_f = L_gam1;
                    } else {
                        L_f = L_gam1 + L_f*norm_f_2;
                    }
                } 
            }
            if (index_var_fit == INDEX_x_SH_2) {
                if (m == Order_SH_L_N3_ijk - l) {
                    L_m = m*L_f;
                } else {
                    // L_m = m*L_f + L_m*x2;
                    L_m = m*L_f + L_m*y2;
                }
            } else {
                if (m == Order_SH_L_N3_ijk - l) {
                    L_m = L_f;
                } else {
                    // L_m = L_f + L_m*x2;
                    L_m = L_f + L_m*y2;
                }
            }
        }
        if (index_var_fit == INDEX_z_SH_2) {   
            if (l == Order_SH_L_N3_ijk) {
                L = l*L_m;
            } else {
                // L = l*L_m + L*z2;
                L = l*L_m + L*x2;
            }
        } else {   
            if (l == Order_SH_L_N3_ijk) {
                L = L_m;
            } else {
                // L = L_m + L*z2;
                L = L_m + L*x2;
            }
        }
    }
    
    // cout << "Update this for INDEX_y_SH_2 !!!!!!!!!!!!!" << endl;
    // exit(0);
    
    L *= Length_Scale;
    
    return L;
}

//***************************************************************************************************
// This routine computes the exponential mapping, E_I0_eta, of the radiative energy density, used for
// the interpolative-based approximation of the third-order closing fluxes, N3_ijk
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Exponential_Mapping(const long double &I0_star, const long double &L_I0_star_N3) {
    long double map_I0_star;
    
    map_I0_star = 1.0 - 2.0*exp(-I0_star/L_I0_star_N3);
    
    if (fabs(map_I0_star) > 1.0) {
        cout << " ********************** map_I0_star out of bounds **********************" << endl;
        cout << "map_I0_star = " << map_I0_star << "   " << "I0_star = " << I0_star << "   " << "L_I0_star_N3 = " << L_I0_star_N3 << endl;
    }
    
    return map_I0_star;
}

//***************************************************************************************************
// This routine computes the exponential mapping, E_I0_eta, of the radiative energy density, used for
// the interpolative-based approximation of the third-order closing fluxes, N3_ijk
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Inverse_Exponential_Mapping(const long double &map_I0_star, const long double &L_I0_star_N3) {
    long double I0_star;
    
    I0_star = -L_I0_star_N3*log((1.0 - map_I0_star)/2.0);
    
    if (fabs(map_I0_star) > 1.0) {
        cout << " ********************** map_I0_star out of bounds **********************" << endl;
        cout << "map_I0_star = " << map_I0_star << "   " << "I0_star = " << I0_star << "   " << "L_I0_star_N3 = " << L_I0_star_N3 << endl;
    }
    
    return I0_star;
}

//*************************************************************************************
// This routine computes derivatives of the exponential mapping of the radiative energy 
// density, with respect to either the length scale or the radiative energy density.
// This will then be used to compute derivatives of our interpolative based approximation 
// of the Eddington factor for the non-gray M2 closure, with respect to either the length 
// scale or the radiative energy density
//*************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Diff_Exponential_Mapping(const long double &I0_star, const long double &L_I0_star_N3, const int &VAR_NUM) {
    long double dmap_I0_star_dvar;
    
    switch (VAR_NUM) {
        case INDEX_IO_STAR: // Then derivatives with respect to the energy density
            dmap_I0_star_dvar = (2.0/L_I0_star_N3)*exp(-I0_star/L_I0_star_N3);
            break;
        case INDEX_LENGTH_SCALE_IO_STAR: // Then derivatives with respect to the length scale
            dmap_I0_star_dvar = -2.0*(I0_star/(L_I0_star_N3*L_I0_star_N3))*exp(-I0_star/L_I0_star_N3);
            break;
        default:
            cout << "Specified value for VAR_NUM not correct !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    return dmap_I0_star_dvar;
}

// ******************************************************************************************
// This routine computes the interpolative-based approximation of the closing flux, N3_111
// ******************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_N3_111(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double N3_Fit;
    
    N3_Fit = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
    
    return N3_Fit;
}

long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_N3_112(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double N3_Fit;
    
    N3_Fit = Evaluate_N3_ijk(ratio_E, N1_2, -N1_1, N1_3, gam2, gam1, N3_122_ENTRY);
    
    return N3_Fit;
}

// ******************************************************************************************
// This routine computes the interpolative-based approximation of the closing flux, N3_122
// ******************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_N3_122(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double N3_Fit;
    long double norm_f_2;
    
    N3_Fit = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
    
    return N3_Fit;
}

long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_N3_123(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double gam3;
    long double N3_Fit;
    long double norm_f_2;
    
    N3_Fit = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_123_ENTRY);
    
    return N3_Fit;
}

long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_N3_222(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double N3_Fit;
    
    N3_Fit = Evaluate_N3_ijk(ratio_E, N1_2, -N1_1, N1_3, gam2, gam1, N3_111_ENTRY);
    
    return N3_Fit;
}

//*************************************************************************************************
// This routine computes the interpolative-based approximation of the third-order closing fluxes, 
// N3_ijk, in the case of the gray or the non-gray M2 closure
//*************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::Evaluate_N3_ijk(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {
    long double N3_ijk;
    long double f_interp; 
    long double norm_f_2;
                                                                                                    
    norm_f_2 = N1_1*N1_1 + N1_2*N1_2 + N1_3*N1_3;
    
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            f_interp = Evaluate_f_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
            N3_ijk = N1_1*(pow(N1_1,2) + f_interp*(1.0 - norm_f_2));
//             N3_ijk = pow(N1_1,3) + f_interp*N1_1/**(1.0 - norm_f_2)*/;
            break;
        case N3_122_ENTRY:
            f_interp = Evaluate_f_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
            N3_ijk = N1_1*(pow(N1_2, 2) + f_interp*(1.0 - norm_f_2));
//             N3_ijk = N1_1*pow(N1_2, 2) + f_interp*N1_1/**(1.0 - norm_f_2)*/;
            break;
        case N3_123_ENTRY:
            f_interp = Evaluate_f_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_123_ENTRY);
            N3_ijk = N1_1*N1_2*N1_3*f_interp;
            break;
        default:
            cout << "Incorrect value for CLOSING_FLUX_INDEX" << endl;
            exit(0);
            break;
    };
    
    // N3_ijk = f_interp;
    
    return N3_ijk;
}

//*************************************************************************************************
// This routine computes the interpolative-based approximation of the third-order closing fluxes, 
// N3_ijk, in the case of the gray or the non-gray M2 closure
//*************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::Evaluate_f_N3_ijk(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {
    long double f_N3_ijk, g_N3_ijk;
    double gam3;
    
    switch (closure_type) {
        case MOMENT_CLOSURE_PROJECTION:
            switch (CLOSING_FLUX_INDEX) {
                case N3_111_ENTRY:
                    f_N3_ijk = gam1;
                    break;
                case N3_122_ENTRY:
                    f_N3_ijk = gam2;
                    break;
                case N3_123_ENTRY:
                    cout << "Update this !!!!!!!!!!" << endl;
                    break;
                default:
                    cout << "Incorrect value for CLOSING_FLUX_INDEX" << endl;
                    exit(0);
                    break;
            };
            break;
        case MOMENT_CLOSURE_M2:
            switch (CLOSING_FLUX_INDEX) {
                case N3_111_ENTRY:
                    switch (flag_basis_type) {
                        case MONOMIAL_BASIS:
                            g_N3_ijk = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
                            break;
                        case CHEBYSHEV_BASIS:
                            g_N3_ijk = Evaluate_g_N3_ijk_Cheby_Basis(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
                            break;
                        default:
                            cout << "flag_basis_type not specified !!!!!!!!!!!!!!!!!!!!!" << endl;
                            exit(0);
                            break;
                    };
                    // cout << "g_N3_ijk = " << g_N3_ijk << endl;
//                     f_N3_ijk = gam1*(/*1.0 + (1.0 - gam1)**/g_N3_ijk);
                    f_N3_ijk = gam1*(1.0 + (1.0 - gam1)*g_N3_ijk);
                    // f_N3_ijk = g_N3_ijk;
                    break;
                case N3_122_ENTRY:
                    switch (flag_basis_type) {
                        case MONOMIAL_BASIS:
                            g_N3_ijk = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
                            break;
                        case CHEBYSHEV_BASIS:
                            g_N3_ijk = Evaluate_g_N3_ijk_Cheby_Basis(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
                            break;
                        default:
                            cout << "flag_basis_type not specified !!!!!!!!!!!!!!!!!!!!!" << endl;
                            exit(0);
                            break;
                    };
                    f_N3_ijk = gam2*(1.0 + (1.0 - gam2)*g_N3_ijk);
                    // f_N3_ijk = gam2*(1.0 + gam1*g_N3_ijk);
                    // f_N3_ijk = g_N3_ijk;
                    break;
                case N3_123_ENTRY:
                    switch (flag_basis_type) {
                        case MONOMIAL_BASIS:
                            g_N3_ijk = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_123_ENTRY);
                            break;
                        case CHEBYSHEV_BASIS:
                            g_N3_ijk = Evaluate_g_N3_ijk_Cheby_Basis(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_123_ENTRY);
                            break;
                        default:
                            cout << "flag_basis_type not specified !!!!!!!!!!!!!!!!!!!!!" << endl;
                            exit(0);
                            break;
                    };
                    gam3 = 1.0 - gam1 - gam2;
                    f_N3_ijk = 1.0 + gam1*gam2*gam3*g_N3_ijk;
                    break;
                default:
                    cout << "Incorrect value for CLOSING_FLUX_INDEX" << endl;
                    exit(0);
                    break;
            };
    
//     cout << "N1_1 = " << N1_1 << "   " << "N1_2 = " << N1_2 << "   " << "gam_1 = " << gam1 << "   " << "gam2 = " << gam2 << "   " << "f_N3_ijk = " << f_N3_ijk << endl;
    
//     f_N3_ijk = g_N3_ijk;
            break;
        default:
            cout << "Incorrect value for Closure Type" << endl;
            exit(0);
            break;
    };
    
    return f_N3_ijk;
}

//*************************************************************************************************
// This routine computes the interpolative-based approximation of the weighting function g_N3_ijk
// for our interpolative-based approximation of the third-order closing fluxes, N3_ijk in the case 
// of the gray M2 closure
//*************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::Evaluate_g_N3_ijk(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {    
    long double N3_Fit, N3_Fit_m, N3_Fit_f, N3_Fit_E, N3_Fit_gam1, N3_Fit_gam2, norm_f;
    long double norm_f_2;
    long double coeff_N3_ijk;
    int index_Coeffs;
    int N_Points_Gray_Cheby_without_Lebed = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
    long double x, y, z, x2, y2, z2;
    int Order_SH_temp = Order_SH/2;
    
    index_Coeffs = N_Points_Gray_Cheby_without_Lebed*N_Coeffs_SH - 1;
    
    norm_f_2 = pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2);
    norm_f = sqrt(norm_f_2);
    
    x = N1_1/norm_f;
    y = N1_2/norm_f;
    z = N1_3/norm_f;
    
    if (norm_f < 1.0e-8) {
        x = 1.0;
        y = 0.0;
        z = 0.0;
    }
    
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    // cout << "Order_SH = " << Order_SH << "  " << "Order_SH_temp = " << Order_SH_temp << "  " << "N_Points_E = " << N_Points_E << "  " << "N_Points_f = " << N_Points_f << "  " << "N_Points_f = " << N_Points_f << "  " << "N_Points_gam1 = " << N_Points_gam1 << endl;
    
    for (int l = Order_SH_temp; l >= 0; l--) {
        for (int m = Order_SH_temp - l; m >= 0; m--) {
            for (int i_fit_E = N_Points_E - 1; i_fit_E >= 0; i_fit_E--) {
                for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
                    for (int i_fit_gam1 = N_Points_gam1 - 1; i_fit_gam1 >= 0; i_fit_gam1--) {
                        for (int i_fit_gam2 = N_Points_gam1 - i_fit_gam1 - 1; i_fit_gam2 >= 0; i_fit_gam2--) {
                            
//                             cout << "l = " << l << "  " << "m = " << m << "  " << "i_fit_E = " << i_fit_E << "  " << "i_fit_f = " << i_fit_f << "  " << "i_fit_gam1 = " << i_fit_gam1 << "  " << "i_fit_gam2 = " << i_fit_gam2 << "  " << "Coefficient_Matrix_Fit_N3_111 = " << Coefficient_Matrix_Fit_N3_111[index_Coeffs] << endl;
                            
                            switch (CLOSING_FLUX_INDEX) {
                                    case N3_111_ENTRY:
                                        coeff_N3_ijk = Coefficient_Matrix_Fit_N3_111[index_Coeffs];
                                        // cout << "Coefficient_Matrix_Fit_N3_111 = " << coeff_N3_ijk << endl;
                                        break;
                                    case N3_122_ENTRY:
                                        coeff_N3_ijk = Coefficient_Matrix_Fit_N3_122[index_Coeffs];
                                        // cout << "Coefficient_Matrix_Fit_N3_122 = " << coeff_N3_ijk << endl;
                                        break;
                                    case N3_123_ENTRY:
                                        coeff_N3_ijk = Coefficient_Matrix_Fit_N3_123[index_Coeffs];
                                        break;
                                    default:
                                        cout << "Incorrect value for closing flux index" << endl;
                                        exit(0);
                                        break;
                            }
                            
                            if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
                                N3_Fit_gam2 = coeff_N3_ijk;
                            } else {
                                N3_Fit_gam2 = coeff_N3_ijk + N3_Fit_gam2*gam2;
                            }
                            
                            //                             cout << "Coefficient_Matrix_Fit_N3_111 = " << Coefficient_Matrix_Fit_N3_111[index_Coeffs] << endl;
                            index_Coeffs--;
                        }
                        if (i_fit_gam1 == N_Points_gam1 - 1) {
                            N3_Fit_gam1 = N3_Fit_gam2;
                        } else {
                            N3_Fit_gam1 = N3_Fit_gam2 + N3_Fit_gam1*gam1;
                        }
                    }
                    if (i_fit_f == N_Points_f - 1) {
                        N3_Fit_f = N3_Fit_gam1;
                    } else {
                        N3_Fit_f = N3_Fit_gam1 + N3_Fit_f*norm_f_2;
                    } 
                }
                if (i_fit_E == N_Points_E - 1) {
                    N3_Fit_E = N3_Fit_f;
                } else {
                    N3_Fit_E = N3_Fit_f + N3_Fit_E*ratio_E;
                } 
            }
            if (m == Order_SH_temp - l) {
                N3_Fit_m = N3_Fit_E;
            } else {
                N3_Fit_m = N3_Fit_E + N3_Fit_m*y2;
                // N3_Fit_m = N3_Fit_E + N3_Fit_m*z2;
                // N3_Fit_m = N3_Fit_E + N3_Fit_m*x2;
            }
        }
        if (l == Order_SH_temp) {
            N3_Fit = N3_Fit_m;
        } else {
            N3_Fit = N3_Fit_m + N3_Fit*x2;
            // N3_Fit = N3_Fit_m + N3_Fit*z2;
        }
    }
    
    
    
    if (index_Coeffs != -1) {
        cout << "index_Coeffs = " << index_Coeffs << endl;
        exit(0);
    }
    
    // cout << "index_Coeffs = " << index_Coeffs << endl;
    
    // cout << "N3_Fit_m = " << N3_Fit_m << "  " << "N3_Fit_E = " << N3_Fit_E << "  " << "N3_Fit_f = " << N3_Fit_f << "  " << "N3_Fit_gam1 = " << N3_Fit_gam1 << "  " << "N3_Fit_gam2 = " << N3_Fit_gam2 << endl;
    
    return N3_Fit;
}

long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_g_N3_ijk_Cheby_Basis(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &VAR_NUM) {
    long double N3_Fit, N3_Fit_m, N3_Fit_f, N3_Fit_E, N3_Fit_gam1, N3_Fit_gam2, norm_f;
    long double norm_f_2;
    int index;
    long double x, y, z;
    long double poly_map_I0_star, poly_norm_f, poly_gam1_gam2, poly_SH, Coeff_val;
    int index_triangle_iter;
    int N_Points_Gray_Cheby_without_Lebed = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
    
    norm_f_2 = pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2);
    norm_f = sqrt(norm_f_2);
    
    x = N1_1/norm_f;
    y = N1_2/norm_f;
    z = N1_3/norm_f;
    
    if (norm_f < 1.0e-8) {
        x = 1.0;
        y = 0.0;
        z = 0.0;
    }
    
    N3_Fit = 0.0;
    for (int index_SH = 0; index_SH < N_Coeffs_SH; index_SH++) {
        for (int i_fit_E = 0; i_fit_E < N_Points_E; i_fit_E++) {
            for (int i_fit_f = 0; i_fit_f < N_Points_f; i_fit_f++) {
                index_triangle_iter = 0;
                for (int i_fit_gam1 = 0; i_fit_gam1 < N_Points_gam1; i_fit_gam1++) {
                    for (int i_fit_gam2 = 0; i_fit_gam2 < N_Points_gam1 - i_fit_gam1; i_fit_gam2++) {
                        index = (index_SH*N_Points_E + i_fit_E)*N_Points_f + i_fit_f;
                        index = index*N_Points_Triangle_gam1_gam2 + index_triangle_iter;
                        switch (VAR_NUM) {
                            case N3_111_ENTRY:
                                Coeff_val = Coefficient_Matrix_Fit_N3_111_Cheby_Basis[index];
                                break;
                            case N3_122_ENTRY:
                                Coeff_val = Coefficient_Matrix_Fit_N3_122_Cheby_Basis[index];
                                break;
                            case N3_123_ENTRY:
                                Coeff_val = Coefficient_Matrix_Fit_N3_123_Cheby_Basis[index];
                                break;
                            default:
                                cout << "Values for VAR_NUM not specified !!!!!!!!!!!!!!!!!!!!!" << endl;
                                exit(0);
                                break;
                        }
                        
                        poly_SH = Cartesian_Spherical_harmonics(x, y, z, Array_l_SH[index_SH], Array_m_SH[index_SH]);
                        poly_map_I0_star = Chebyshev_Polynomial_Basis(ratio_E, i_fit_E);
                        poly_norm_f = Chebyshev_Polynomial_Basis(norm_f, 2*i_fit_f);
                        poly_gam1_gam2 = Proriol_Polynomials(gam1, gam2, i_fit_gam1, i_fit_gam2);
                            
                        N3_Fit += Coeff_val*poly_map_I0_star*poly_norm_f*poly_gam1_gam2*poly_SH;
                        index_triangle_iter++;
                    }
                }
            }
        }
    }
    
    if (index != N_Points_Gray_Cheby_without_Lebed*N_Coeffs_SH - 1) {
        cout << "index = " << index << "  " << "N_Points_Gray_Cheby_without_Lebed*N_Coeffs_SH - 1 = " << N_Points_Gray_Cheby_without_Lebed*N_Coeffs_SH - 1 << endl;
        exit(0);
    }
    
    return N3_Fit;
}

//*************************************************************************************************
// Derivaties of thirs-order closing fluxes
//*************************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby::dN3_111_2D_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_111_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double dgam1_dIn, dgam2_dIn, dB_dIn;
    long double norm_f_2;
                                                                                                    
    norm_f_2 = N1_1*N1_1 + N1_2*N1_2 + N1_3*N1_3;
    
     if (fabs(1.0 - norm_f_2) < TOLER) {
        // Then N3_111 = N1_1 N2_11
        d_N3_RT_dU[1] = N1_1*N1_1;
        d_N3_RT_dU[2] = 0.0;
        d_N3_RT_dU[3] = N1_1;
        d_N3_RT_dU[4] = 0.0;
        d_N3_RT_dU[5] = 0.0;
    } else {
        dN3_ijk_2D_dU(dN3_111_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
        // dN3_ijk_2D_dU_Finite_Difference(dN3_111_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
        
        // Now compute full derivatives using product rule
        // Compute {d N3}{d N1_1}
        dgam1_dIn = dgam1_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        dgam2_dIn = dgam2_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        dB_dIn = dB_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[1] = dN3_111_RT.dN1_1 + dN3_111_RT.dgam1*dgam1_dIn + dN3_111_RT.dgam2*dgam2_dIn + dN3_111_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N1_2}
        dgam1_dIn = dgam1_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        dgam2_dIn = dgam2_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        dB_dIn = dB_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[2] = dN3_111_RT.dN1_2 + dN3_111_RT.dgam1*dgam1_dIn + dN3_111_RT.dgam2*dgam2_dIn + dN3_111_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N2_11}
        dgam1_dIn = dgam1_dU(INDEX_N2_11, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[3] = dN3_111_RT.dgam1*dgam1_dIn;
        
        // Compute {d N3}{d N2_12}
        dB_dIn = dB_dU(INDEX_N2_12, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[4] = dN3_111_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N2_22}
        dgam2_dIn = dgam2_dU(INDEX_N2_22, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[5] = dN3_111_RT.dgam2*dgam2_dIn;
    }
}

//***************************************************************************************************
// This routine computes the derivatives of the interpolative-based approximation of the third-order 
// closing fluxes, N3_112, with respect to the conserved solutions in the case of the gray M2 closure
//***************************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby::dN3_112_2D_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_112_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double dgam1_dIn, dgam2_dIn, dB_dIn;
    record_d_N3_ijk_RT dN3_122_RT;
    long double norm_f_2;
                                                                                                    
    norm_f_2 = N1_1*N1_1 + N1_2*N1_2 + N1_3*N1_3;
    
     if (fabs(1.0 - norm_f_2) < TOLER) {
        // Then N3_112 = N1_2 N2_11
        d_N3_RT_dU[1] = 0.0;
        d_N3_RT_dU[2] = N1_1*N1_1;
        d_N3_RT_dU[3] = N1_2;
        d_N3_RT_dU[4] = 0.0;
        d_N3_RT_dU[5] = 0.0;
    } else {
        dN3_ijk_2D_dU(dN3_122_RT, d_r_I0_star_N3, ratio_E, N1_2, -N1_1, N1_3, gam2, gam1, N3_122_ENTRY);
        // dN3_ijk_2D_dU_Finite_Difference(dN3_122_RT, d_r_I0_star_N3, ratio_E, N1_2, -N1_1, N1_3, gam2, gam1, N3_122_ENTRY);
        
        dN3_112_RT.dI0 = dN3_122_RT.dI0;
        dN3_112_RT.dN1_1 = -dN3_122_RT.dN1_2;
        dN3_112_RT.dN1_2 = dN3_122_RT.dN1_1;
        dN3_112_RT.dgam1 = dN3_122_RT.dgam2;
        dN3_112_RT.dgam2 = dN3_122_RT.dgam1;
        
        // Now compute full derivatives using product rule
        // Compute {d N3}{d N1_1}
        dgam1_dIn = dgam1_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        dgam2_dIn = dgam2_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        dB_dIn = dB_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[1] = dN3_112_RT.dN1_1 + dN3_112_RT.dgam1*dgam1_dIn + dN3_112_RT.dgam2*dgam2_dIn + dN3_112_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N1_2}
        dgam1_dIn = dgam1_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        dgam2_dIn = dgam2_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        dB_dIn = dB_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[2] = dN3_112_RT.dN1_2 + dN3_112_RT.dgam1*dgam1_dIn + dN3_112_RT.dgam2*dgam2_dIn + dN3_112_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N2_11}
        dgam1_dIn = dgam1_dU(INDEX_N2_11, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[3] = dN3_112_RT.dgam1*dgam1_dIn;
        
        // Compute {d N3}{d N2_12}
        dB_dIn = dB_dU(INDEX_N2_12, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[4] =  dN3_112_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N2_22}
        dgam2_dIn = dgam2_dU(INDEX_N2_22, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[5] = dN3_112_RT.dgam2*dgam2_dIn;
    }
}

//***************************************************************************************************
// This routine computes the derivatives of the interpolative-based approximation of the third-order 
// closing fluxes, N3_122, with respect to the conserved solutions in the case of the gray M2 closure
//***************************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby::dN3_122_2D_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_122_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double dgam1_dIn, dgam2_dIn, dB_dIn;
    long double norm_f_2;
                                                                                                    
    norm_f_2 = N1_1*N1_1 + N1_2*N1_2 + N1_3*N1_3;
    
    if (fabs(1.0 - norm_f_2) < TOLER) {
        // Then N3_122 = N1_1 N2_22
        d_N3_RT_dU[1] = N1_2*N1_2;
        d_N3_RT_dU[2] = 0.0;
        d_N3_RT_dU[3] = 0.0;
        d_N3_RT_dU[4] = 0.0;
        d_N3_RT_dU[5] = N1_1;
    } else {
        dN3_ijk_2D_dU(dN3_122_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
        // dN3_ijk_2D_dU_Finite_Difference(dN3_122_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
        
        // Now compute full derivatives using product rule
        // Compute {d N3}{d N1_1}
        dgam1_dIn = dgam1_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        dgam2_dIn = dgam2_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        dB_dIn = dB_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[1] = dN3_122_RT.dN1_1 + dN3_122_RT.dgam1*dgam1_dIn + dN3_122_RT.dgam2*dgam2_dIn + dN3_122_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N1_2}
        dgam1_dIn = dgam1_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        dgam2_dIn = dgam2_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        dB_dIn = dB_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[2] = dN3_122_RT.dN1_2 + dN3_122_RT.dgam1*dgam1_dIn + dN3_122_RT.dgam2*dgam2_dIn + dN3_122_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N2_11}
        dgam1_dIn = dgam1_dU(INDEX_N2_11, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[3] = dN3_122_RT.dgam1*dgam1_dIn;
        
        // Compute {d N3}{d N2_12}
        dB_dIn = dB_dU(INDEX_N2_12, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[4] = dN3_122_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N2_22}
        dgam2_dIn = dgam2_dU(INDEX_N2_22, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[5] = dN3_122_RT.dgam2*dgam2_dIn;
    }
}

//***************************************************************************************************
// This routine computes the derivatives of the interpolative-based approximation of the third-order 
// closing fluxes, N3_222, with respect to the conserved solutions in the case of the gray M2 closure
//***************************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby::dN3_222_2D_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_222_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double dgam1_dIn, dgam2_dIn, dB_dIn;
    record_d_N3_ijk_RT dN3_111_RT;
    long double norm_f_2;
                                                                                                    
    norm_f_2 = N1_1*N1_1 + N1_2*N1_2 + N1_3*N1_3;
    
    if (fabs(1.0 - norm_f_2) < TOLER) {
        // Then N3_222 = N1_2 N2_22
        d_N3_RT_dU[1] = 0.0;
        d_N3_RT_dU[2] = N1_2*N1_2;
        d_N3_RT_dU[3] = 0.0;
        d_N3_RT_dU[4] = 0.0;
        d_N3_RT_dU[5] = N1_2;
    } else {
        dN3_ijk_2D_dU(dN3_111_RT, d_r_I0_star_N3, ratio_E, N1_2, -N1_1, N1_3, gam2, gam1, N3_111_ENTRY);
        // dN3_ijk_2D_dU_Finite_Difference(dN3_111_RT, d_r_I0_star_N3, ratio_E, N1_2, -N1_1, N1_3, gam2, gam1, N3_111_ENTRY);
        
        dN3_222_RT.dI0 = dN3_111_RT.dI0;
        dN3_222_RT.dN1_1 = -dN3_111_RT.dN1_2;
        dN3_222_RT.dN1_2 = dN3_111_RT.dN1_1;
        dN3_222_RT.dgam1 = dN3_111_RT.dgam2;
        dN3_222_RT.dgam2 = dN3_111_RT.dgam1;
        
        // Now compute full derivatives using product rule
        // Compute {d N3}{d N1_1}
        dgam1_dIn = dgam1_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        dgam2_dIn = dgam2_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        dB_dIn = dB_dU(INDEX_N1_1, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[1] = dN3_222_RT.dN1_1 + dN3_222_RT.dgam1*dgam1_dIn + dN3_222_RT.dgam2*dgam2_dIn + dN3_222_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N1_2}
        dgam1_dIn = dgam1_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        dgam2_dIn = dgam2_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        dB_dIn = dB_dU(INDEX_N1_2, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[2] = dN3_222_RT.dN1_2 + dN3_222_RT.dgam1*dgam1_dIn + dN3_222_RT.dgam2*dgam2_dIn + dN3_222_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N2_11}
        dgam1_dIn = dgam1_dU(INDEX_N2_11, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[3] = dN3_222_RT.dgam1*dgam1_dIn;
        
        // Compute {d N3}{d N2_12}
        dB_dIn = dB_dU(INDEX_N2_12, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[4] = dN3_222_RT.dB*dB_dIn;
        
        // Compute {d N3}{d N2_22}
        dgam2_dIn = dgam2_dU(INDEX_N2_22, N1_1, N1_2, N1_3, gam1, gam2);
        d_N3_RT_dU[5] = dN3_222_RT.dgam2*dgam2_dIn;
    }
}

//***************************************************************************************************
// This routine computes the derivatives of the interpolative-based approximation of the third-order 
// closing fluxes, N3_ijk, with respect to the conserved solutions in the case of the gray M2 closure
//***************************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby::dN3_ijk_2D_dU(record_d_N3_ijk_RT &dN3_ijk_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {
    long double f_interp, d_f_interp, d_gam1_dU;
    long double norm_f_2;
    long double dgam1_dIn, dgam2_dIn;
    long double df_N3_ijk_dratioI0, df_N3_ijk_dnorm_f_2, df_N3_ijk_dx2, df_N3_ijk_dy2, df_N3_ijk_dgam1, df_N3_ijk_dgam2;
    long double dnorm_f2_dN1_1, dnorm_f2_dN1_2, dx2_dN1_1, dx2_dN1_2, dy2_dN1_1, dy2_dN1_2;
    
    norm_f_2 = pow(N1_1,2) + pow(N1_2,2) + pow(N1_3,2);
    
    switch(Problem_Type) {
        case GRAY:
            df_N3_ijk_dratioI0 = 0.0;
            break;
        case NON_GRAY:
            df_N3_ijk_dratioI0 = df_N3_ijk_interp_dU(INDEX_RATIO_I0, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
            break;
    }
    
    df_N3_ijk_dnorm_f_2 = df_N3_ijk_interp_dU(INDEX_NORM_f_2, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    df_N3_ijk_dx2 = df_N3_ijk_interp_dU(INDEX_x_SH_2, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    df_N3_ijk_dy2 = df_N3_ijk_interp_dU(INDEX_y_SH_2, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    df_N3_ijk_dgam1 = df_N3_ijk_interp_dU(INDEX_GAM1, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    df_N3_ijk_dgam2 = df_N3_ijk_interp_dU(INDEX_GAM2, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    
    dnorm_f2_dN1_1 = dnorm_f2_dU(INDEX_N1_1, N1_1, N1_2, N1_3);
    dnorm_f2_dN1_2 = dnorm_f2_dU(INDEX_N1_2, N1_1, N1_2, N1_3);
    
    dx2_dN1_1 = dx2_dU(INDEX_N1_1, N1_1, N1_2, N1_3);
    dx2_dN1_2 = dx2_dU(INDEX_N1_2, N1_1, N1_2, N1_3);
    
    dy2_dN1_1 = dy2_dU(INDEX_N1_1, N1_1, N1_2, N1_3);
    dy2_dN1_2 = dy2_dU(INDEX_N1_2, N1_1, N1_2, N1_3);
    
    f_interp = Evaluate_f_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    
    // cout << endl;
    // cout << "dnorm_f2_dN1_1 = " << dnorm_f2_dN1_1 << "  " << "dx2_dN1_1 = " << dx2_dN1_1 << "  " << "dy2_dN1_1 = " << dy2_dN1_1 << endl;
    // cout << endl;
    
    // cout << "df_N3_ijk_dx2 = " << df_N3_ijk_dx2 << "  " << "df_N3_ijk_dy2 = " << df_N3_ijk_dy2 << endl;
    
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            // Compute {\partial N3}{\partial I0}
            dN3_ijk_RT.dI0 = N1_1*(1.0 - norm_f_2) * df_N3_ijk_dratioI0;
            dN3_ijk_RT.dI0 *= d_r_I0_star_N3.dI0_star;
            
            // Compute {\partial N3}{\partial N1_1}
            d_f_interp = df_N3_ijk_dnorm_f_2*dnorm_f2_dN1_1 + df_N3_ijk_dx2*dx2_dN1_1 + df_N3_ijk_dy2*dy2_dN1_1;
            dN3_ijk_RT.dN1_1 = 3.0*pow(N1_1,2) + f_interp*(1.0 - norm_f_2 - 2.0*pow(N1_1,2)) + N1_1*(1.0 - norm_f_2)*d_f_interp;
            
            // Compute {\partial N3}{\partial N1_2}
            d_f_interp = df_N3_ijk_dnorm_f_2*dnorm_f2_dN1_2 + df_N3_ijk_dx2*dx2_dN1_2 + df_N3_ijk_dy2*dy2_dN1_2;
            dN3_ijk_RT.dN1_2 = N1_1*(-2.0*N1_2*f_interp + (1.0 - norm_f_2)*d_f_interp);
            
            // Compute {\partial N3}{\partial gam1}
            d_f_interp = df_N3_ijk_dgam1;
            dN3_ijk_RT.dgam1 = N1_1*(1.0 - norm_f_2)*d_f_interp;
            
            // Compute {\partial N3}{\partial gam2}
            d_f_interp = df_N3_ijk_dgam2;
            dN3_ijk_RT.dgam2 = N1_1*(1.0 - norm_f_2)*d_f_interp;
            break;
        case N3_122_ENTRY:
            // Compute {\partial N3}{\partial I0}
            dN3_ijk_RT.dI0 = N1_1*(1.0 - norm_f_2) * df_N3_ijk_dratioI0;
            dN3_ijk_RT.dI0 *= d_r_I0_star_N3.dI0_star;
            
            // Compute {\partial N3}{\partial N1_1}
            d_f_interp = df_N3_ijk_dnorm_f_2*dnorm_f2_dN1_1 + df_N3_ijk_dx2*dx2_dN1_1 + df_N3_ijk_dy2*dy2_dN1_1;
            dN3_ijk_RT.dN1_1 = pow(N1_2,2) + f_interp*(1.0 - norm_f_2 - 2.0*pow(N1_1,2)) + N1_1*(1.0 - norm_f_2)*d_f_interp;
            
            // Compute {\partial N3}{\partial N1_2}
            d_f_interp = df_N3_ijk_dnorm_f_2*dnorm_f2_dN1_2 + df_N3_ijk_dx2*dx2_dN1_2 + df_N3_ijk_dy2*dy2_dN1_2;
            dN3_ijk_RT.dN1_2 = N1_1*(2.0*N1_2*(1.0 - f_interp) + (1.0 - norm_f_2)*d_f_interp);
            
            // Compute {\partial N3}{\partial gam1}
            d_f_interp = df_N3_ijk_dgam1;
            dN3_ijk_RT.dgam1 = N1_1*(1.0 - norm_f_2)*d_f_interp;
            
            // Compute {\partial N3}{\partial gam2}
            d_f_interp = df_N3_ijk_dgam2;
            dN3_ijk_RT.dgam2 = N1_1*(1.0 - norm_f_2)*d_f_interp;
            break;
        case N3_123_ENTRY:
            cout << "Derivatives for N3_123 not implemented yet !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
        default:
            cout << "Invalid value for CLOSING_FLUX_INDEX in dN3_ijk_2D_dU !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
}

void N3_Non_Gray_M2_3D_RT_Cheby::dN3_ijk_2D_dU_Finite_Difference(record_d_N3_ijk_RT &dN3_ijk_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {
    long double epsilon = 1.0e-6;
    long double N3_ijk_L, N3_ijk_R;
    
    N3_ijk_L = Evaluate_N3_ijk(ratio_E - epsilon, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    N3_ijk_R = Evaluate_N3_ijk(ratio_E + epsilon, N1_1, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    dN3_ijk_RT.dI0 = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    dN3_ijk_RT.dI0 *= d_r_I0_star_N3.dI0_star;
    
    N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1 - epsilon, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1 + epsilon, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    dN3_ijk_RT.dN1_1 = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    
    N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1, N1_2 - epsilon, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1, N1_2 + epsilon, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    dN3_ijk_RT.dN1_2 = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    
    N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3 - epsilon, gam1, gam2, CLOSING_FLUX_INDEX);
    N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3 + epsilon, gam1, gam2, CLOSING_FLUX_INDEX);
    dN3_ijk_RT.dN1_3 = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    
    N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1 - epsilon, gam2, CLOSING_FLUX_INDEX);
    N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1 + epsilon, gam2, CLOSING_FLUX_INDEX);
    dN3_ijk_RT.dgam1 = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    
    N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 - epsilon, CLOSING_FLUX_INDEX);
    N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 + epsilon, CLOSING_FLUX_INDEX);
    dN3_ijk_RT.dgam2 = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    
    // cout << "dq_dN1_1 = " << dN3_ijk_RT.dN1_1 << "  " << "dq_dN1_2 = " << dN3_ijk_RT.dN1_2 << "  " << "dq_dN1_3 = " << dN3_ijk_RT.dN1_3 << endl;
}

//*************************************************************************************************
// This routine computes the derivatives of the interpolative-based approximation of the weighting 
// function f_N3_ijk for our interpolative-based approximation of the third-order closing fluxes, 
// N3_ijk, in the case of the gray M2 closure
//*************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::df_N3_ijk_interp_dU(const int &index_var_fit, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {
    long double g_N3_ijk, df_N3_ijk, dg_N3_ijk;
    
    switch (closure_type) {
        case MOMENT_CLOSURE_PROJECTION:
            switch (CLOSING_FLUX_INDEX) {
                case N3_111_ENTRY:
                    if (index_var_fit == INDEX_GAM1) {
                        df_N3_ijk = 1.0;
                    } else {
                        df_N3_ijk = 0.0;
                    }
                    break;
                case N3_122_ENTRY:
                    if (index_var_fit == INDEX_GAM2) {
                        df_N3_ijk = 1.0;
                    } else {
                        df_N3_ijk = 0.0;
                    }
                    break;
                case N3_123_ENTRY:
                    
                    break;
                default:
                    cout << "Incorrect value for CLOSING_FLUX_INDEX" << endl;
                    exit(0);
                    break;
            };
            break;
        case MOMENT_CLOSURE_M2:
            switch (CLOSING_FLUX_INDEX) {
                case N3_111_ENTRY:
                    g_N3_ijk = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
                    dg_N3_ijk = dg_N3_ijk_interp_dU(index_var_fit, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
                    if (index_var_fit == INDEX_GAM1) {
                        df_N3_ijk = 1.0 + (1.0 - 2.0*gam1)*g_N3_ijk + gam1*(1.0 - gam1)*dg_N3_ijk;
                    } else {
                        df_N3_ijk = gam1*(1.0 - gam1)*dg_N3_ijk;
                    }
                    break;
                case N3_122_ENTRY:
                    g_N3_ijk = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
                    dg_N3_ijk = dg_N3_ijk_interp_dU(index_var_fit, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
                    if (index_var_fit == INDEX_GAM1) {
                        df_N3_ijk = gam2*(g_N3_ijk + gam1*dg_N3_ijk);
                    } else if (index_var_fit == INDEX_GAM2) {
                        df_N3_ijk = 1.0 + gam1*g_N3_ijk + gam1*gam2*dg_N3_ijk;
                    } else {
                        df_N3_ijk = gam1*gam2*dg_N3_ijk;
                    }
                    break;
                case N3_123_ENTRY:
                    
                    break;
                default:
                    cout << "Incorrect value for CLOSING_FLUX_INDEX" << endl;
                    exit(0);
                    break;
            };
            break;
        default:
            cout << "Incorrect value for Closure Type" << endl;
            exit(0);
            break;
    };
    
    return df_N3_ijk;
}

//*************************************************************************************************
// This routine computes the derivatives of the interpolative-based approximation of the weighting 
// function g_N3_ijk for our interpolative-based approximation of the third-order closing fluxes, 
// N3_ijk, in the case of the gray M2 closure
//*************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::dg_N3_ijk_interp_dU(const int &index_var_fit, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {
    long double f_interp, f_interp_m, f_interp_E, f_interp_f, f_interp_gam1, f_interp_gam2;
    long double norm_f_2, norm_f;
    long double coeff_N3_ijk;
    int index_Coeffs;
    int N_Points_Gray_Cheby_without_Lebed = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
    long double x, y, z, x2, y2, z2;
    int Order_SH_temp = Order_SH/2;
    
    index_Coeffs = N_Points_Gray_Cheby_without_Lebed*N_Coeffs_SH - 1;
    
    norm_f_2 = pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2);
    norm_f = sqrt(norm_f_2);
    
    x = N1_1/norm_f;
    y = N1_2/norm_f;
    z = N1_3/norm_f;
    
    if (norm_f < 1.0e-8) {
        x = 1.0;
        y = 0.0;
        z = 0.0;
    }
    
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    for (int l = Order_SH_temp; l >= 0; l--) {
        for (int m = Order_SH_temp - l; m >= 0; m--) {
            for (int i_fit_E = N_Points_E - 1; i_fit_E >= 0; i_fit_E--) {
                for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
                    for (int i_fit_gam1 = N_Points_gam1 - 1; i_fit_gam1 >= 0; i_fit_gam1--) {
                        for (int i_fit_gam2 = N_Points_gam1 - i_fit_gam1 - 1; i_fit_gam2 >= 0; i_fit_gam2--) {
                            // Choose appropriate coefficients depending of what third-order closing flux is of interest
                            switch (CLOSING_FLUX_INDEX) {
                                    case N3_111_ENTRY:
                                        coeff_N3_ijk = Coefficient_Matrix_Fit_N3_111[index_Coeffs];
                                        break;
                                    case N3_122_ENTRY:
                                        coeff_N3_ijk = Coefficient_Matrix_Fit_N3_122[index_Coeffs];
                                        break;
                                    case N3_123_ENTRY:
                                        coeff_N3_ijk = Coefficient_Matrix_Fit_N3_123[index_Coeffs];
                                        break;
                                    default:
                                        cout << "Incorrect value for closing flux index" << endl;
                                        exit(0);
                                        break;
                            }
                            
                            if (index_var_fit == INDEX_GAM2) {
                                if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
                                    f_interp_gam2 = i_fit_gam2*coeff_N3_ijk;
                                } else if (i_fit_gam2 >= 1) {
                                    f_interp_gam2 = i_fit_gam2*coeff_N3_ijk + f_interp_gam2*gam2;
                                }
                            } else {
                                if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
                                    f_interp_gam2 = coeff_N3_ijk;
                                } else {
                                    f_interp_gam2 = coeff_N3_ijk + f_interp_gam2*gam2;
                                }
                            }
                            index_Coeffs--;
                        }
                        
                        if (index_var_fit == INDEX_GAM1) {
                            if (i_fit_gam1 == N_Points_gam1 - 1) {
                                f_interp_gam1 = i_fit_gam1*f_interp_gam2;
                            } else if (i_fit_gam1 >= 1) {
                                f_interp_gam1 = i_fit_gam1*f_interp_gam2 + f_interp_gam1*gam1;
                            }
                        } else {
                            if (i_fit_gam1 == N_Points_gam1 - 1) {
                                f_interp_gam1 = f_interp_gam2;
                            } else {
                                f_interp_gam1 = f_interp_gam2 + f_interp_gam1*gam1;
                            }
                        }
                    }
                    
                    if (index_var_fit == INDEX_NORM_f_2) {
                        if (i_fit_f == N_Points_f - 1) {
                            f_interp_f = i_fit_f*f_interp_gam1;
                        } else if (i_fit_f >= 1) {
                            f_interp_f = i_fit_f*f_interp_gam1 + f_interp_f*norm_f_2;
                        }
                    } else {
                        if (i_fit_f == N_Points_f - 1) {
                            f_interp_f = f_interp_gam1;
                        } else {
                            f_interp_f = f_interp_gam1 + f_interp_f*norm_f_2;
                        }
                    } 
                }
                if (index_var_fit == INDEX_RATIO_I0) {
                    if (i_fit_E == N_Points_E - 1) {
                        f_interp_E = i_fit_E*f_interp_f;
                    } else if (i_fit_E >= 1) {
                        f_interp_E = i_fit_E*f_interp_f + f_interp_E*ratio_E;
                    }
                } else {
                    if (i_fit_E == N_Points_E - 1) {
                        f_interp_E = f_interp_f;
                    } else {
                        f_interp_E = f_interp_f + f_interp_E*ratio_E;
                    }
                }
            }
            
            if (index_var_fit == INDEX_y_SH_2) {
                if (m == Order_SH_temp - l) {
                    f_interp_m = m*f_interp_E;
                } else if (m >= 1) {
                    // f_interp_m = m*f_interp_E + f_interp_m*x2;
                    f_interp_m = m*f_interp_E + f_interp_m*y2;
                }
            } else {
                if (m == Order_SH_temp - l) {
                    f_interp_m = f_interp_E;
                } else {
                    // f_interp_m = f_interp_E + f_interp_m*x2;
                    f_interp_m = f_interp_E + f_interp_m*y2;
                }
            }
        }
        
        if (index_var_fit == INDEX_x_SH_2) {   
            if (l == Order_SH_temp) {
                f_interp = l*f_interp_m;
            } else if (l >= 1) {
                // f_interp = l*f_interp_m + f_interp*z2;
                f_interp = l*f_interp_m + f_interp*x2;
            }
        } else {   
            if (l == Order_SH_temp) {
                f_interp = f_interp_m;
            } else {
                 // f_interp = f_interp_m + f_interp*z2;
                 f_interp = f_interp_m + f_interp*x2;
            }
        }
    }
    
    if (index_Coeffs != -1) {
        cout << "index_Coeffs = " << index_Coeffs << endl;
        exit(0);
    }
    
    return f_interp;
}

//***************************************************************************************************
// This routine computes the derivatives of the norm of the first-order angular moments with respect 
// to components of the latter
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::dnorm_f2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3) {
    long double d_norm_f_dU = 0.0;
    switch (index_U) {
        case INDEX_N1_1: // Then N1_1
            d_norm_f_dU = 2.0*N1_1;
            break;
        case INDEX_N1_2: // Then N1_2
            d_norm_f_dU = 2.0*N1_2;
            break;
        default:
            cout << "Ivalid value for index_U in dnorm_f2_dU !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    return d_norm_f_dU;
}

//***************************************************************************************************
// This routine computes the derivatives of the direction cosines of the first-order angular moments
// with respect to components of the latter
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::dx2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3) {
    long double d_x2_dU;
    long double norm_f_2 = N1_1*N1_1 + N1_2*N1_2 + N1_3*N1_3;
    switch (index_U) {
        case INDEX_N1_1: // Then N1_1
            if (norm_f_2 < 1.0e-7) {
                d_x2_dU = 0.0;
            } else {
                d_x2_dU = 2.0*N1_1*(N1_2*N1_2 + N1_3*N1_3)/pow(norm_f_2, 2);
            }
            break;
        case INDEX_N1_2: // Then N1_2
            if (norm_f_2 < 1.0e-7) {
                d_x2_dU = 0.0;
            } else {
                d_x2_dU = -2.0*N1_2*(N1_1*N1_1)/pow(norm_f_2, 2);
            }
            break;
        default:
            cout << "Ivalid value for index_U in dx2_dU !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    return d_x2_dU;
}

//***************************************************************************************************
// This routine computes the derivatives of the direction cosines of the first-order angular moments
// with respect to components of the latter
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::dy2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3) {
    long double d_y2_dU;
    long double norm_f_2 = N1_1*N1_1 + N1_2*N1_2 + N1_3*N1_3;
    switch (index_U) {
        case INDEX_N1_1: // Then N1_1
            if (norm_f_2 < 1.0e-7) {
                d_y2_dU = 0.0;
            } else {
                d_y2_dU = -2.0*N1_1*(N1_2*N1_2)/pow(norm_f_2, 2);
            }
            break;
        case INDEX_N1_2: // Then N1_2
            if (norm_f_2 < 1.0e-7) {
                d_y2_dU = 0.0;
            } else {
                d_y2_dU = 2.0*N1_2*(N1_1*N1_1 + N1_3*N1_3)/pow(norm_f_2, 2);
            }
            break;
        default:
            cout << "Ivalid value for index_U in dy2_dU !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    return d_y2_dU;
}

//***************************************************************************************************
// This routine computes the derivatives of the direction cosines of the first-order angular moments
// with respect to components of the latter
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::dz2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3) {
    long double d_z2_dU;
    long double norm_f_2 = N1_1*N1_1 + N1_2*N1_2 + N1_3*N1_3;
    switch (index_U) {
        case INDEX_N1_1: // Then N1_1
            if (norm_f_2 < 1.0e-7) {
                d_z2_dU = 0.0;
            } else {
                d_z2_dU = -2.0*N1_1*(N1_3*N1_3)/pow(norm_f_2, 2);
            }
            break;
        case INDEX_N1_2: // Then N1_2
            if (norm_f_2 < 1.0e-7) {
                d_z2_dU = 0.0;
            } else {
                d_z2_dU = -2.0*N1_2*(N1_3*N1_3)/pow(norm_f_2, 2);
            }
            break;
        default:
            cout << "Ivalid value for index_U in dz2_dU !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    return d_z2_dU;
}

//***************************************************************************************************
// This routine computes the derivatives of the eigenvalues of the covariance matrix with respect to 
// the conserved solution variables
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::dgam1_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double d_gam1_dU;
    long double norm_f_2;
    norm_f_2 = pow(N1_1,2) + pow(N1_2,2) + pow(N1_3,2);
    
    switch (index_U) {
        case INDEX_I0:
            d_gam1_dU = (pow(N1_1, 2) - gam1*(1.0 + norm_f_2))/(1.0 - norm_f_2);
            cout << "Fix this !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
        case INDEX_N1_1: // Then N1_1
            d_gam1_dU = -2.0*N1_1*(1.0 - gam1)/(1.0 - norm_f_2);
            break;
        case INDEX_N1_2: // Then N1_2
            d_gam1_dU = 2.0*N1_2*gam1/(1.0 - norm_f_2);
            break;
        case INDEX_N2_11: // Then N2_11
            d_gam1_dU = 1.0/(1.0 - norm_f_2);
            break;
        case INDEX_N2_12: // Then N2_12
            d_gam1_dU = 0.0;
            break;
        case INDEX_N2_22: // Then N2_22
            d_gam1_dU = 0.0;
            break;
        default:
            cout << "Ivalid value for index_U in dgam1_dU !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    if (fabs(1.0 - norm_f_2) < 1.0e-8) {
        d_gam1_dU = 0.0;
    }
    
    return d_gam1_dU;
}

//***************************************************************************************************
// This routine computes the derivatives of the eigenvalues of the covariance matrix with respect to 
// the conserved solution variables
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::dgam2_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double d_gam2_dU;
    long double norm_f_2;
    norm_f_2 = pow(N1_1,2) + pow(N1_2,2) + pow(N1_3,2);
    
    switch (index_U) {
        case INDEX_I0:
            d_gam2_dU = (pow(N1_2, 2) - gam2*(1.0 + norm_f_2))/(1.0 - norm_f_2);
            cout << "Fix this !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
        case INDEX_N1_1: // Then N1_1
            d_gam2_dU = 2.0*N1_1*gam2/(1.0 - norm_f_2);
            break;
        case INDEX_N1_2: // Then N1_2
            d_gam2_dU = -2.0*N1_2*(1.0 - gam2)/(1.0 - norm_f_2);
            break;
        case INDEX_N2_11: // Then N2_11
            d_gam2_dU = 0.0;
            break;
        case INDEX_N2_12: // Then N2_12
            d_gam2_dU = 0.0;
            break;
        case INDEX_N2_22: // Then N2_22
            d_gam2_dU = 1.0/(1.0 - norm_f_2);
            break;
        default:
            cout << "Ivalid value for index_U in dgam2_dU !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    if (fabs(1.0 - norm_f_2) < 1.0e-8) {
        d_gam2_dU = 0.0;
    }
    
    return d_gam2_dU;
}

long double N3_Non_Gray_M2_3D_RT_Cheby::dB_dU(const int &index_U, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double d_B_dU;
    long double norm_f_2;
    norm_f_2 = pow(N1_1,2) + pow(N1_2,2) + pow(N1_3,2);
    
    // B = N2_12 - N1_1 N1_2
    
    switch (index_U) {
        case INDEX_I0:
            d_B_dU = 0.0;
            cout << "Fix this !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
        case INDEX_N1_1: // Then N1_1
            d_B_dU = -N1_2;
            break;
        case INDEX_N1_2: // Then N1_2
            d_B_dU = -N1_1;
            break;
        case INDEX_N2_11: // Then N2_11
            d_B_dU = 0.0;
            break;
        case INDEX_N2_12: // Then N2_12
            d_B_dU = 1.0;
            break;
        case INDEX_N2_22: // Then N2_22
            d_B_dU = 0.0;
            break;
        default:
            cout << "Ivalid value for index_U in dB_dU !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    return d_B_dU;
}

//***********************************************************************
// This routine computes the interpolative-based approximations of the 
// third-order closing fluxes, N3_ijk, in the frame where the covariance 
// matrix is diagonal
//***********************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby::set_Closure_RT() {
    Jacobian_M2.setup_Frame_Rotation();
    
    long double ratio_E, N1_1, N1_2, N1_3, gam1, gam2;
    long double L_I0_star_N3;
    
    N1_1 = Jacobian_M2.fx_RT();
    N1_2 = Jacobian_M2.fy_RT();
    N1_3 = 0.0;
    gam1 = Jacobian_M2.gam_1();
    gam2 = Jacobian_M2.gam_2();
    
    switch(Problem_Type) {
        case GRAY:
            ratio_E = 0.0;
            break;
        case NON_GRAY:
            L_I0_star_N3 = Evaluate_Length_Scale_N3(N1_1, N1_2, N1_3, gam1, gam2);
            ratio_E = Exponential_Mapping(Jacobian_M2.e_RT(), L_I0_star_N3);
            break;
    }
    
//     if (flag_compute_closing_fluxes) {
        Jacobian_M2.qxxx_RT_val = qxxx_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
        Jacobian_M2.qxxy_RT_val = qxxy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
        Jacobian_M2.qxyy_RT_val = qxyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
        Jacobian_M2.qyyy_RT_val = qyyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
        
        Jacobian_M2.qxxx_val = qxxx();
        Jacobian_M2.qxxy_val = qxxy();
        Jacobian_M2.qxyy_val = qxyy();
        Jacobian_M2.qyyy_val = qyyy();
//         
//         flag_compute_closing_fluxes = false;
//     }
    
//         cout << "E = " << Jacobian_M2.e_RT() << "   " << "fx = " << Jacobian_M2.fx_RT() << "   " << "fy = " << Jacobian_M2.fy_RT() << "   " << "gam_1 = " << Jacobian_M2.gam_1() << "   " << "gam_2 = " << Jacobian_M2.gam_2() << endl;
    
    
//     if (qxxx_RT_val != qxxx_RT_val || qxxy_RT_val != qxxy_RT_val || qxyy_RT_val != qxyy_RT_val || qyyy_RT_val != qyyy_RT_val) {
//         cout << "N1_1 = " << Vars_N3_RT.N1_1 << "   " << "N1_2 = " << Vars_N3_RT.N1_2 << "   " << "N1_3 = " << Vars_N3_RT.N1_3 << "   " << "gam1 = " << Vars_N3_RT.gam1 << "   " << "gam2 = " << Vars_N3_RT.gam2 << endl;
//         
//         cout << "x_SH = " << Vars_N3_RT.x_SH << "   " << "y_SH = " << Vars_N3_RT.y_SH << "   " << "z_SH = " << Vars_N3_RT.z_SH << endl;
//         
//         cout << "qxxx_RT_val = " << Jacobian_M2.qxxx_RT_val << "   " << "qxxy_RT_val = " << Jacobian_M2.qxxy_RT_val << "   " << "qxyy_RT_val = " << Jacobian_M2.qxyy_RT_val << "   " << "qyyy_RT_val = " << Jacobian_M2.qyyy_RT_val << endl;
//         
//         cout << endl;
//     }
//     
//     R_rot_mat_2D[0][0] = Cos;
//     R_rot_mat_2D[0][1] = -Sin;
//     R_rot_mat_2D[1][0] = Sin;
//     R_rot_mat_2D[1][1] = Cos;
}

//***********************************************************************
// This routine computes the interpolative-based approximations of the 
// third-order closing fluxes, N3_ijk, in the frame where the covariance 
// matrix is diagonal
//***********************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby::set_Closure_RT_Derivatives() {
    long double ratio_E, N1_1, N1_2, N1_3, gam1, gam2;
    long double L_I0_star_N3;
    set_Closure_RT();
    
    N1_1 = Jacobian_M2.fx_RT();
    N1_2 = Jacobian_M2.fy_RT();
    N1_3 = 0.0;
    gam1 = Jacobian_M2.gam_1();
    gam2 = Jacobian_M2.gam_2();
    
    switch(Problem_Type) {
        case GRAY:
            ratio_E = 0.0;
            break;
        case NON_GRAY:
            L_I0_star_N3 = Evaluate_Length_Scale_N3(N1_1, N1_2, N1_3, gam1, gam2);
            ratio_E = Exponential_Mapping(Jacobian_M2.e_RT(), L_I0_star_N3);
            break;
    }
    
    N1_1 = Jacobian_M2.fx_RT();
    N1_2 = Jacobian_M2.fy_RT();
    N1_3 = 0.0;
    gam1 = Jacobian_M2.gam_1();
    gam2 = Jacobian_M2.gam_2();
    
    Compute_dN3_ijk_dN2_12_Finite_Difference(Jacobian_M2.e_RT(), N1_1, N1_2, gam1, gam2);
    
    d_q_RT_dU(Jacobian_M2.dqxxx_RT, Jacobian_M2.dN3_111_RT, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
    d_q_RT_dU(Jacobian_M2.dqxxy_RT, Jacobian_M2.dN3_112_RT, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_112_ENTRY);
    d_q_RT_dU(Jacobian_M2.dqxyy_RT, Jacobian_M2.dN3_122_RT, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
    d_q_RT_dU(Jacobian_M2.dqyyy_RT, Jacobian_M2.dN3_222_RT, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_222_ENTRY);
    
    Jacobian_M2.Evaluate_d_q_RT_dU0();
    
//     cout << endl;
//     cout << "N1_1 = " << N1_1 << "   " << "N1_2 = " << N1_2 << "   " << "gam1 = " << gam1 << "   " << "gam2 = " << gam2 << endl;
//     
//     cout << "dqxxx_RT[0] = " << Jacobian_M2.dqxxx_RT[0] << "  " << "dqxxx_RT[1] = " << Jacobian_M2.dqxxx_RT[1] << "  " << "dqxxx_RT[2] = " << Jacobian_M2.dqxxx_RT[2] << "  " << "dqxxx_RT[3] = " << Jacobian_M2.dqxxx_RT[3] << "  " << "dqxxx_RT[4] = " << Jacobian_M2.dqxxx_RT[4] << "  " << "dqxxx_RT[5] = " << Jacobian_M2.dqxxx_RT[5] << endl;
//     
//     cout << "dqxxy_RT[0] = " << Jacobian_M2.dqxxy_RT[0] << "  " << "dqxxy_RT[1] = " << Jacobian_M2.dqxxy_RT[1] << "  " << "dqxxy_RT[2] = " << Jacobian_M2.dqxxy_RT[2] << "  " << "dqxxy_RT[3] = " << Jacobian_M2.dqxxy_RT[3] << "  " << "dqxxy_RT[4] = " << Jacobian_M2.dqxxy_RT[4] << "  " << "dqxxy_RT[5] = " << Jacobian_M2.dqxxy_RT[5] << endl;
//     
//     cout << "dqxyy_RT[0] = " << Jacobian_M2.dqxyy_RT[0] << "  " << "dqxyy_RT[1] = " << Jacobian_M2.dqxyy_RT[1] << "  " << "dqxyy_RT[2] = " << Jacobian_M2.dqxyy_RT[2] << "  " << "dqxyy_RT[3] = " << Jacobian_M2.dqxyy_RT[3] << "  " << "dqxyy_RT[4] = " << Jacobian_M2.dqxyy_RT[4] << "  " << "dqxyy_RT[5] = " << Jacobian_M2.dqxyy_RT[5] << endl;
//     cout << endl;
//     
//     cout << "qxxx_RT_val = " << Jacobian_M2.qxxx_RT_val << "  " << "qxxy_RT_val = " << Jacobian_M2.qxxy_RT_val << "  " << "qyyy_RT_val = " << Jacobian_M2.qyyy_RT_val << endl;
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Compute_dN3_ijk_dN2_12_Finite_Difference(const long double &I0_star_val, const long double &N1_1, const long double &N1_2, const long double &gam1, const long double &gam2) {
    long double epsilon = 1.0e-8;
    double N1_1_temp, N1_2_temp, N2_11, N2_12, N2_22;
    
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_L(N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_R(N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_Cheby(N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    
    Copy_to(&N3_M2_3D_RT_L);
    Copy_to(&N3_M2_3D_RT_R);
    Copy_to(&N3_M2_3D_RT_Cheby);
    
    N3_M2_3D_RT_Cheby.Jacobian_M2.Compute_W_array(I0_star_val, N1_1, N1_2, gam1, gam2);
//     N3_M2_3D_RT_Cheby.set_Closure_RT_Derivatives();
//     N3_M2_3D_RT_Cheby.Jacobian_M2.set_Closure_Derivatives_Original_Basis();
    
    N1_1_temp = N3_M2_3D_RT_Cheby.Jacobian_M2.W_array[1];
    N1_2_temp = N3_M2_3D_RT_Cheby.Jacobian_M2.W_array[2];
    N2_11 = N3_M2_3D_RT_Cheby.Jacobian_M2.W_array[3];
    N2_12 = N3_M2_3D_RT_Cheby.Jacobian_M2.W_array[4];
    N2_22 = N3_M2_3D_RT_Cheby.Jacobian_M2.W_array[5];
    
    // For N2_12
    N3_M2_3D_RT_L.Jacobian_M2.Set_W_array(I0_star_val, N1_1_temp, N1_2_temp, N2_11, N2_12 - epsilon, N2_22);
    N3_M2_3D_RT_R.Jacobian_M2.Set_W_array(I0_star_val, N1_1_temp, N1_2_temp, N2_11, N2_12 + epsilon, N2_22);
    N3_M2_3D_RT_L.set_Closure_RT();
    N3_M2_3D_RT_R.set_Closure_RT();
    
    Jacobian_M2.dN3_111_RT.dB = (N3_M2_3D_RT_R.Jacobian_M2.qxxx_val - N3_M2_3D_RT_L.Jacobian_M2.qxxx_val)/(2.0*epsilon);
    Jacobian_M2.dN3_112_RT.dB = (N3_M2_3D_RT_R.Jacobian_M2.qxxy_val - N3_M2_3D_RT_L.Jacobian_M2.qxxy_val)/(2.0*epsilon);
    Jacobian_M2.dN3_122_RT.dB = (N3_M2_3D_RT_R.Jacobian_M2.qxyy_val - N3_M2_3D_RT_L.Jacobian_M2.qxyy_val)/(2.0*epsilon);
    Jacobian_M2.dN3_222_RT.dB = (N3_M2_3D_RT_R.Jacobian_M2.qyyy_val - N3_M2_3D_RT_L.Jacobian_M2.qyyy_val)/(2.0*epsilon);
    
    Jacobian_M2.dN3_111_RT.dB = 0.0;
    Jacobian_M2.dN3_112_RT.dB = 0.0;
    Jacobian_M2.dN3_122_RT.dB = 0.0;
    Jacobian_M2.dN3_222_RT.dB = 0.0;
    
    // cout << "N2_12 = " << N2_12 << "  " << "dN3_111_RT.dB = " << Jacobian_M2.dN3_111_RT.dB << "  " << "dN3_112_RT.dB = " << Jacobian_M2.dN3_112_RT.dB << "  " << "dN3_122_RT.dB = " << Jacobian_M2.dN3_122_RT.dB << "  " << "dN3_222_RT.dB = " << Jacobian_M2.dN3_222_RT.dB << endl;
}

void N3_Non_Gray_M2_3D_RT_Cheby::d_q_RT_dU(long double *d_N3_RT_dU, record_d_N3_ijk_RT &dN3_ijk_RT, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {
    long double L_I0_star_N3, I0_star_val;
    long double dratioI0_d_L_I0_Star;
    record_d_L_I0_star_N3_ijk_RT d_r_I0_star_N3;
    
    switch (Problem_Type) {
        case GRAY:
            
            break;
        case NON_GRAY:
            switch (closure_type) {
                case MOMENT_CLOSURE_PROJECTION:
                    
                    break;
                case MOMENT_CLOSURE_M2:
                    L_I0_star_N3 = Evaluate_Length_Scale_N3(N1_1, N1_2, N1_3, gam1, gam2);
                    I0_star_val = Inverse_Exponential_Mapping(ratio_E, L_I0_star_N3);
                    
                    d_r_I0_star_N3.dI0_star = Diff_Exponential_Mapping(I0_star_val, L_I0_star_N3, INDEX_IO_STAR);
                    
                    dratioI0_d_L_I0_Star = Diff_Exponential_Mapping(I0_star_val, L_I0_star_N3, INDEX_LENGTH_SCALE_IO_STAR);
                    
                    d_r_I0_star_N3.d_norm_f_2 = Evaluate_diff_Length_Scale_N3_ijk(INDEX_NORM_f_2, N1_1, N1_2, N1_3, gam1, gam2);
                    d_r_I0_star_N3.d_norm_f_2 *= dratioI0_d_L_I0_Star;
                    
                    d_r_I0_star_N3.d_x2 = Evaluate_diff_Length_Scale_N3_ijk(INDEX_x_SH_2, N1_1, N1_2, N1_3, gam1, gam2);
                    d_r_I0_star_N3.d_x2 *= dratioI0_d_L_I0_Star;
                    
                    d_r_I0_star_N3.d_y2 = Evaluate_diff_Length_Scale_N3_ijk(INDEX_y_SH_2, N1_1, N1_2, N1_3, gam1, gam2);
                    d_r_I0_star_N3.d_y2 *= dratioI0_d_L_I0_Star;
                    
                    d_r_I0_star_N3.dgam1 = Evaluate_diff_Length_Scale_N3_ijk(INDEX_GAM1, N1_1, N1_2, N1_3, gam1, gam2);
                    d_r_I0_star_N3.dgam1 *= dratioI0_d_L_I0_Star;
                    
                    d_r_I0_star_N3.dgam2 = Evaluate_diff_Length_Scale_N3_ijk(INDEX_GAM2, N1_1, N1_2, N1_3, gam1, gam2);
                    d_r_I0_star_N3.dgam2 *= dratioI0_d_L_I0_Star;
                    break;
                default:
                    cout << "Incorrect value for Closure Type" << endl;
                    exit(0);
                    break;
            }
            break;
        default:
            cout << "Problem type not specified !!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
    
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            dN3_111_2D_dU(d_N3_RT_dU, dN3_ijk_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
            break;
        case N3_112_ENTRY:
            dN3_112_2D_dU(d_N3_RT_dU, dN3_ijk_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
            break;
        case N3_122_ENTRY:
            dN3_122_2D_dU(d_N3_RT_dU, dN3_ijk_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
            break;
        case N3_222_ENTRY:
            dN3_222_2D_dU(d_N3_RT_dU, dN3_ijk_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
            break;
        default:
            cout << "Invalid value for CLOSING_FLUX_INDEX in d_q_RT_dU !!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
}

void N3_Non_Gray_M2_3D_RT_Cheby::Setup_Flux_Jacobian_Matrix() {
    int NVARS = STATIC_NUM_VAR_RT;
    long double c = 1.0; // SPEED_OF_LIGHT;
    
    // set_Closure_RT_Derivatives();
    // Jacobian_M2.set_Closure_Derivatives_Original_Basis();
    
//     set_Closure_RT();
//     Compute_Closure_Derivatives_Original_Basis_Adept(m_values, io_band);
    
    for (int i = 0; i < NVARS; i++) {
        for (int j = 0; j < NVARS; j++) {
            Jacobian_M2.dFdU[i*NVARS+j] = 0.0;
        }
    }
    
// //     cout << endl;
//     
//     cout << "E = " << Jacobian_M2.e_RT() << "   " << "fx = " << Jacobian_M2.fx_RT() << "   " << "fy = " << Jacobian_M2.fy_RT() << "   " << "gam_1 = " << Jacobian_M2.gam_1() << "   " << "gam_2 = " << Jacobian_M2.gam_2() << endl;
//     
//     cout << "qxxx_RT_val = " << Jacobian_M2.qxxx_RT_val << "  " << "qxxy_RT_val = " << Jacobian_M2.qxxy_RT_val << "  " << "qxyy_RT_val = " << Jacobian_M2.qxyy_RT_val << "  " << "qyyy_RT_val = " << Jacobian_M2.qyyy_RT_val << endl;
//     
//     cout << "qxxx_val = " << Jacobian_M2.qxxx_val << "  " << "qxxy_val = " << Jacobian_M2.qxxy_val << "  " << "qxyy_val = " << Jacobian_M2.qxyy_val << "  " << "qyyy_val = " << Jacobian_M2.qyyy_val << endl;
//     
//     cout << "dqxxx_dU[0] = " << Jacobian_M2.dqxxx_dU[0] << "  " << "dqxxx_dU[1] = " << Jacobian_M2.dqxxx_dU[1] << "  " << "dqxxx_dU[2] = " << Jacobian_M2.dqxxx_dU[2] << "  " << "dqxxx_dU[3] = " << Jacobian_M2.dqxxx_dU[3] << "  " << "dqxxx_dU[4] = " << Jacobian_M2.dqxxx_dU[4] << "  " << "dqxxx_dU[5] = " << Jacobian_M2.dqxxx_dU[5] << endl;
//     
//     cout << "dqxxy_dU[0] = " << Jacobian_M2.dqxxy_dU[0] << "  " << "dqxxy_dU[1] = " << Jacobian_M2.dqxxy_dU[1] << "  " << "dqxxy_dU[2] = " << Jacobian_M2.dqxxy_dU[2] << "  " << "dqxxy_dU[3] = " << Jacobian_M2.dqxxy_dU[3] << "  " << "dqxxy_dU[4] = " << Jacobian_M2.dqxxy_dU[4] << "  " << "dqxxy_dU[5] = " << Jacobian_M2.dqxxy_dU[5] << endl;
//     
//     cout << "dqxyy_dU[0] = " << Jacobian_M2.dqxyy_dU[0] << "  " << "dqxyy_dU[1] = " << Jacobian_M2.dqxyy_dU[1] << "  " << "dqxyy_dU[2] = " << Jacobian_M2.dqxyy_dU[2] << "  " << "dqxyy_dU[3] = " << Jacobian_M2.dqxyy_dU[3] << "  " << "dqxyy_dU[4] = " << Jacobian_M2.dqxyy_dU[4] << "  " << "dqxyy_dU[5] = " << Jacobian_M2.dqxyy_dU[5] << endl;
//     cout << endl;
    
    Jacobian_M2.dFdU[0*NVARS+1] = c;
    Jacobian_M2.dFdU[1*NVARS+3] = c;
    Jacobian_M2.dFdU[2*NVARS+4] = c;
    Jacobian_M2.dFdU[3*NVARS+0] = c*Jacobian_M2.dqxxx_dU[0];
    Jacobian_M2.dFdU[3*NVARS+1] = c*Jacobian_M2.dqxxx_dU[1];
    Jacobian_M2.dFdU[3*NVARS+2] = c*Jacobian_M2.dqxxx_dU[2];
    Jacobian_M2.dFdU[3*NVARS+3] = c*Jacobian_M2.dqxxx_dU[3];
    Jacobian_M2.dFdU[3*NVARS+4] = c*Jacobian_M2.dqxxx_dU[4];
    Jacobian_M2.dFdU[3*NVARS+5] = c*Jacobian_M2.dqxxx_dU[5];
    Jacobian_M2.dFdU[4*NVARS+0] = c*Jacobian_M2.dqxxy_dU[0];
    Jacobian_M2.dFdU[4*NVARS+1] = c*Jacobian_M2.dqxxy_dU[1];
    Jacobian_M2.dFdU[4*NVARS+2] = c*Jacobian_M2.dqxxy_dU[2];
    Jacobian_M2.dFdU[4*NVARS+3] = c*Jacobian_M2.dqxxy_dU[3];
    Jacobian_M2.dFdU[4*NVARS+4] = c*Jacobian_M2.dqxxy_dU[4];
    Jacobian_M2.dFdU[4*NVARS+5] = c*Jacobian_M2.dqxxy_dU[5];
    Jacobian_M2.dFdU[5*NVARS+0] = c*Jacobian_M2.dqxyy_dU[0];
    Jacobian_M2.dFdU[5*NVARS+1] = c*Jacobian_M2.dqxyy_dU[1];
    Jacobian_M2.dFdU[5*NVARS+2] = c*Jacobian_M2.dqxyy_dU[2];
    Jacobian_M2.dFdU[5*NVARS+3] = c*Jacobian_M2.dqxyy_dU[3];
    Jacobian_M2.dFdU[5*NVARS+4] = c*Jacobian_M2.dqxyy_dU[4];
    Jacobian_M2.dFdU[5*NVARS+5] = c*Jacobian_M2.dqxyy_dU[5];
    
//     Jacobian_M2.dFdU[0*NVARS+1] = c;
//     Jacobian_M2.dFdU[1*NVARS+3] = c;
//     Jacobian_M2.dFdU[2*NVARS+4] = c;
//     Jacobian_M2.dFdU[3*NVARS+0] = c*Jacobian_M2.dqxxx_RT[0];
//     Jacobian_M2.dFdU[3*NVARS+1] = c*Jacobian_M2.dqxxx_RT[1];
//     Jacobian_M2.dFdU[3*NVARS+2] = c*Jacobian_M2.dqxxx_RT[2];
//     Jacobian_M2.dFdU[3*NVARS+3] = c*Jacobian_M2.dqxxx_RT[3];
//     Jacobian_M2.dFdU[3*NVARS+4] = c*Jacobian_M2.dqxxx_RT[4];
//     Jacobian_M2.dFdU[3*NVARS+5] = c*Jacobian_M2.dqxxx_RT[5];
//     Jacobian_M2.dFdU[4*NVARS+0] = c*Jacobian_M2.dqxxy_RT[0];
//     Jacobian_M2.dFdU[4*NVARS+1] = c*Jacobian_M2.dqxxy_RT[1];
//     Jacobian_M2.dFdU[4*NVARS+2] = c*Jacobian_M2.dqxxy_RT[2];
//     Jacobian_M2.dFdU[4*NVARS+3] = c*Jacobian_M2.dqxxy_RT[3];
//     Jacobian_M2.dFdU[4*NVARS+4] = c*Jacobian_M2.dqxxy_RT[4];
//     Jacobian_M2.dFdU[4*NVARS+5] = c*Jacobian_M2.dqxxy_RT[5];
//     Jacobian_M2.dFdU[5*NVARS+0] = c*Jacobian_M2.dqxyy_RT[0];
//     Jacobian_M2.dFdU[5*NVARS+1] = c*Jacobian_M2.dqxyy_RT[1];
//     Jacobian_M2.dFdU[5*NVARS+2] = c*Jacobian_M2.dqxyy_RT[2];
//     Jacobian_M2.dFdU[5*NVARS+3] = c*Jacobian_M2.dqxyy_RT[3];
//     Jacobian_M2.dFdU[5*NVARS+4] = c*Jacobian_M2.dqxyy_RT[4];
//     Jacobian_M2.dFdU[5*NVARS+5] = c*Jacobian_M2.dqxyy_RT[5];
}

// void N3_Non_Gray_M2_3D_RT_Cheby :: Compute_Closure_Derivatives_Original_Basis_Adept() {
//     int NUM_CLOSING_FLUXES = 4;
//     double jac[NUM_CLOSING_FLUXES*STATIC_NUM_VAR_RT];
//     Stack stack; // Where the derivative information is stored
//     
//     N3_Non_Gray_M2_3D_RT_Cheby<adouble> N3_RT_Adept;
//     // Setup static variables
//     N3_RT_Adept.v_freq = v_freq;
//     
//      // Vector of active input variables
//     N3_RT_Adept.Jacobian_M2.W_array[0].set_value(m_values[io_band]);   // Fill vector x
//     N3_RT_Adept.Jacobian_M2.W_array[1].set_value(m_values[io_band + 1]);   // Fill vector x
//     N3_RT_Adept.Jacobian_M2.W_array[2].set_value(m_values[io_band + 2]);   // Fill vector x
//     N3_RT_Adept.Jacobian_M2.W_array[3].set_value(m_values[io_band + 3]);   // Fill vector x
//     N3_RT_Adept.Jacobian_M2.W_array[4].set_value(m_values[io_band + 4]);   // Fill vector x
//     N3_RT_Adept.Jacobian_M2.W_array[5].set_value(m_values[io_band + 5]);   // Fill vector x
//     
//     stack.new_recording(); // Start recording
//     
//     vector<adouble> N3_ijk(NUM_CLOSING_FLUXES); // Create vector of active output variables
//     
//     N3_RT_Adept.set_Closure_RT();
//     N3_ijk[0] = N3_RT_Adept.Jacobian_M2.qxxx_val;    // Compute dependent variable
//     N3_ijk[1] = N3_RT_Adept.Jacobian_M2.qxxy_val;    // Compute dependent variable
//     N3_ijk[2] = N3_RT_Adept.Jacobian_M2.qxyy_val;    // Compute dependent variable
//     N3_ijk[3] = N3_RT_Adept.Jacobian_M2.qyyy_val;    // Compute dependent variable
//     
//     stack.independent(&N3_RT_Adept.Jacobian_M2.W_array[0], STATIC_NUM_VAR_RT); // Identify independent variables
//     stack.dependent(&N3_ijk[0], NUM_CLOSING_FLUXES); // Identify dependent variables
//     stack.jacobian(jac); // Compute & store Jacobian in jac
//     
//     for (int i = 0; i < STATIC_NUM_VAR_RT; i++) {
//         Jacobian_M2.dqxxx_dU[i] = jac[i*NUM_CLOSING_FLUXES + 0];
//         Jacobian_M2.dqxxy_dU[i] = jac[i*NUM_CLOSING_FLUXES + 1];
//         Jacobian_M2.dqxyy_dU[i] = jac[i*NUM_CLOSING_FLUXES + 2];
//         Jacobian_M2.dqyyy_dU[i] = jac[i*NUM_CLOSING_FLUXES + 3];
//     }
//     
//     Jacobian_M2.Evaluate_d_q_Orig_Basis_dU0_Adept();
//     
// //     cout << "dqxxx_dU[0] = " << Jacobian_M2.dqxxx_dU[0] << "  " << "dqxxx_dU[1] = " << Jacobian_M2.dqxxx_dU[1] << "  " << "dqxxx_dU[2] = " << Jacobian_M2.dqxxx_dU[2] << "  " << "dqxxx_dU[3] = " << Jacobian_M2.dqxxx_dU[3] << "  " << "dqxxx_dU[4] = " << Jacobian_M2.dqxxx_dU[4] << "  " << "dqxxx_dU[5] = " << Jacobian_M2.dqxxx_dU[5] << endl;
// //     
// //     cout << "dqxxy_dU[0] = " << Jacobian_M2.dqxxy_dU[0] << "  " << "dqxxy_dU[1] = " << Jacobian_M2.dqxxy_dU[1] << "  " << "dqxxy_dU[2] = " << Jacobian_M2.dqxxy_dU[2] << "  " << "dqxxy_dU[3] = " << Jacobian_M2.dqxxy_dU[3] << "  " << "dqxxy_dU[4] = " << Jacobian_M2.dqxxy_dU[4] << "  " << "dqxxy_dU[5] = " << Jacobian_M2.dqxxy_dU[5] << endl;
// //     
// //     cout << "dqxyy_dU[0] = " << Jacobian_M2.dqxyy_dU[0] << "  " << "dqxyy_dU[1] = " << Jacobian_M2.dqxyy_dU[1] << "  " << "dqxyy_dU[2] = " << Jacobian_M2.dqxyy_dU[2] << "  " << "dqxyy_dU[3] = " << Jacobian_M2.dqxyy_dU[3] << "  " << "dqxyy_dU[4] = " << Jacobian_M2.dqxyy_dU[4] << "  " << "dqxyy_dU[5] = " << Jacobian_M2.dqxyy_dU[5] << endl;
// //     cout << endl;
// //     
// //     cout << "dqyyy_dU[0] = " << Jacobian_M2.dqyyy_dU[0] << "  " << "dqyyy_dU[1] = " << Jacobian_M2.dqyyy_dU[1] << "  " << "dqyyy_dU[2] = " << Jacobian_M2.dqyyy_dU[2] << "  " << "dqyyy_dU[3] = " << Jacobian_M2.dqyyy_dU[3] << "  " << "dqyyy_dU[4] = " << Jacobian_M2.dqyyy_dU[4] << "  " << "dqyyy_dU[5] = " << Jacobian_M2.dqyyy_dU[5] << endl;
// //     cout << endl;
// }

//***************************************************************************************************
// This routine computes the third-order closing flux, N3_111, in the frame where the covariance 
// matrix is diagonal, based on the proposed interpolative-based approximation of the closing fluxes
// in such a frame
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::qxxx_RT(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double q_xxx;
    
    q_xxx = Evaluate_N3_111(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
    
    return q_xxx;
}

//***************************************************************************************************
// This routine computes the third-order closing flux, N3_112, in the frame where the covariance 
// matrix is diagonal, based on the proposed interpolative-based approximation of the closing fluxes
// in such a frame
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::qxxy_RT(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double q_xxy;
    
    q_xxy = Evaluate_N3_112(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
    
    return q_xxy;
}

//***************************************************************************************************
// This routine computes the third-order closing flux, N3_122, in the frame where the covariance 
// matrix is diagonal, based on the proposed interpolative-based approximation of the closing fluxes
// in such a frame
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::qxyy_RT(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double q_xyy;
    
    q_xyy = Evaluate_N3_122(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
    
    return q_xyy;
}

//***************************************************************************************************
// This routine computes the third-order closing flux, N3_222, in the frame where the covariance 
// matrix is diagonal, based on the proposed interpolative-based approximation of the closing fluxes
// in such a frame
//***************************************************************************************************
long double N3_Non_Gray_M2_3D_RT_Cheby::qyyy_RT(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
    long double q_yyy;
    
    q_yyy = Evaluate_N3_222(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
    
    return q_yyy;
}

long double N3_Non_Gray_M2_3D_RT_Cheby::qxxx() {
    long double q_xxx;
    
    q_xxx = Jacobian_M2.Rotate_N3(Jacobian_M2.rotation, N3_111_ENTRY);
    
    return q_xxx;
}

long double N3_Non_Gray_M2_3D_RT_Cheby::qxxy() {
    long double q_xxy;
    
    q_xxy = Jacobian_M2.Rotate_N3(Jacobian_M2.rotation, N3_112_ENTRY);
    
    return q_xxy;
}

long double N3_Non_Gray_M2_3D_RT_Cheby::qxyy() {
    long double q_xyy;
    
    q_xyy = Jacobian_M2.Rotate_N3(Jacobian_M2.rotation, N3_122_ENTRY);
    
    return q_xyy;
}

long double N3_Non_Gray_M2_3D_RT_Cheby::qyyy() {
    long double q_yyy;
    
    q_yyy = Jacobian_M2.Rotate_N3(Jacobian_M2.rotation, N3_222_ENTRY);
    
    return q_yyy;
}

void N3_Non_Gray_M2_3D_RT_Cheby::Test_Derivatives_N3_ijk_RT(record_d_N3_ijk_RT &dN3_ijk_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &CLOSING_FLUX_INDEX) {
    record_d_N3_ijk_RT dN3_ijk_RT_temp;
    long double diff_N3_ijk_num, diff_N3_ijk_analy;
    long double N3_ijk_R, N3_ijk_L;
    long double epsilon = 1.0e-6;
    
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            dN3_ijk_2D_dU(dN3_ijk_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
            break;
        case N3_112_ENTRY:
            dN3_ijk_2D_dU(dN3_ijk_RT_temp, d_r_I0_star_N3, ratio_E, N1_2, -N1_1, N1_3, gam2, gam1, N3_122_ENTRY);
            dN3_ijk_RT.dI0 = dN3_ijk_RT_temp.dI0;
            dN3_ijk_RT.dN1_1 = -dN3_ijk_RT_temp.dN1_2;
            dN3_ijk_RT.dN1_2 = dN3_ijk_RT_temp.dN1_1;
            dN3_ijk_RT.dgam1 = dN3_ijk_RT_temp.dgam2;
            dN3_ijk_RT.dgam2 = dN3_ijk_RT_temp.dgam1;
            break;
        case N3_122_ENTRY:
            dN3_ijk_2D_dU(dN3_ijk_RT, d_r_I0_star_N3, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
            break;
        case N3_222_ENTRY:
            dN3_ijk_2D_dU(dN3_ijk_RT_temp, d_r_I0_star_N3, ratio_E, N1_2, -N1_1, N1_3, gam2, gam1, N3_111_ENTRY);
            dN3_ijk_RT.dI0 = dN3_ijk_RT_temp.dI0;
            dN3_ijk_RT.dN1_1 = -dN3_ijk_RT_temp.dN1_2;
            dN3_ijk_RT.dN1_2 = dN3_ijk_RT_temp.dN1_1;
            dN3_ijk_RT.dgam1 = dN3_ijk_RT_temp.dgam2;
            dN3_ijk_RT.dgam2 = dN3_ijk_RT_temp.dgam1;
            break;
    }
    
    // For N1_1
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            N3_ijk_L = qxxx_RT(ratio_E, N1_1 - epsilon, N1_2, N1_3, gam1, gam2);
            N3_ijk_R = qxxx_RT(ratio_E, N1_1 + epsilon, N1_2, N1_3, gam1, gam2);
            break;
        case N3_112_ENTRY:
            N3_ijk_L = qxxy_RT(ratio_E, N1_1 - epsilon, N1_2, N1_3, gam1, gam2);
            N3_ijk_R = qxxy_RT(ratio_E, N1_1 + epsilon, N1_2, N1_3, gam1, gam2);
            break;
        case N3_122_ENTRY:
            N3_ijk_L = qxyy_RT(ratio_E, N1_1 - epsilon, N1_2, N1_3, gam1, gam2);
            N3_ijk_R = qxyy_RT(ratio_E, N1_1 + epsilon, N1_2, N1_3, gam1, gam2);
            break;
        case N3_222_ENTRY:
            N3_ijk_L = qyyy_RT(ratio_E, N1_1 - epsilon, N1_2, N1_3, gam1, gam2);
            N3_ijk_R = qyyy_RT(ratio_E, N1_1 + epsilon, N1_2, N1_3, gam1, gam2);
            break;
    }
    // N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1 - epsilon, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    // N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1 + epsilon, N1_2, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    diff_N3_ijk_num = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    diff_N3_ijk_analy = dN3_ijk_RT.dN1_1;
    
    cout << "For N1_1, we have diff_N3_ijk_num = " << diff_N3_ijk_num << "  " << "diff_N3_ijk_analy = " << diff_N3_ijk_analy << "  " << "diff_N3_ijk_num - diff_N3_ijk_analy = " << diff_N3_ijk_num - diff_N3_ijk_analy << endl;
    
    if (fabs(diff_N3_ijk_num - diff_N3_ijk_analy) > 1.0e-8) {
        exit(0);
    }
    
    // For N1_2
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            N3_ijk_L = qxxx_RT(ratio_E, N1_1, N1_2 - epsilon, N1_3, gam1, gam2);
            N3_ijk_R = qxxx_RT(ratio_E, N1_1, N1_2 + epsilon, N1_3, gam1, gam2);
            break;
        case N3_112_ENTRY:
            N3_ijk_L = qxxy_RT(ratio_E, N1_1, N1_2 - epsilon, N1_3, gam1, gam2);
            N3_ijk_R = qxxy_RT(ratio_E, N1_1, N1_2 + epsilon, N1_3, gam1, gam2);
            break;
        case N3_122_ENTRY:
            N3_ijk_L = qxyy_RT(ratio_E, N1_1, N1_2 - epsilon, N1_3, gam1, gam2);
            N3_ijk_R = qxyy_RT(ratio_E, N1_1, N1_2 + epsilon, N1_3, gam1, gam2);
            break;
        case N3_222_ENTRY:
            N3_ijk_L = qyyy_RT(ratio_E, N1_1, N1_2 - epsilon, N1_3, gam1, gam2);
            N3_ijk_R = qyyy_RT(ratio_E, N1_1, N1_2 + epsilon, N1_3, gam1, gam2);
            break;
    }
    // N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1, N1_2 - epsilon, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    // N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1, N1_2 + epsilon, N1_3, gam1, gam2, CLOSING_FLUX_INDEX);
    diff_N3_ijk_num = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    diff_N3_ijk_analy = dN3_ijk_RT.dN1_2;
    
    cout << "For N1_2, we have diff_N3_ijk_num = " << diff_N3_ijk_num << "  " << "diff_N3_ijk_analy = " << diff_N3_ijk_analy << "  " << "diff_N3_ijk_num - diff_N3_ijk_analy = " << diff_N3_ijk_num - diff_N3_ijk_analy << endl;
    
    if (fabs(diff_N3_ijk_num - diff_N3_ijk_analy) > 1.0e-8) {
        exit(0);
    }
    
    // For gam1
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            N3_ijk_L = qxxx_RT(ratio_E, N1_1, N1_2, N1_3, gam1 - epsilon, gam2);
            N3_ijk_R = qxxx_RT(ratio_E, N1_1, N1_2, N1_3, gam1 + epsilon, gam2);
            break;
        case N3_112_ENTRY:
            N3_ijk_L = qxxy_RT(ratio_E, N1_1, N1_2, N1_3, gam1 - epsilon, gam2);
            N3_ijk_R = qxxy_RT(ratio_E, N1_1, N1_2, N1_3, gam1 + epsilon, gam2);
            break;
        case N3_122_ENTRY:
            N3_ijk_L = qxyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1 - epsilon, gam2);
            N3_ijk_R = qxyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1 + epsilon, gam2);
            break;
        case N3_222_ENTRY:
            N3_ijk_L = qyyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1 - epsilon, gam2);
            N3_ijk_R = qyyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1 + epsilon, gam2);
            break;
    }
    // N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1 - epsilon, gam2, CLOSING_FLUX_INDEX);
    // N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1 + epsilon, gam2, CLOSING_FLUX_INDEX);
    diff_N3_ijk_num = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    diff_N3_ijk_analy = dN3_ijk_RT.dgam1;
    
    cout << "For gam1, we have diff_N3_ijk_num = " << diff_N3_ijk_num << "  " << "diff_N3_ijk_analy = " << diff_N3_ijk_analy << "  " << "diff_N3_ijk_num - diff_N3_ijk_analy = " << diff_N3_ijk_num - diff_N3_ijk_analy << endl;
    
    if (fabs(diff_N3_ijk_num - diff_N3_ijk_analy) > 1.0e-8) {
        exit(0);
    }
    
    // For gam2
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            N3_ijk_L = qxxx_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 - epsilon);
            N3_ijk_R = qxxx_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 + epsilon);
            break;
        case N3_112_ENTRY:
            N3_ijk_L = qxxy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 - epsilon);
            N3_ijk_R = qxxy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 + epsilon);
            break;
        case N3_122_ENTRY:
            N3_ijk_L = qxyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 - epsilon);
            N3_ijk_R = qxyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 + epsilon);
            break;
        case N3_222_ENTRY:
            N3_ijk_L = qyyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 - epsilon);
            N3_ijk_R = qyyy_RT(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 + epsilon);
            break;
    }
    // N3_ijk_L = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 - epsilon, CLOSING_FLUX_INDEX);
    // N3_ijk_R = Evaluate_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2 + epsilon, CLOSING_FLUX_INDEX);
    diff_N3_ijk_num = (N3_ijk_R - N3_ijk_L)/(2.0*epsilon);
    diff_N3_ijk_analy = dN3_ijk_RT.dgam2;
    
    cout << "For gam2, we have diff_N3_ijk_num = " << diff_N3_ijk_num << "  " << "diff_N3_ijk_analy = " << diff_N3_ijk_analy << "  " << "diff_N3_ijk_num - diff_N3_ijk_analy = " << diff_N3_ijk_num - diff_N3_ijk_analy << endl;
    
    if (fabs(diff_N3_ijk_num - diff_N3_ijk_analy) > 1.0e-8) {
        exit(0);
    }
}

void N3_Non_Gray_M2_3D_RT_Cheby::Test_Derivatives_N3_ijk(record_d_N3_ijk_RT &dN3_ijk_RT, const record_d_L_I0_star_N3_ijk_RT &d_r_I0_star_N3, const long double &I0_star_val, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &N2_11, const long double &N2_12, const long double &N2_22) {
    record_d_N3_ijk_RT dN3_ijk_RT_temp;
    long double d_N3_111_dN1_1_num, d_N3_111_dN1_2_num, d_N3_111_dN2_11_num, d_N3_111_dN2_12_num, d_N3_111_dN2_22_num;
    long double d_N3_112_dN1_1_num, d_N3_112_dN1_2_num, d_N3_112_dN2_11_num, d_N3_112_dN2_12_num, d_N3_112_dN2_22_num;
    long double d_N3_122_dN1_1_num, d_N3_122_dN1_2_num, d_N3_122_dN2_11_num, d_N3_122_dN2_12_num, d_N3_122_dN2_22_num;
    long double d_N3_222_dN1_1_num, d_N3_222_dN1_2_num, d_N3_222_dN2_11_num, d_N3_222_dN2_12_num, d_N3_222_dN2_22_num;
    
    long double d_N3_111_dN1_1_analy, d_N3_111_dN1_2_analy, d_N3_111_dN2_11_analy, d_N3_111_dN2_12_analy, d_N3_111_dN2_22_analy;
    long double d_N3_112_dN1_1_analy, d_N3_112_dN1_2_analy, d_N3_112_dN2_11_analy, d_N3_112_dN2_12_analy, d_N3_112_dN2_22_analy;
    long double d_N3_122_dN1_1_analy, d_N3_122_dN1_2_analy, d_N3_122_dN2_11_analy, d_N3_122_dN2_12_analy, d_N3_122_dN2_22_analy;
    long double d_N3_222_dN1_1_analy, d_N3_222_dN1_2_analy, d_N3_222_dN2_11_analy, d_N3_222_dN2_12_analy, d_N3_222_dN2_22_analy;
    
    long double d_N3_111_dI0_num, d_N3_112_dI0_num, d_N3_122_dI0_num, d_N3_222_dI0_num;
    
    long double epsilon = 1.0e-6;
    
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_L(N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_R(N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_Cheby(N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    
    Copy_to(&N3_M2_3D_RT_L);
    Copy_to(&N3_M2_3D_RT_R);
    Copy_to(&N3_M2_3D_RT_Cheby);
    
    N3_M2_3D_RT_Cheby.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2, N2_11, N2_12, N2_22);
    N3_M2_3D_RT_Cheby.set_Closure_RT_Derivatives();
    N3_M2_3D_RT_Cheby.Jacobian_M2.set_Closure_Derivatives_Original_Basis();
    
    // For N1_1
    N3_M2_3D_RT_L.Jacobian_M2.Set_W_array(I0_star_val, N1_1 - epsilon, N1_2, N2_11, N2_12, N2_22);
    N3_M2_3D_RT_R.Jacobian_M2.Set_W_array(I0_star_val, N1_1 + epsilon, N1_2, N2_11, N2_12, N2_22);
    N3_M2_3D_RT_L.set_Closure_RT();
    N3_M2_3D_RT_R.set_Closure_RT();
    
    d_N3_111_dN1_1_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxx_val - N3_M2_3D_RT_L.Jacobian_M2.qxxx_val)/(2.0*epsilon);
    d_N3_112_dN1_1_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxy_val - N3_M2_3D_RT_L.Jacobian_M2.qxxy_val)/(2.0*epsilon);
    d_N3_122_dN1_1_num = (N3_M2_3D_RT_R.Jacobian_M2.qxyy_val - N3_M2_3D_RT_L.Jacobian_M2.qxyy_val)/(2.0*epsilon);
    d_N3_222_dN1_1_num = (N3_M2_3D_RT_R.Jacobian_M2.qyyy_val - N3_M2_3D_RT_L.Jacobian_M2.qyyy_val)/(2.0*epsilon);
    
    // For N1_2
    N3_M2_3D_RT_L.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2 - epsilon, N2_11, N2_12, N2_22);
    N3_M2_3D_RT_R.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2 + epsilon, N2_11, N2_12, N2_22);
    N3_M2_3D_RT_L.set_Closure_RT();
    N3_M2_3D_RT_R.set_Closure_RT();
    
    d_N3_111_dN1_2_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxx_val - N3_M2_3D_RT_L.Jacobian_M2.qxxx_val)/(2.0*epsilon);
    d_N3_112_dN1_2_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxy_val - N3_M2_3D_RT_L.Jacobian_M2.qxxy_val)/(2.0*epsilon);
    d_N3_122_dN1_2_num = (N3_M2_3D_RT_R.Jacobian_M2.qxyy_val - N3_M2_3D_RT_L.Jacobian_M2.qxyy_val)/(2.0*epsilon);
    d_N3_222_dN1_2_num = (N3_M2_3D_RT_R.Jacobian_M2.qyyy_val - N3_M2_3D_RT_L.Jacobian_M2.qyyy_val)/(2.0*epsilon);
    
    // For N2_11
    N3_M2_3D_RT_L.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2, N2_11 - epsilon, N2_12, N2_22);
    N3_M2_3D_RT_R.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2, N2_11 + epsilon, N2_12, N2_22);
    N3_M2_3D_RT_L.set_Closure_RT();
    N3_M2_3D_RT_R.set_Closure_RT();
    
    d_N3_111_dN2_11_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxx_val - N3_M2_3D_RT_L.Jacobian_M2.qxxx_val)/(2.0*epsilon);
    d_N3_112_dN2_11_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxy_val - N3_M2_3D_RT_L.Jacobian_M2.qxxy_val)/(2.0*epsilon);
    d_N3_122_dN2_11_num = (N3_M2_3D_RT_R.Jacobian_M2.qxyy_val - N3_M2_3D_RT_L.Jacobian_M2.qxyy_val)/(2.0*epsilon);
    d_N3_222_dN2_11_num = (N3_M2_3D_RT_R.Jacobian_M2.qyyy_val - N3_M2_3D_RT_L.Jacobian_M2.qyyy_val)/(2.0*epsilon);
    
    // For N2_12
    N3_M2_3D_RT_L.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2, N2_11, N2_12 - epsilon, N2_22);
    N3_M2_3D_RT_R.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2, N2_11, N2_12 + epsilon, N2_22);
    N3_M2_3D_RT_L.set_Closure_RT();
    N3_M2_3D_RT_R.set_Closure_RT();
    
    // cout << "d_cs Finite Diff = " << (N3_M2_3D_RT_R.Jacobian_M2.rotation.x - N3_M2_3D_RT_L.Jacobian_M2.rotation.x)/(2.0*epsilon) << "  " << "d_sn Finite Diff = " << (N3_M2_3D_RT_R.Jacobian_M2.rotation.y - N3_M2_3D_RT_L.Jacobian_M2.rotation.y)/(2.0*epsilon) << endl; 
    
    d_N3_111_dN2_12_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxx_val - N3_M2_3D_RT_L.Jacobian_M2.qxxx_val)/(2.0*epsilon);
    d_N3_112_dN2_12_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxy_val - N3_M2_3D_RT_L.Jacobian_M2.qxxy_val)/(2.0*epsilon);
    d_N3_122_dN2_12_num = (N3_M2_3D_RT_R.Jacobian_M2.qxyy_val - N3_M2_3D_RT_L.Jacobian_M2.qxyy_val)/(2.0*epsilon);
    d_N3_222_dN2_12_num = (N3_M2_3D_RT_R.Jacobian_M2.qyyy_val - N3_M2_3D_RT_L.Jacobian_M2.qyyy_val)/(2.0*epsilon);
    
    // For N2_22
    N3_M2_3D_RT_L.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2, N2_11, N2_12, N2_22 - epsilon);
    N3_M2_3D_RT_R.Jacobian_M2.Set_W_array(I0_star_val, N1_1, N1_2, N2_11, N2_12, N2_22 + epsilon);
    N3_M2_3D_RT_L.set_Closure_RT();
    N3_M2_3D_RT_R.set_Closure_RT();
    
    d_N3_111_dN2_22_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxx_val - N3_M2_3D_RT_L.Jacobian_M2.qxxx_val)/(2.0*epsilon);
    d_N3_112_dN2_22_num = (N3_M2_3D_RT_R.Jacobian_M2.qxxy_val - N3_M2_3D_RT_L.Jacobian_M2.qxxy_val)/(2.0*epsilon);
    d_N3_122_dN2_22_num = (N3_M2_3D_RT_R.Jacobian_M2.qxyy_val - N3_M2_3D_RT_L.Jacobian_M2.qxyy_val)/(2.0*epsilon);
    d_N3_222_dN2_22_num = (N3_M2_3D_RT_R.Jacobian_M2.qyyy_val - N3_M2_3D_RT_L.Jacobian_M2.qyyy_val)/(2.0*epsilon);
    
    d_N3_111_dI0_num = N3_M2_3D_RT_Cheby.Jacobian_M2.qxxx_val - N1_1*d_N3_111_dN1_1_num - N1_2*d_N3_111_dN1_2_num - N2_11*d_N3_111_dN2_11_num - N2_12*d_N3_111_dN2_12_num - N2_22*d_N3_111_dN2_22_num;
    d_N3_112_dI0_num = N3_M2_3D_RT_Cheby.Jacobian_M2.qxxy_val - N1_1*d_N3_112_dN1_1_num - N1_2*d_N3_112_dN1_2_num - N2_11*d_N3_112_dN2_11_num - N2_12*d_N3_112_dN2_12_num - N2_22*d_N3_112_dN2_22_num;
    d_N3_122_dI0_num = N3_M2_3D_RT_Cheby.Jacobian_M2.qxyy_val - N1_1*d_N3_122_dN1_1_num - N1_2*d_N3_122_dN1_2_num - N2_11*d_N3_122_dN2_11_num - N2_12*d_N3_122_dN2_12_num - N2_22*d_N3_122_dN2_22_num;
    d_N3_222_dI0_num = N3_M2_3D_RT_Cheby.Jacobian_M2.qyyy_val - N1_1*d_N3_222_dN1_1_num - N1_2*d_N3_222_dN1_2_num - N2_11*d_N3_222_dN2_11_num - N2_12*d_N3_222_dN2_12_num - N2_22*d_N3_222_dN2_22_num;
    
    cout << "Finite Differencing !!!!!!!!!!!!!!!!" << endl;
    cout << "dqxxx_dU[0] = " << d_N3_111_dI0_num << "  " << "dqxxx_dU[1] = " << d_N3_111_dN1_1_num << "  " << "dqxxx_dU[2] = " << d_N3_111_dN1_2_num << "  " << "dqxxx_dU[3] = " << d_N3_111_dN2_11_num << "  " << "dqxxx_dU[4] = " << d_N3_111_dN2_12_num << "  " << "dqxxx_dU[5] = " << d_N3_111_dN2_22_num << endl;
    
    cout << "dqxxy_dU[0] = " << d_N3_112_dI0_num << "  " << "dqxxy_dU[1] = " << d_N3_112_dN1_1_num << "  " << "dqxxy_dU[2] = " << d_N3_112_dN1_2_num << "  " << "dqxxy_dU[3] = " << d_N3_112_dN2_11_num << "  " << "dqxxy_dU[4] = " << d_N3_112_dN2_12_num << "  " << "dqxxy_dU[5] = " << d_N3_112_dN2_22_num << endl;
    
    cout << "dqxyy_dU[0] = " << d_N3_122_dI0_num << "  " << "dqxyy_dU[1] = " << d_N3_122_dN1_1_num << "  " << "dqxyy_dU[2] = " << d_N3_122_dN1_2_num << "  " << "dqxyy_dU[3] = " << d_N3_122_dN2_11_num << "  " << "dqxyy_dU[4] = " << d_N3_122_dN2_12_num << "  " << "dqxyy_dU[5] = " << d_N3_122_dN2_22_num << endl;
    
    cout << "dqyyy_dU[0] = " << d_N3_222_dI0_num << "  " << "dqyyy_dU[1] = " << d_N3_222_dN1_1_num << "  " << "dqyyy_dU[2] = " << d_N3_222_dN1_2_num << "  " << "dqyyy_dU[3] = " << d_N3_222_dN2_11_num << "  " << "dqyyy_dU[4] = " << d_N3_222_dN2_12_num << "  " << "dqyyy_dU[5] = " << d_N3_222_dN2_22_num << endl;
    
    cout << "Analytical !!!!!!!!!!!!!!!!" << endl;
//     switch (DOMAIN_TYPE) {
//         case R2_RT:
//             cout << "dqxxx_RT[0] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_RT[0] << "  " << "dqxxx_RT[1] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_RT[1] << "  " << "dqxxx_RT[2] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_RT[2] << "  " << "dqxxx_RT[3] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_RT[3] << "  " << "dqxxx_RT[4] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_RT[4] << "  " << "dqxxx_RT[5] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_RT[5] << endl;
//             
//             cout << "dqxxy_RT[0] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_RT[0] << "  " << "dqxxy_RT[1] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_RT[1] << "  " << "dqxxy_RT[2] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_RT[2] << "  " << "dqxxy_RT[3] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_RT[3] << "  " << "dqxxy_RT[4] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_RT[4] << "  " << "dqxxy_RT[5] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_RT[5] << endl;
//             
//             cout << "dqxyy_RT[0] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_RT[0] << "  " << "dqxyy_RT[1] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_RT[1] << "  " << "dqxyy_RT[2] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_RT[2] << "  " << "dqxyy_RT[3] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_RT[3] << "  " << "dqxyy_RT[4] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_RT[4] << "  " << "dqxyy_RT[5] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_RT[5] << endl;
//             
//             cout << "dqyyy_RT[0] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_RT[0] << "  " << "dqyyy_RT[1] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_RT[1] << "  " << "dqyyy_RT[2] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_RT[2] << "  " << "dqyyy_RT[3] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_RT[3] << "  " << "dqyyy_RT[4] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_RT[4] << "  " << "dqyyy_RT[5] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_RT[5] << endl;
//             break;
//         case R2_ORIGINAL:
            cout << "dqxxx_dU[0] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_dU[0] << "  " << "dqxxx_dU[1] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_dU[1] << "  " << "dqxxx_dU[2] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_dU[2] << "  " << "dqxxx_dU[3] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_dU[3] << "  " << "dqxxx_dU[4] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_dU[4] << "  " << "dqxxx_dU[5] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxx_dU[5] << endl;
            
            cout << "dqxxy_dU[0] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_dU[0] << "  " << "dqxxy_dU[1] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_dU[1] << "  " << "dqxxy_dU[2] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_dU[2] << "  " << "dqxxy_dU[3] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_dU[3] << "  " << "dqxxy_dU[4] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_dU[4] << "  " << "dqxxy_dU[5] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxxy_dU[5] << endl;
            
            cout << "dqxyy_dU[0] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_dU[0] << "  " << "dqxyy_dU[1] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_dU[1] << "  " << "dqxyy_dU[2] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_dU[2] << "  " << "dqxyy_dU[3] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_dU[3] << "  " << "dqxyy_dU[4] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_dU[4] << "  " << "dqxyy_dU[5] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqxyy_dU[5] << endl;
            
            cout << "dqyyy_dU[0] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_dU[0] << "  " << "dqyyy_dU[1] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_dU[1] << "  " << "dqyyy_dU[2] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_dU[2] << "  " << "dqyyy_dU[3] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_dU[3] << "  " << "dqyyy_dU[4] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_dU[4] << "  " << "dqyyy_dU[5] = " << N3_M2_3D_RT_Cheby.Jacobian_M2.dqyyy_dU[5] << endl;
//             break;
//         default:
//             cout << "Invalid choice for DOMAIN_TYPE!!!!!!!!!!!!!!!!!!!" << endl;
//             exit(0);
//             break;
//     }
}


// long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_dratio_E_dLength_Scale(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2) {
//     long double dratio_dL;
//     long double Length_Scale, Moment;
//     
//     Length_Scale = Evaluate_Length_Scale_N3_ijk(N1_1, N1_2, N1_3, gam1, gam2);
//     
//     Moment = Inverse_Mobius_Transformation(ratio_E, Length_Scale);
//     
//     // [-(Moment + Length_Scale) - (Moment - Length_Scale)]/(Moment + Length_Scale)^2
//     // (-2*Moment)/(Moment + Length_Scale)^2
//     
//     if (fabs(fabs(ratio_E) - 1) < 1.0e-6) {
//         dratio_dL = 0.0;
//     } else {
//         dratio_dL = -2.0*Moment/pow(Moment + Length_Scale, 2);
//     }
//     
//     if (dratio_dL != dratio_dL) {
//         cout << "dratio_dL = " << dratio_dL << "   " << "Length_Scale = " << Length_Scale << "   " << "Moment = " << Moment << "   " << "ratio_E = " << ratio_E << endl;
//     }
//     
//     return dratio_dL;
// }

// ******************************************************************************************
// This routine computes the weighting function, g_N3_111, g_N3_122, or g_N3_123, for the 
// interpolative-based approximation of the closing flux, N3_111
// ******************************************************************************************
// long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_g_N3_ijk_Least_Squares_L_I0_star(const long double &ratio_E, const int &VAR_NUM) {
//     long double N3_Fit_E;
//     int index_Coeffs;
//     
//     index_Coeffs = N_Points_E - 1;
//     
//     for (int i_fit_E = N_Points_E - 1; i_fit_E >= 0; i_fit_E--) {
//         switch (VAR_NUM) {
//             case N3_111_ENTRY:
//                 if (i_fit_E == N_Points_E - 1) {
//                     N3_Fit_E = Coefficient_Matrix_Fit_N3_111[index_Coeffs];
//                 } else {
//                     N3_Fit_E = Coefficient_Matrix_Fit_N3_111[index_Coeffs] + N3_Fit_E*ratio_E;
//                 }
//                 break;
//             case N3_122_ENTRY:
//                 if (i_fit_E == N_Points_E - 1) {
//                     N3_Fit_E = Coefficient_Matrix_Fit_N3_122[index_Coeffs];
//                 } else {
//                     N3_Fit_E = Coefficient_Matrix_Fit_N3_122[index_Coeffs] + N3_Fit_E*ratio_E;
//                 }
//                 break;
//             case N3_123_ENTRY:
//                 if (i_fit_E == N_Points_E - 1) {
//                     N3_Fit_E = Coefficient_Matrix_Fit_N3_123[index_Coeffs];
//                 } else {
//                     N3_Fit_E = Coefficient_Matrix_Fit_N3_123[index_Coeffs] + N3_Fit_E*ratio_E;
//                 }
//                 break;
//             default:
//                 cout << "Values for VAR_NUM not specified !!!!!!!!!!!!!!!!!!!!!" << endl;
//                 exit(0);
//                 break;
//         }
//         index_Coeffs--;
//     }
//     
//     return N3_Fit_E;
// }
// 
// long double N3_Non_Gray_M2_3D_RT_Cheby :: Evaluate_dg_N3_ijk_dratio_E(const long double &ratio_E, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2, const int &VAR_NUM) {
//     long double N3_Fit, N3_Fit_m, N3_Fit_f, N3_Fit_E, N3_Fit_gam1, N3_Fit_gam2, norm_f;
//     long double norm_f_2;
//     int index_Coeffs;
//     int N_Points_Gray_Cheby_without_Lebed = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
//     long double x, z, x2, z2;
//     long double gam3;
//     
//     index_Coeffs = N_Points_Gray_Cheby_without_Lebed*N_Coeffs_SH - 1;
//     
//     norm_f_2 = pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2);
//     norm_f = sqrt(norm_f_2);
//     
//     x = N1_1/norm_f;
//     z = N1_3/norm_f;
//     
//     if (norm_f < 1.0e-8) {
//         x = 0.0;
//         z = 1.0;
//     }
//     
//     x2 = x*x;
//     z2 = z*z;
//     
//     for (int l = Order_SH; l >= 0; l-=2) {
//         for (int m = Order_SH - l; m >= 0; m-=2) {
//             for (int i_fit_E = N_Points_E - 1; i_fit_E >= 0; i_fit_E--) {
//                 for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
//                     for (int i_fit_gam1 = N_Points_gam1 - 1; i_fit_gam1 >= 0; i_fit_gam1--) {
//                         for (int i_fit_gam2 = N_Points_gam1 - i_fit_gam1 - 1; i_fit_gam2 >= 0; i_fit_gam2--) {
//                             switch (VAR_NUM) {
//                                 case N3_111_ENTRY:
//                                     if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
//                                         N3_Fit_gam2 = Coefficient_Matrix_Fit_N3_123[index_Coeffs];
//                                     } else {
//                                         N3_Fit_gam2 = Coefficient_Matrix_Fit_N3_123[index_Coeffs] + N3_Fit_gam2*gam2;
//                                     }
//                                     break;
//                                 case N3_122_ENTRY:
//                                     if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
//                                         N3_Fit_gam2 = Coefficient_Matrix_Fit_N3_122[index_Coeffs];
//                                     } else {
//                                         N3_Fit_gam2 = Coefficient_Matrix_Fit_N3_122[index_Coeffs] + N3_Fit_gam2*gam2;
//                                     }
//                                     break;
//                                 case N3_123_ENTRY:
//                                     if (i_fit_gam2 == N_Points_gam1 - i_fit_gam1 - 1) {
//                                         N3_Fit_gam2 = Coefficient_Matrix_Fit_N3_123[index_Coeffs];
//                                     } else {
//                                         N3_Fit_gam2 = Coefficient_Matrix_Fit_N3_123[index_Coeffs] + N3_Fit_gam2*gam2;
//                                     }
//                                     break;
//                                 default:
//                                     cout << "Values for VAR_NUM not specified !!!!!!!!!!!!!!!!!!!!!" << endl;
//                                     exit(0);
//                                     break;
//                             }
//                             index_Coeffs--;
//                         }
//                         if (i_fit_gam1 == N_Points_gam1 - 1) {
//                             N3_Fit_gam1 = N3_Fit_gam2;
//                         } else {
//                             N3_Fit_gam1 = N3_Fit_gam2 + N3_Fit_gam1*gam1;
//                         }
//                     }
//                     if (i_fit_f == N_Points_f - 1) {
//                         N3_Fit_f = N3_Fit_gam1;
//                     } else {
//                         N3_Fit_f = N3_Fit_gam1 + N3_Fit_f*norm_f_2;
//                     } 
//                 }
//                 if (i_fit_E == N_Points_E - 1) {
//                     N3_Fit_E = i_fit_E*N3_Fit_f;
//                 } else if (i_fit_E >= 1) {
//                     N3_Fit_E = i_fit_E*N3_Fit_f + N3_Fit_E*ratio_E;
//                 } 
//             }
//             if (m == Order_SH - l) {
//                 N3_Fit_m = N3_Fit_E;
//             } else {
//                 N3_Fit_m = N3_Fit_E + N3_Fit_m*x2;
//             }
//         }
//         if (l == Order_SH) {
//             N3_Fit = N3_Fit_m;
//         } else {
//             N3_Fit = N3_Fit_m + N3_Fit*z2;
//         }
//     }
//     
//     return N3_Fit;
// }
