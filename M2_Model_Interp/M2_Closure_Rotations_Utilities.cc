#ifndef _M2_CLOSURE_FIT_EIGENSTRUCTURE_H_INCLUDED
#include "M2_Closure_Fit_Eigenstructure.h"
#endif // _M2_CLOSURE_FIT_EIGENSTRUCTURE_H_INCLUDED

int    sqr(int x)
                { return x*x; }
double sqr(const double &x)
                { return x*x; }
long double sqr(const long double &x)
                { return x*x; }
                
long double cube(const long double &x)
                { return x*x*x; }

//*************************************************************************
// This routine computes the direction cosines characterizing the rotation 
// matrix that transforms the covariance matrix into a diagonal matrix
//*************************************************************************
void Closure_RT::setup_Frame_Rotation() {
    long double A_11_diag, A_22_diag;
    long double Cos, Sin;
    
    // Now compute the direction cosines characterizing the rotation matrix that transforms
    // the covariance matrix into a diagonal matrix
    Rotation_Matrix_Diag_Covariance(A_11_diag, A_22_diag, W_array, Cos, Sin);
    rotation.x = Cos;
    rotation.y = Sin;
    
    // cout << "Cos = " << Cos << "  " << "Sin = " << Sin << endl;
    
    // Compute the angular moments up to second-order in the new frame
    Rotate(A_11_diag, A_22_diag);
}

void Closure_RT :: Rotate(const long double &A_11_diag, const long double &A_22_diag) {
    // Apply the rotational transformation: U' = T U
    long double E_val, fx_val, fy_val, pxx_val, pxy_val, pyy_val;
    long double pxx_rot, pxy_rot, pyy_rot;
    long double cos, sin;
    
    cos = rotation.x;
    sin = rotation.y;
    
    E_val = W_array[0];
    fx_val = W_array[1];
    fy_val = W_array[2];
    pxx_val = W_array[3];
    pxy_val = W_array[4];
    pyy_val = W_array[5];
    
//     cout << "E_val = " << E_val << "  " << "fx_val = " << fx_val << "  " << "fy_val = " << fy_val << "  " << "pxx_val = " << pxx_val << "  " << "pxy_val = " << pxy_val << "  " << "pyy_val = " << pyy_val << endl;
    // cout << "cos = " << cos << "  " << "sin = " << sin << endl;
    
    // Radiative energy density
    m_values_RT[0] = E_val;
    
    // Rotate flux vector
    m_values_RT[1] = fx_val*cos + fy_val*sin;
    m_values_RT[2] = -fx_val*sin + fy_val*cos;
    
    // Rotate Pressure tensor
    pxx_rot = pxx_val*sqr(cos) + 2.0*cos*sin*pxy_val + pyy_val*sqr(sin);
    pxy_rot = -cos*sin*pxx_val + (sqr(cos)-sqr(sin))*pxy_val + cos*sin*pyy_val;
    pyy_rot = pxx_val*sqr(sin) - 2.0*cos*sin*pxy_val + pyy_val*sqr(cos);
    
    m_values_RT[3] = gamma_1(pxx_rot);
    m_values_RT[4] = gamma_2(pyy_rot);
    
    // cout << "E_val rot = " << m_values_RT[0] << "  " << "fx_val rot = " << m_values_RT[1] << "  " << "fy_val rot = " << m_values_RT[2] << "  " << "gam1 = " << m_values_RT[3] << "  " << "gam2 = " << m_values_RT[4] << endl;
    
    if (fabs(pxy_rot - fx_RT()*fy_RT()) > TOLER_ROT_M2) {
        cout << "pxy() - fx_RT()*fy_RT() = " << pxy_rot - fx_RT()*fy_RT()<< endl;
        cout << "dsyev2 failed !!!!!!!!!!!!!!!!!!!!!!!!" << endl;    
    } 
//     else if (fabs(pxx_rot - sqr(fx_RT()) - A_11_diag) > TOLER_ROT_M2 ) {
//         cout << "pxx() - sqr(fx_RT()) = " << pxx_rot - sqr(fx_RT()) << "  " << "A_11_diag = " << A_11_diag << "  " << "A_22_diag = " << A_22_diag << endl;
//         cout << "dsyev2 failed !!!!!!!!!!!!!!!!!!!!!!!!" << endl;    
//     } else if (fabs(pyy_rot - sqr(fy_RT()) - A_22_diag) > TOLER_ROT_M2 ) {
//         cout << "pyy() - sqr(fy_RT()) = " << pyy_rot - sqr(fy_RT()) << "  " << "A_11_diag = " << A_11_diag << "  " << "A_22_diag = " << A_22_diag << endl;
//         cout << "dsyev2 failed !!!!!!!!!!!!!!!!!!!!!!!!" << endl;    
//     }
    
    R_rot_mat_2D[0][0] = cos;
    R_rot_mat_2D[0][1] = -sin;
    R_rot_mat_2D[1][0] = sin;
    R_rot_mat_2D[1][1] = cos;
}

void Closure_RT :: set_Closure_Derivatives_Original_Basis() {
    long double temp_val, temp_val_0;
    int max_index;
    long double fx_val, fy_val, pxx_val, pxy_val, pyy_val;
    long double d_R_il_R_jm_R_kn;
    
    fx_val = W_array[1];
    fy_val = W_array[2];
    pxx_val = W_array[3];
    pxy_val = W_array[4];
    pyy_val = W_array[5];
        
    for (int index_U = 1; index_U < STATIC_NUM_VAR_RT; index_U++) {
        set_Rotation_Matrix_RT_Derivatives(index_U);
        
        // Only choose non0.0 entries of the matrix mat to improve efficiency
        // and avoid unecessary computations
        if (index_U == 0) {
            max_index = 1;
        } else if (index_U < 3) {
            max_index = 3;
        } else {
            max_index = STATIC_NUM_VAR_RT;
        }
        
        dqxxx_dU[index_U] = 0.0;
        dqxxy_dU[index_U] = 0.0;
        dqxyy_dU[index_U] = 0.0;
        dqyyy_dU[index_U] = 0.0;
        
        for (int l = 0; l < 2; l++) {
            for (int m = 0; m < 2; m++) {
                for (int n = 0; n < 2; n++) {
                    // Calculate {\partial I^{'(3)}_{lmn}}{\partial U_{q}}
                    // q ==> index_U in this case
                    temp_val = 0.0;
                    for (int p = 0; p < STATIC_NUM_VAR_RT; p++) {
                    //for (int p = index_U; p < max_index; p++) {
                        // Calculate T_{pq} + U_{r} {\partial T_{pr}}{\partial U_{q}}
                        temp_val_0 = mat[p][index_U];
                        for (int r = 0; r < STATIC_NUM_VAR_RT; r++) {
                        // for (int r = index_U; r < max_index; r++) {
                            if (r == 0) {
                                temp_val_0 += d_mat[p][r];
                            } else {
                                temp_val_0 += W_array[r]*d_mat[p][r];
                            }
                        }
                        // Now multiply with {\partial I^{'(3)}_{lmn}}{\partial U'_{p}}
                        temp_val += d_q_RT_index(l,m,n,p)*temp_val_0;
                        
                        if (d_q_RT_index(l,m,n,p) != d_q_RT_index(l,m,n,p)) {
                            cout << "d_q_RT_index(l,m,n,p) = " << d_q_RT_index(l,m,n,p) << endl;
                            exit(0);
                        }
                    }
                    
                    if (temp_val != temp_val) {
                        for (int r = 0; r < STATIC_NUM_VAR_RT; r++) {
                            cout << "r = " << r << "  " << "m_values = " << W_array[r] << endl;
                        }
                        
                        for (int p = 0; p < STATIC_NUM_VAR_RT; p++) {
                            for (int r = 0; r < STATIC_NUM_VAR_RT; r++) {
                                cout << "p = " << p << "  " << "r = " << r << "  " << "mat[p][r] = " << mat[p][r] << "  " << "d_mat[p][r] = " << d_mat[p][r] << endl;
                            }
                        }
                        
                        cout << "temp_val = " << temp_val << "  " << "temp_val_0 = " << temp_val_0 << endl;
                        exit(0);
                    }
                    
                    // Compute d_R_0l_R_0m_R_0n
                    d_R_il_R_jm_R_kn = R_rot_mat_2D[0][m]*R_rot_mat_2D[0][n]*d_R_rot_mat_2D[0][l] + R_rot_mat_2D[0][l]*R_rot_mat_2D[0][n]*d_R_rot_mat_2D[0][m] + R_rot_mat_2D[0][l]*R_rot_mat_2D[0][m]*d_R_rot_mat_2D[0][n];
                    
                    dqxxx_dU[index_U] += d_R_il_R_jm_R_kn*q_RT_index(l,m,n) + R_rot_mat_2D[0][l]*R_rot_mat_2D[0][m]*R_rot_mat_2D[0][n]*temp_val;
                    
                    // Compute d_R_0l_R_0m_R_1n
                    d_R_il_R_jm_R_kn = R_rot_mat_2D[0][m]*R_rot_mat_2D[1][n]*d_R_rot_mat_2D[0][l] + R_rot_mat_2D[0][l]*R_rot_mat_2D[1][n]*d_R_rot_mat_2D[0][m] + R_rot_mat_2D[0][l]*R_rot_mat_2D[0][m]*d_R_rot_mat_2D[1][n];
                    
                    dqxxy_dU[index_U] += d_R_il_R_jm_R_kn*q_RT_index(l,m,n) + R_rot_mat_2D[0][l]*R_rot_mat_2D[0][m]*R_rot_mat_2D[1][n]*temp_val;
                    
                    // Compute d_R_0l_R_1m_R_1n
                    d_R_il_R_jm_R_kn = R_rot_mat_2D[1][m]*R_rot_mat_2D[1][n]*d_R_rot_mat_2D[0][l] + R_rot_mat_2D[0][l]*R_rot_mat_2D[1][n]*d_R_rot_mat_2D[1][m] + R_rot_mat_2D[0][l]*R_rot_mat_2D[1][m]*d_R_rot_mat_2D[1][n];
                    
                    dqxyy_dU[index_U] += d_R_il_R_jm_R_kn*q_RT_index(l,m,n) + R_rot_mat_2D[0][l]*R_rot_mat_2D[1][m]*R_rot_mat_2D[1][n]*temp_val;
                    
                    // Compute d_R_1l_R_1m_R_1n
                    d_R_il_R_jm_R_kn = R_rot_mat_2D[1][m]*R_rot_mat_2D[1][n]*d_R_rot_mat_2D[1][l] + R_rot_mat_2D[1][l]*R_rot_mat_2D[1][n]*d_R_rot_mat_2D[1][m] + R_rot_mat_2D[1][l]*R_rot_mat_2D[1][m]*d_R_rot_mat_2D[1][n];
                    
                    dqyyy_dU[index_U] += d_R_il_R_jm_R_kn*q_RT_index(l,m,n) + R_rot_mat_2D[1][l]*R_rot_mat_2D[1][m]*R_rot_mat_2D[1][n]*temp_val;
                }
            }
        }
    }
    
    // Now calculate dN3_dI0
    dqxxx_dU[0] = qxxx_val - fx_val*dqxxx_dU[1] - fy_val*dqxxx_dU[2] - pxx_val*dqxxx_dU[3] - pxy_val*dqxxx_dU[4] - pyy_val*dqxxx_dU[5];
    dqxxy_dU[0] = qxxy_val - fx_val*dqxxy_dU[1] - fy_val*dqxxy_dU[2] - pxx_val*dqxxy_dU[3] - pxy_val*dqxxy_dU[4] - pyy_val*dqxxy_dU[5];
    dqxyy_dU[0] = qxyy_val - fx_val*dqxyy_dU[1] - fy_val*dqxyy_dU[2] - pxx_val*dqxyy_dU[3] - pxy_val*dqxyy_dU[4] - pyy_val*dqxyy_dU[5];
    dqyyy_dU[0] = qyyy_val - fx_val*dqyyy_dU[1] - fy_val*dqyyy_dU[2] - pxx_val*dqyyy_dU[3] - pxy_val*dqyyy_dU[4] - pyy_val*dqyyy_dU[5];
    cout << "Fix this for non-gray !!!!!!!!!!!!!!!" << endl;
    
//     for (int i = 0; i < STATIC_NUM_VAR_RT; i++) {
//         dqxxx_dU[i] = dqxxx_RT[i];
//         dqxxy_dU[i] = dqxxy_RT[i];
//         dqxyy_dU[i] = dqxyy_RT[i];
//         dqyyy_dU[i] = dqyyy_RT[i];
//     }
    
//     for (int i = 0; i < 2; i++) {
//         for (int j = 0; j < 2; j++) {
//             cout << "i = " << i << "  " << "j = " << j << "  " << "R_rot_mat_2D[i][j] = " << R_rot_mat_2D[i][j] << "  " << "d_mat[i][j] = " << d_mat[i][j] << endl;
//         }
//     }
    
    for (int index_U = 0; index_U < STATIC_NUM_VAR_RT; index_U++) {
        if (dqxxx_dU[index_U] != dqxxx_dU[index_U]) {
            cout << endl;
            cout << "Nan in dqxxx_dU[index_U]" << endl;
            cout << "index_U = " << index_U << endl;
            cout << "dqxxx_dU[0] = " << dqxxx_dU[0] << "  " << "dqxxx_dU[1] = " << dqxxx_dU[1] << "  " << "dqxxx_dU[2] = " << dqxxx_dU[2] << "  " << "dqxxx_dU[3] = " << dqxxx_dU[3] << "  " << "dqxxx_dU[4] = " << dqxxx_dU[4] << endl;
            
        cout << "qxxx_RT_val = " << qxxx_RT_val << "   " << "qxxy_RT_val = " << qxxy_RT_val << "   " << "qxyy_RT_val = " << qxyy_RT_val << "   " << "qyyy_RT_val = " << qyyy_RT_val << endl;
        
        cout << "qxxx_val = " << qxxx_val << "   " << "qxxy_val = " << qxxy_val << "   " << "qxyy_val = " << qxyy_val << "   " << "qyyy_val = " << qyyy_val << endl;
            exit(0);
        }
        if (dqxxy_dU[index_U] != dqxxy_dU[index_U]) {
            cout << "Nan in dqxxy_dU[index_U]" << endl;
            exit(0);
        }
        if (dqxyy_dU[index_U] != dqxyy_dU[index_U]) {
            cout << "Nan in dqxyy_dU[index_U]" << endl;
            exit(0);
        }
        if (dqyyy_dU[index_U] != dqyyy_dU[index_U]) {
            cout << "Nan in dqyyy_dU[index_U]" << endl;
            exit(0);
        }
    }
    
//     cout << "fx_val = " << fx_val << " " << "fy_val = " << fy_val << " " << "pxx_val = " << pxx_val << " " << "pxy_val = " << pxy_val << " " << "pyy_val = " << pyy_val << endl;
//     
//     cout << "qxxx_RT_val = " << qxxx_RT_val << "   " << "qxxy_RT_val = " << qxxy_RT_val << "   " << "qxyy_RT_val = " << qxyy_RT_val << "   " << "qyyy_RT_val = " << qyyy_RT_val << endl;
//     
//     cout << "qxxx_val = " << qxxx_val << "   " << "qxxy_val = " << qxxy_val << "   " << "qxyy_val = " << qxyy_val << "   " << "qyyy_val = " << qyyy_val << endl;
//         
//     cout << endl;
//     cout << "dqxxx_RT[0] = " << dqxxx_RT[0] << "  " << "dqxxx_RT[1] = " << dqxxx_RT[1] << "  " << "dqxxx_RT[2] = " << dqxxx_RT[2] << "  " << "dqxxx_RT[3] = " << dqxxx_RT[3] << "  " << "dqxxx_RT[4] = " << dqxxx_RT[4] << "  " << "dqxxx_RT[5] = " << dqxxx_RT[5] << endl;
//     
//     cout << "dqxxy_RT[0] = " << dqxxy_RT[0] << "  " << "dqxxy_RT[1] = " << dqxxy_RT[1] << "  " << "dqxxy_RT[2] = " << dqxxy_RT[2] << "  " << "dqxxy_RT[3] = " << dqxxy_RT[3] << "  " << "dqxxy_RT[4] = " << dqxxy_RT[4] << "  " << "dqxxy_RT[5] = " << dqxxy_RT[5] << endl;
//     
//     cout << "dqxyy_RT[0] = " << dqxyy_RT[0] << "  " << "dqxyy_RT[1] = " << dqxyy_RT[1] << "  " << "dqxyy_RT[2] = " << dqxyy_RT[2] << "  " << "dqxyy_RT[3] = " << dqxyy_RT[3] << "  " << "dqxyy_RT[4] = " << dqxyy_RT[4] << "  " << "dqxyy_RT[5] = " << dqxyy_RT[5] << endl;
//     
//     cout << "dqyyy_RT[0] = " << dqyyy_RT[0] << "  " << "dqyyy_RT[1] = " << dqyyy_RT[1] << "  " << "dqyyy_RT[2] = " << dqyyy_RT[2] << "  " << "dqyyy_RT[3] = " << dqyyy_RT[3] << "  " << "dqyyy_RT[4] = " << dqyyy_RT[4] << "  " << "dqyyy_RT[5] = " << dqyyy_RT[5] << endl;
//     
//     cout << endl;
//     cout << "dqxxx_dU[0] = " << dqxxx_dU[0] << "  " << "dqxxx_dU[1] = " << dqxxx_dU[1] << "  " << "dqxxx_dU[2] = " << dqxxx_dU[2] << "  " << "dqxxx_dU[3] = " << dqxxx_dU[3] << "  " << "dqxxx_dU[4] = " << dqxxx_dU[4] << "  " << "dqxxx_dU[5] = " << dqxxx_dU[5] << endl;
//     
//     
//     cout << "dqxxy_dU[0] = " << dqxxy_dU[0] << "  " << "dqxxy_dU[1] = " << dqxxy_dU[1] << "  " << "dqxxy_dU[2] = " << dqxxy_dU[2] << "  " << "dqxxy_dU[3] = " << dqxxy_dU[3] << "  " << "dqxxy_dU[4] = " << dqxxy_dU[4] << "  " << "dqxxy_dU[5] = " << dqxxy_dU[5] << endl;
//     
//     cout << "dqxyy_dU[0] = " << dqxyy_dU[0] << "  " << "dqxyy_dU[1] = " << dqxyy_dU[1] << "  " << "dqxyy_dU[2] = " << dqxyy_dU[2] << "  " << "dqxyy_dU[3] = " << dqxyy_dU[3] << "  " << "dqxyy_dU[4] = " << dqxyy_dU[4] << "  " << "dqxyy_dU[5] = " << dqxyy_dU[5] << endl;
//     cout << endl;
//     
//     cout << "dqyyy_dU[0] = " << dqyyy_dU[0] << "  " << "dqyyy_dU[1] = " << dqyyy_dU[1] << "  " << "dqyyy_dU[2] = " << dqyyy_dU[2] << "  " << "dqyyy_dU[3] = " << dqyyy_dU[3] << "  " << "dqyyy_dU[4] = " << dqyyy_dU[4] << "  " << "dqyyy_dU[5] = " << dqyyy_dU[5] << endl;
//     cout << endl;
}

// Rotate frame for N3
long double Closure_RT :: Rotate_N3(const Vector2D &norm_dir, const int &CLOSING_FLUX_INDEX) {
    long double q_temp;
    long double sin, cos;
    cos = norm_dir.x;
    sin = norm_dir.y;
    
//     if (flag_compute_closing_fluxes) {
//         cout << "qxxx_RT_val, qxxy_RT_val, qxyy_RT_val, and qyyy_RT_val not been computed yet !!!!!!" << endl;
//         exit(0);
//     }
    
//     qxxx_RT_val = qxxx_RT();
//     qxxy_RT_val = qxxy_RT();
//     qxyy_RT_val = qxyy_RT();
//     qyyy_RT_val = qyyy_RT();
    
    // Rotate third order tensor
    switch (CLOSING_FLUX_INDEX) {
        case N3_111_ENTRY:
            // q_temp = cube(cos)*qxxx_RT_val-3.0*sqr(cos)*sin*qxxy_RT_val+3.0*cos*sqr(sin)*qxyy_RT_val-cube(sin)*qyyy_RT_val;
            q_temp = 0.0;
            for (int l = 0; l < 2; l++) {
                for (int m = 0; m < 2; m++) {
                    for (int n = 0; n < 2; n++) {
                        q_temp += R_rot_mat_2D[0][l]*R_rot_mat_2D[0][m]*R_rot_mat_2D[0][n]*q_RT_index(l,m,n);
                        // cout << "l = " << l << "  " << "m = " << m << "  " << "n = " << n << "  " << "R_rot_mat_2D[0][l] = " << R_rot_mat_2D[0][l] << "  " << "R_rot_mat_2D[0][m] = " << R_rot_mat_2D[0][m] << "  " << "R_rot_mat_2D[0][n] = " << R_rot_mat_2D[0][n] << endl;
                    }
                }
            }
            break;
        case N3_112_ENTRY:
            // q_temp = sqr(cos)*sin*qxxx_RT_val+(cube(cos)-2.0*cos*sqr(sin))*qxxy_RT_val+(cube(sin)-2.0*sqr(cos)*sin)*qxyy_RT_val+cos*sqr(sin)*qyyy_RT_val;
            q_temp = 0.0;
            for (int l = 0; l < 2; l++) {
                for (int m = 0; m < 2; m++) {
                    for (int n = 0; n < 2; n++) {
                        q_temp += R_rot_mat_2D[0][l]*R_rot_mat_2D[0][m]*R_rot_mat_2D[1][n]*q_RT_index(l,m,n);
                    }
                }
            }
            break;
        case N3_122_ENTRY:
            // q_temp = cos*sqr(sin)*qxxx_RT_val+(2.0*sqr(cos)*sin-cube(sin))*qxxy_RT_val+(cube(cos)-2.0*cos*sqr(sin))*qxyy_RT_val-sqr(cos)*sin*qyyy_RT_val;
            q_temp = 0.0;
            for (int l = 0; l < 2; l++) {
                for (int m = 0; m < 2; m++) {
                    for (int n = 0; n < 2; n++) {
                        q_temp += R_rot_mat_2D[0][l]*R_rot_mat_2D[1][m]*R_rot_mat_2D[1][n]*q_RT_index(l,m,n);
                    }
                }
            }
            break;
        case N3_222_ENTRY:
            // q_temp = cube(sin)*qxxx_RT_val+3.0*cos*sqr(sin)*qxxy_RT_val+3.0*sqr(cos)*sin*qxyy_RT_val+cube(cos)*qyyy_RT_val;
            q_temp = 0.0;
            for (int l = 0; l < 2; l++) {
                for (int m = 0; m < 2; m++) {
                    for (int n = 0; n < 2; n++) {
                        q_temp += R_rot_mat_2D[1][l]*R_rot_mat_2D[1][m]*R_rot_mat_2D[1][n]*q_RT_index(l,m,n);
                    }
                }
            }
            break;
        default:
            cout << "Invalid choice for CLOSING_FLUX_INDEX !!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    return q_temp;
}

//  long double Closure_RT :: Rotate_N3(const Vector2D &norm_dir, const int &Entry) {
//      long double q_temp;
// 
// //     if (flag_compute_closing_fluxes) {
// //         cout << "qxxx_RT_val, qxxy_RT_val, qxyy_RT_val, and qyyy_RT_val not been computed yet !!!!!!" << endl;
// //         exit(0);
// //     }
//     
// //     qxxx_RT_val = qxxx_RT();
// //     qxxy_RT_val = qxxy_RT();
// //     qxyy_RT_val = qxyy_RT();
// //     qyyy_RT_val = qyyy_RT();
//      // Rotate third order tensor
//      switch (Entry) {
//          case 1 :
//              q_temp = cube(norm_dir.x)*qxxx_RT_val+3.0*sqr(norm_dir.x)*norm_dir.y*qxxy_RT_val+3.0*norm_dir.x*sqr(norm_dir.y)*qxyy_RT_val+cube(norm_dir.y)*qyyy_RT_val;
//              break;
//          case 2 :
//              q_temp = -sqr(norm_dir.x)*norm_dir.y*qxxx_RT_val+(cube(norm_dir.x)-2.0*norm_dir.x*sqr(norm_dir.y))*qxxy_RT_val+(2.0*sqr(norm_dir.x)*norm_dir.y-cube(norm_dir.y))*qxyy_RT_val+norm_dir.x*sqr(norm_dir.y)*qyyy_RT_val;
//              break;
//          case 3 :
//              q_temp =  norm_dir.x*sqr(norm_dir.y)*qxxx_RT_val+(cube(norm_dir.y)-2.0*sqr(norm_dir.x)*norm_dir.y)*qxxy_RT_val+(cube(norm_dir.x)-2.0*norm_dir.x*sqr(norm_dir.y))*qxyy_RT_val+sqr(norm_dir.x)*norm_dir.y*qyyy_RT_val;
//              break;
//          case 4 :
//              q_temp = -cube(norm_dir.y)*qxxx_RT_val+3.0*norm_dir.x*sqr(norm_dir.y)*qxxy_RT_val-3.0*sqr(norm_dir.x)*norm_dir.y*qxyy_RT_val+cube(norm_dir.x)*qyyy_RT_val;
//              break;
//      }
//      return q_temp;
// }

//*********************************************************************************
// This routine computes the eigenvalues and eigenvectors of the covariance matrix
//*********************************************************************************
void Closure_RT :: Rotation_Matrix_Diag_Covariance(long double &rt1, long double &rt2, const long double *m_values, long double &cs, long double &sn)
// ----------------------------------------------------------------------------
// Calculates the eigensystem of a real symmetric 2x2 matrix
//    [ A  B ]
//    [ B  C ]
// in the form
//    [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
//    [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
// where rt1 >= rt2. Note that this convention is different from the 1.0 used
// in the LAPACK routine DLAEV2, where |rt1| >= |rt2|.
// ----------------------------------------------------------------------------
{
    long double A, B, C;
    
    A = m_values[3] - sqr(m_values[1]);           //pxx() - sqr(fx());
    B = m_values[4] - m_values[1]*m_values[2];    //pxy() - fx()*fy();
    C = m_values[5] - sqr(m_values[2]);           //pyy() - sqr(fy());
    
    double sm, df, delta;
    double temp_val;
    sm = A + C;
    df = A - C;
    delta = sqrt(sqr(df) + 4.0*B*B); // discriminant of characteristic polynomial
    
    rt1 = 0.5*(sm + delta);
    rt2 = 0.5*(sm - delta);
    
    // cout << "fx = " << m_values[1] << "  " << "fy = " << m_values[2] << "  " << "pxy = " << m_values[4] << endl;
    // cout << "B = " << B << endl;
    
    // Calculate eigenvectors
    if (fabs(B) < TOLER_ROT_M2) {
        if (A >= C) {
            temp_val = 1.0;
        } else {
            temp_val = -1.0; 
        }
//         cout << "A = " << A << "  " << "C = " << C << "  " << "A - C = " << A - C << endl;
    } else {
        temp_val = df/delta;
        // cout << "temp_val = " << temp_val << "  " << "1 + temp_val = " << 1 + temp_val << "  " << "1 - temp_val = " << 1 - temp_val << endl;
    }
    
    cs = sqrt(0.5*(1.0 + temp_val));
    sn = sqrt(1.0 - cs*cs);
    // sn = sqrt(0.5*(1.0 - temp_val));
    
    if (B >= 0.0) {
        cs = cs;
        sn = sn;
    } else if (B < 0.0) {
        cs = -cs;
        sn = sn;
    }
    
    if (cs != cs || sn != sn || isinf(cs) || isinf(sn) /*|| (cs*cs + sn*sn != 1)*/) {
        cout << "cs = " << cs << "  " << "sn = " << sn << endl;
        exit(0);
    }
    
//         cout << "cs = " << cs << "  " << "sn = " << sn << endl;
//     cs = 1.0;
//     sn = 0.0;
}

//*********************************************************************************
// This routine computes derivatives of the eigenvectors of the covariance matrix 
// with respect to the conserved solutions
//*********************************************************************************
void Closure_RT :: Rotation_Matrix_Diag_Covariance_Diff(long double &rt1, long double &rt2, const long double *m_values, long double &d_cs, long double &d_sn, const int &Index_U)
// ----------------------------------------------------------------------------
// Calculates the eigensystem of a real symmetric 2x2 matrix
//    [ A  B ]
//    [ B  C ]
// in the form
//    [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
//    [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
// where rt1 >= rt2. Note that this convention is different from the 1.0 used
// in the LAPACK routine DLAEV2, where |rt1| >= |rt2|.
// ----------------------------------------------------------------------------
{
    long double A, B, C;
    
    A = m_values[3] - sqr(m_values[1]);         //pxx() - sqr(fx());
    B = m_values[4] - m_values[1]*m_values[2];  //pxy() - fx()*fy();
    C = m_values[5] - sqr(m_values[2]);         //pyy() - sqr(fy());
    
    long double sm, df, delta;
    long double d_sm, d_df, d_delta;
    long double t, dt, dA, dC, dB;
    long double cs, sn;
    double temp_val, d_temp_val;
    
    dA = 0.0, dB = 0.0, dC = 0.0;
    switch (Index_U) {
        case 0:
            dA = -m_values[3] + 2.0*sqr(m_values[1]);            // -(pxx() - 2.0*sqr(fx()));
            dB = -m_values[4] + 2.0*m_values[1]*m_values[2];         // -(pxy() - 2.0*fx()*fy());
            dC = -m_values[5] + 2.0*sqr(m_values[2]);            // -(pyy() - 2.0*sqr(fy()));
            break;
        case 1:
            dA = -2.0*m_values[1];
            dB = -m_values[2];
            break;
        case 2:
            dB = -m_values[1];
            dC = -2.0*m_values[2];
            break;
        case 3:
            dA = 1.0;
            break;
        case 4:
            dB = 1.0;
            break;
        case 5:
            dC = 1.0;
            break;
        default:
            cout << "Invalid value for Index_U !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
    
    sm = A + C;
    d_sm = dA + dC;
    
    df = A - C;
    d_df = dA - dC;
    
    delta = sqrt(sqr(df) + 4.0*B*B);
    
    if (fabs(B) < TOLER_ROT_M2) {
        if (A >= C) {
            temp_val = 1.0;
        } else {
            temp_val = -1.0; 
        }
        d_temp_val = 0.0;
        d_cs = 0.0;
        d_sn = 0.0;
    } else {
        d_delta = 0.5*(2.0*d_df*df + 8.0*dB*B)/delta;
        temp_val = df/delta;
        d_temp_val = (d_df*delta - df*d_delta)/pow(delta, 2);
        
        // Calculate eigenvectors
        cs = sqrt(0.5*(1.0 + temp_val));
        sn = sqrt(1.0 - cs*cs);
        // sn = sqrt(0.5*(1.0 - temp_val));
        
        d_cs = 0.5*d_temp_val/(2.0*cs);
        d_sn = -0.5*d_temp_val/(2.0*sn);
        
        if (B >= 0.0) {
            d_cs = d_cs;
            d_sn = d_sn;
        } else if (B < 0.0) {
            d_cs = -d_cs;
            d_sn = d_sn;
        } 
    }
    
    if (d_cs != d_cs || d_sn != d_sn || isinf(d_cs) || isinf(d_sn)) {
//         cout << "f2_RT() = " << f2_RT() << "  " << "fx() = " << fx_RT() << "  " << "fy() = " << fy_RT() << endl;
//         cout << "cs = " << cs << "  " << "sn = " << sn << endl;
//         cout << "d_cs = " << d_cs << "  " << "d_sn = " << d_sn << endl;
        exit(0);
    }
//         cout << "d_cs = " << d_cs << "  " << "d_sn = " << d_sn << endl;
        d_cs = 0.0;
        d_sn = 0.0;
}



void Closure_RT::set_Rotation_Matrix_RT_Derivatives(const int &index_U) {
    long double A_11_diag, A_22_diag;
    long double Cos, Sin;
    
    // In the case where the frame is rotated about an angle theta, the angular moments up to second-order
    // in the new frame can be expressed in terms of the angular moments in the original frame as follows
    // 0.0th-order moment: I^('(0)) = I^((0))
    // First-order moments: I^('(1))_{i} = R_{ji} I^((1))_{j}
    // Second-order moments: I^('(2))_{ij} = R_{pi} R_{qj} I^((2))_{pq}
    
    // For the third-order tensor, we also have: I^('(3))_{ijk} = R_{pi} R_{qj} R_{rk} I^((3))_{pqr}
    // or I^((3))_{ijk} = R_{ip} R_{jq} R_{kr} I^('(3))_{pqr}
    // where R is standard rotation matrix
    //     [ cos theta   -sin theta ]
    // R = [                        ]
    //     [ sin theta    cos theta ]
    
    // Now compute the direction cosines characterizing the rotation matrix that transforms
    // the covariance matrix into a diagonal matrix
    Rotation_Matrix_Diag_Covariance(A_11_diag, A_22_diag, W_array, Cos, Sin);
    rotation.x = Cos;
    rotation.y = Sin;
    
    long double rt1, rt2, cs, sn, d_cs, d_sn;
    
    cs = rotation.x;
    sn = rotation.y;
    
    Rotation_Matrix(mat, cs, sn);
    
    R_rot_mat_2D[0][0] = Cos;
    R_rot_mat_2D[0][1] = -Sin;
    R_rot_mat_2D[1][0] = Sin;
    R_rot_mat_2D[1][1] = Cos;
    
    Rotation_Matrix_Diag_Covariance_Diff(rt1, rt2, W_array, d_cs, d_sn, index_U);
    d_Rotation_Matrix(d_mat, cs, sn, d_cs, d_sn);
    d_R_rot_mat_2D[0][0] = d_cs;
    d_R_rot_mat_2D[0][1] = -d_sn;
    d_R_rot_mat_2D[1][0] = d_sn;
    d_R_rot_mat_2D[1][1] = d_cs;
    
    // cout << "cs = " << cs<< "   " << "sn = " << sn << "   " << "d_cs = " << d_cs << "   " << "d_sn = " << d_sn << endl;
}

void Closure_RT :: Rotation_Matrix(long double mat[][STATIC_NUM_VAR_RT], const long double &cos_angle, const long double &sin_angle) {
    for (int i = 0; i < STATIC_NUM_VAR_RT; i++) {
        for (int j = 0; j < STATIC_NUM_VAR_RT; j++) {
            mat[i][j] = 0.0;
        }
    }
    // Transformation matrix from U to U'
    
    mat[0][0] = 1.0;
    mat[1][1] = cos_angle;
    mat[1][2] = sin_angle;
    mat[2][1] = -sin_angle;
    mat[2][2] = cos_angle;
    mat[3][3] = sqr(cos_angle);
    mat[3][4] = 2.0*cos_angle*sin_angle;
    mat[3][5] = sqr(sin_angle);
    mat[4][3] = -cos_angle*sin_angle;
    mat[4][4] = sqr(cos_angle) - sqr(sin_angle);
    mat[4][5] = cos_angle*sin_angle;
    mat[5][3] = sqr(sin_angle);
    mat[5][4] = -2.0*cos_angle*sin_angle;
    mat[5][5] = sqr(cos_angle);
    
//     mat[0][0] = 1.0;
//     mat[1][1] = cos_angle;
//     mat[1][2] = -sin_angle;
//     mat[2][1] = sin_angle;
//     mat[2][2] = cos_angle;
//     mat[3][3] = sqr(cos_angle);
//     mat[3][4] = -2.0*cos_angle*sin_angle;
//     mat[3][5] = sqr(sin_angle);
//     mat[4][3] = cos_angle*sin_angle;
//     mat[4][4] = sqr(cos_angle) - sqr(sin_angle);
//     mat[4][5] = -cos_angle*sin_angle;
//     mat[5][3] = sqr(sin_angle);
//     mat[5][4] = 2.0*cos_angle*sin_angle;
//     mat[5][5] = sqr(cos_angle);
} 

void Closure_RT :: d_Rotation_Matrix(long double mat[][STATIC_NUM_VAR_RT],const long double &cos_angle, const long double &sin_angle, const long double &d_cos_angle, const long double &d_sin_angle) {
    for (int i = 0; i < STATIC_NUM_VAR_RT; i++) {
        for (int j = 0; j < STATIC_NUM_VAR_RT; j++) {
            mat[i][j] = 0.0;
        }
    }
    mat[0][0] = 0.0;
    mat[1][1] = d_cos_angle;
    mat[1][2] = d_sin_angle;
    mat[2][1] = -d_sin_angle;
    mat[2][2] = d_cos_angle;
    mat[3][3] = 2.0*d_cos_angle*cos_angle;
    mat[3][4] = 2.0*(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
    mat[3][5] = 2.0*d_sin_angle*sin_angle;
    mat[4][3] = -(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
    mat[4][4] = 2.0*d_cos_angle*cos_angle - 2.0*d_sin_angle*sin_angle;
    mat[4][5] = (d_cos_angle*sin_angle + cos_angle*d_sin_angle);
    mat[5][3] = 2.0*d_sin_angle*sin_angle;
    mat[5][4] = -2.0*(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
    mat[5][5] = 2.0*d_cos_angle*cos_angle;

//     mat[0][0] = 0.0;
//     mat[1][1] = d_cos_angle;
//     mat[1][2] = -d_sin_angle;
//     mat[2][1] = d_sin_angle;
//     mat[2][2] = d_cos_angle;
//     mat[3][3] = 2.0*d_cos_angle*cos_angle;
//     mat[3][4] = -2.0*(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
//     mat[3][5] = 2.0*d_sin_angle*sin_angle;
//     mat[4][3] = (d_cos_angle*sin_angle + cos_angle*d_sin_angle);
//     mat[4][4] = 2.0*d_cos_angle*cos_angle - 2.0*d_sin_angle*sin_angle;
//     mat[4][5] = -(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
//     mat[5][3] = 2.0*d_sin_angle*sin_angle;
//     mat[5][4] = 2.0*(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
//     mat[5][5] = 2.0*d_cos_angle*cos_angle;
}    

// DenseMatrix Rotation_Matrix(const long double &cos_angle, const long double &sin_angle) {
//     static DenseMatrix mat(STATIC_NUM_VAR_RT,STATIC_NUM_VAR_RT);
//     mat(0,0) = 1.0;
//     mat(1,1) = cos_angle;
//     mat(1,2) = sin_angle;
//     mat(2,1) = -sin_angle;
//     mat(2,2) = cos_angle;
//     mat(3,3) = sqr(cos_angle);
//     mat(3,4) = 2.0*cos_angle*sin_angle;
//     mat(3,5) = sqr(sin_angle);
//     mat(4,3) = -cos_angle*sin_angle;
//     mat(4,4) = sqr(cos_angle) - sqr(sin_angle);
//     mat(4,5) = cos_angle*sin_angle;
//     mat(5,3) = sqr(sin_angle);
//     mat(5,4) = -2.0*cos_angle*sin_angle;
//     mat(5,5) = sqr(cos_angle);
//     return mat;
// } 
// 
// DenseMatrix d_Rotation_Matrix(const long double &cos_angle, const long double &sin_angle, const long double &d_cos_angle, const long double &d_sin_angle) {
//     static DenseMatrix mat(STATIC_NUM_VAR_RT,STATIC_NUM_VAR_RT);
//     mat(0,0) = 0.0;
//     mat(1,1) = d_cos_angle;
//     mat(1,2) = d_sin_angle;
//     mat(2,1) = -d_sin_angle;
//     mat(2,2) = d_cos_angle;
//     mat(3,3) = 2.0*d_cos_angle*cos_angle;
//     mat(3,4) = 2.0*(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
//     mat(3,5) = 2.0*d_sin_angle*sin_angle;
//     mat(4,3) = -(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
//     mat(4,4) = 2.0*d_cos_angle*cos_angle - 2.0*d_sin_angle*sin_angle;
//     mat(4,5) = (d_cos_angle*sin_angle + cos_angle*d_sin_angle);
//     mat(5,3) = 2.0*d_sin_angle*sin_angle;
//     mat(5,4) = -2.0*(d_cos_angle*sin_angle + cos_angle*d_sin_angle);
//     mat(5,5) = 2.0*d_cos_angle*cos_angle;
//     return mat;
// }
