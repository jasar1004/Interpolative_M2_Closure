#ifndef _M2_CLOSURE_FIT_EIGENSTRUCTURE_H_INCLUDED
#define _M2_CLOSURE_FIT_EIGENSTRUCTURE_H_INCLUDED

#ifdef FORTRAN_TRAILING_UNDERSCORE
#define F77NAME(x) x##_
#else
#define F77NAME(x) x
#endif

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

using namespace std;

#define STATIC_NUM_VAR_RT   6

#define MOMENT_CLOSURE_M2                         1
#define MOMENT_CLOSURE_PROJECTION                 2

#define TOLER_DENSITY_M2                          1e-8
#define TOLER_FLUX_M2                             1e-8
#define TOLER_M2_FIT                              1e-8

#define N3_111_ENTRY                              10
#define N3_112_ENTRY                              11
#define N3_122_ENTRY                              12
#define N3_222_ENTRY                              13
#define N3_123_ENTRY                              14

#define INDEX_RATIO_I0                            20
#define INDEX_NORM_f_2                            21
#define INDEX_x_SH_2                              22
#define INDEX_y_SH_2                              23
#define INDEX_z_SH_2                              24
#define INDEX_GAM1                                25
#define INDEX_GAM2                                26

#define INDEX_I0                                  30
#define INDEX_N1_1                                31
#define INDEX_N1_2                                32
#define INDEX_N2_11                               33
#define INDEX_N2_12                               34
#define INDEX_N2_22                               35

#define INDEX_IO_STAR                             40
#define INDEX_LENGTH_SCALE_IO_STAR                41

#define R2_ORIGINAL                               50
#define R2_RT                                     51

#define HALF        0.5

#define TOLER_ROT_M2                               1.0e-6

int    sqr(int x);
double sqr(const double &x);
long double sqr(const long double &x);
long double cube(const long double &x);

void Eigenvalues_overwrite( long double *dFdU, long double *REAL, long double *IMAG);

void Rotation_Matrix(long double mat[][STATIC_NUM_VAR_RT], const long double &cos_angle, const long double &sin_angle);
void d_Rotation_Matrix(long double mat[][STATIC_NUM_VAR_RT], const long double &cos_angle, const long double &sin_angle, const long double &d_cos_angle, const long double &d_sin_angle);

// This structure contains the independent variables for our interpolative-based approximations
// of the third-order closing fluxes. The variables are computed once and stored in the structure
// for use whenever needed within the interpolation procedure, instead of the repeatedly computing
// such parameters.
struct Independent_Vars_Closing_Fluxes_N3 {
    long double I0_star;
    long double N1_1, N1_2, N1_3;
    long double norm_f, x_SH, y_SH, z_SH;
    long double gam1, gam2;
};

struct record_d_N3_ijk_RT {
    long double dI0 = 0.0;
    long double dN1_1 = 0.0, dN1_2 = 0.0, dN1_3 = 0.0;
    long double dgam1 = 0.0, dgam2 = 0.0, dB = 0.0;
};

struct record_d_L_I0_star_N3_ijk_RT {
    long double dI0_star = 0.0;
    long double d_norm_f_2 = 0.0, d_x2 = 0.0, d_y2 = 0.0, d_z2 = 0.0;
    long double dgam1 = 0.0, dgam2 = 0.0;
};

struct Jacobian_Independent_Vars_to_Conserved_Vars {
    long double d_norm_f2_dN1_1, d_norm_f2_dN1_2, d_norm_f2_dN1_3;
    long double dx2_dN1_1, dx2_dN1_2, dx2_dN1_3;
    long double dz2_dN1_1, dz2_dN1_2, dz2_dN1_3;
    
    long double dgam1_dN1_1, dgam1_dN1_2, dgam1_dN2_11, dgam1_dN2_12, dgam1_dN2_22;
    long double dgam2_dN1_1, dgam2_dN1_2, dgam2_dN2_11, dgam2_dN2_12, dgam2_dN2_22;
};

struct Vector2D {
    double x,y;  //!< x- and y-components of 2D vector.
                 // Made public so can access them.
		
    //@{ @name Constructors.
    //! Creation constructor.
    Vector2D(void) {
       x = 0.0;
       y = 0.0;
    }
};

class Closure_RT {
protected:

private:

public:
    long double W_array[STATIC_NUM_VAR_RT];
    long double m_values_RT[STATIC_NUM_VAR_RT];
    
    long double dFdU[STATIC_NUM_VAR_RT*STATIC_NUM_VAR_RT];
    
    long double dqxxx_RT[STATIC_NUM_VAR_RT], dqxxy_RT[STATIC_NUM_VAR_RT], dqxyy_RT[STATIC_NUM_VAR_RT], dqyyy_RT[STATIC_NUM_VAR_RT];
    
    long double mat[STATIC_NUM_VAR_RT][STATIC_NUM_VAR_RT];
    long double d_mat[STATIC_NUM_VAR_RT][STATIC_NUM_VAR_RT];
    long double R_rot_mat_2D[2][2];
    long double d_R_rot_mat_2D[2][2];
    
    long double A_11_diag, A_22_diag;
    // Closure in rotated frame RT
    long double qxxx_RT_val, qxxy_RT_val, qxyy_RT_val, qyyy_RT_val;
    Vector2D rotation;
    
    // Closure in Original frame
    long double qxxx_val, qxxy_val, qxyy_val, qyyy_val;
    long double dqxxx_dU[STATIC_NUM_VAR_RT], dqxxy_dU[STATIC_NUM_VAR_RT], dqxyy_dU[STATIC_NUM_VAR_RT], dqyyy_dU[STATIC_NUM_VAR_RT];
    
    record_d_N3_ijk_RT dN3_111_RT, dN3_112_RT, dN3_122_RT, dN3_222_RT;
    
    // Compute W_array
    void Compute_W_array(const long double &I0_star_val, const long double &N1_1, const long double &N1_2, const long double &gam1, const long double &gam2);
    
    void Set_W_array(const long double &I0_star_val, const long double &N1_1, const long double &N1_2, const long double &N2_11, const long double &N2_12, const long double &N2_22);
    
    void Rotate_Moments(const long double &phi);
    
    long double q_index(const int &i, const int &j, const int &k);
    long double d_q_index(const int &i, const int &j, const int &k, const int &index_U);
    long double q_RT_index(const int &i, const int &j, const int &k);
    long double d_q_RT_index(const int &i, const int &j, const int &k, const int &index_U);
    
    //////////////////////////////////////////////////////////////
    // Rotation Utilities
    //////////////////////////////////////////////////////////////
    void setup_Frame_Rotation();
    void Rotate(const long double &A_11_diag, const long double &A_22_diag);
    void set_Rotation_Matrix_RT_Derivatives(const int &index_U);
    void set_Closure_Derivatives_Original_Basis();
    long double Rotate_N3(const Vector2D &norm_dir, const int &CLOSING_FLUX_INDEX);
    void Rotation_Matrix_Diag_Covariance(long double &rt1, long double &rt2, const long double *m_values, long double &cs, long double &sn);
    void Rotation_Matrix_Diag_Covariance_Diff(long double &rt1, long double &rt2, const long double *m_values, long double &d_cs, long double &d_sn, const int &Index_U);
    void Rotation_Matrix(long double mat[][STATIC_NUM_VAR_RT], const long double &cos_angle, const long double &sin_angle);
    void d_Rotation_Matrix(long double mat[][STATIC_NUM_VAR_RT],const long double &cos_angle, const long double &sin_angle, const long double &d_cos_angle, const long double &d_sin_angle);
    
    // 
    void Evaluate_d_q_RT_dU0();
    
    void Evaluate_d_q_Orig_Basis_dU0_Adept();
    
    // Constructor
    Closure_RT(void) {
        for ( int i = 0; i < STATIC_NUM_VAR_RT; i++){
            W_array[i] = 0.0;
            m_values_RT[i] = 0.0;
        }
    }
    
    long double e_RT() {return m_values_RT[0];}
    long double e_RT() const {return m_values_RT[0];}
    long double fx_RT() {return m_values_RT[1];}
    long double fx_RT() const {return m_values_RT[1];}
    long double fy_RT() {return m_values_RT[2];}
    long double fy_RT() const {return m_values_RT[2];}
    long double gam_1() {return m_values_RT[3];}
    long double gam_1() const {return m_values_RT[3];}
    long double gam_2() {return m_values_RT[4];}
    long double gam_2() const {return m_values_RT[4];}
    
    long double fsca();
    long double f2_RT();
    
    long double gamma_1(const long double &N2_xx);
    long double gamma_2(const long double &N2_yy);
    
    long double N2_xx_RT(const long double &N1_1, const long double &N1_2, const long double &gam1, const long double &gam2);
    long double N2_yy_RT(const long double &N1_1, const long double &N1_2, const long double &gam1, const long double &gam2);
    
    //! Index operator.
    long double &operator[](int index) { 
        assert( index >= 1 && index <= STATIC_NUM_VAR_RT );
        return m_values_RT[index-1]; 
    }
    
    const long double &operator[](int index) const {
        assert( index >= 1 && index <= STATIC_NUM_VAR_RT );
        return m_values_RT[index-1]; 
    }
};

inline long double Closure_RT::fsca() {
    long double norm_f;
    
    norm_f = sqrt(f2_RT());
    
    return norm_f;
}

inline long double Closure_RT::f2_RT() {
    long double norm_f_2;
    
    norm_f_2 = fx_RT()*fx_RT() + fy_RT()*fy_RT();
    
    return norm_f_2;
}

inline long double Closure_RT::gamma_1(const long double &N2_xx){
    long double gam1;
    long double norm_f_2 = f2_RT();
    
    if (fabs(1.0 - norm_f_2) < TOLER_M2_FIT) {
        gam1 = 1.0;
    } else {
        gam1 = (N2_xx - pow(fx_RT(), 2))/(1.0 - norm_f_2);
    }
    
    if (gam1 < 0.0) {
        cout << "gam1 = " << gam1 << "  "  << "N2_xx = " << N2_xx << "  "  << "fx_RT = " << fx_RT() << endl;
        exit(0);
    }
    
    return gam1;
}

inline long double Closure_RT::gamma_2(const long double &N2_yy) {
    long double gam2;
    long double norm_f_2 = f2_RT();
    
    if (fabs(1.0 - norm_f_2) < TOLER_M2_FIT) {
        gam2 = 0.0;
    } else {    
        gam2 = (N2_yy - pow(fy_RT(), 2))/(1.0 - norm_f_2);
    }
    
    if (gam2 < 0.0) {
        cout << "gam2 = " << gam2 << "  "  << "N2_yy = " << N2_yy << "  "  << "fy_RT = " << fy_RT() << endl;
        exit(0);
    }
    
    return gam2;
}

inline long double Closure_RT::N2_xx_RT(const long double &N1_1, const long double &N1_2, const long double &gam1, const long double &gam2) {
    long double pxx_RT;
    long double norm_f_2;
    norm_f_2 = pow(N1_1, 2) + pow(N1_2, 2);
    
    pxx_RT = pow(N1_1, 2) + gam1*(1.0 - norm_f_2);
    
    return pxx_RT;
}

inline long double Closure_RT::N2_yy_RT(const long double &N1_1, const long double &N1_2, const long double &gam1, const long double &gam2) {
    long double pyy_RT;
    long double norm_f_2;
    norm_f_2 = pow(N1_1, 2) + pow(N1_2, 2);
    
    pyy_RT = pow(N1_2, 2) + gam2*(1.0 - norm_f_2);
    
    return pyy_RT;
}

inline void Closure_RT :: Compute_W_array(const long double &I0_star_val, const long double &N1_1, const long double &N1_2, const long double &gam1, const long double &gam2) {
    W_array[0] = I0_star_val;
    W_array[1] = N1_1;
    W_array[2] = N1_2;
    W_array[3] = N2_xx_RT(N1_1, N1_2, gam1, gam2);
    W_array[4] = N1_1*N1_2;
    W_array[5] = N2_yy_RT(N1_1, N1_2, gam1, gam2);
}

inline void Closure_RT :: Set_W_array(const long double &I0_star_val, const long double &N1_1, const long double &N1_2, const long double &N2_11, const long double &N2_12, const long double &N2_22) {
    W_array[0] = I0_star_val;
    W_array[1] = N1_1;
    W_array[2] = N1_2;
    W_array[3] = N2_11;
    W_array[4] = N2_12;
    W_array[5] = N2_22;
}

inline void Closure_RT :: Rotate_Moments(const long double &phi) {
    long double fx_val, fy_val, pxx_val, pxy_val, pyy_val;
    long double cos_tmp, sin_tmp;
    cos_tmp = cos(phi);
    sin_tmp = sin(phi);
    
    fx_val = W_array[1];
    fy_val = W_array[2];
    pxx_val = W_array[3];
    pxy_val = W_array[4];
    pyy_val = W_array[5];
    
    // Radiative energy density remains the same
    
    // Rotate flux vector
    W_array[1] = fx_val*cos_tmp + fy_val*sin_tmp;
    W_array[2] = -fx_val*sin_tmp + fy_val*cos_tmp;
    
    // Rotate Pressure tensor
    W_array[3] = pxx_val*sqr(cos_tmp) + 2.0*cos_tmp*sin_tmp*pxy_val + pyy_val*sqr(sin_tmp);
    W_array[4] = -cos_tmp*sin_tmp*pxx_val + (sqr(cos_tmp)-sqr(sin_tmp))*pxy_val + cos_tmp*sin_tmp*pyy_val;
    W_array[5] = pxx_val*sqr(sin_tmp) - 2.0*cos_tmp*sin_tmp*pxy_val + pyy_val*sqr(cos_tmp);
}

inline void Closure_RT :: Evaluate_d_q_RT_dU0() {
    double N2_xx_RT_val, N2_xy_RT_val, N2_yy_RT_val;
    N2_xx_RT_val = W_array[3];
    N2_xy_RT_val = W_array[4];
    N2_yy_RT_val = W_array[5];
    
//     cout << "dN3_111_RT.dI0 = " << dN3_111_RT.dI0 << "   " << "dN3_112_RT.dI0  = " << dN3_112_RT.dI0 << "   " << "dN3_122_RT.dI0  = " << dN3_122_RT.dI0 << "   " << "dN3_222_RT.dI0  = " << dN3_222_RT.dI0 << endl;
    
//     cout << "e = " << e_RT() << "   " << "N1_1 = " << fx_RT() << "   " << "N1_2 = " << fy_RT() << "   " << "N2_xx_RT_val = " << N2_xx_RT_val << "   " << "N2_yy_RT_val = " << N2_yy_RT_val << "   " << "qxxx_RT_val = " << qxxx_RT_val << endl;
    
    //     if (!flag_compute_closing_fluxes) {
    dqxxx_RT[0] = qxxx_RT_val + e_RT()*dN3_111_RT.dI0 - fx_RT()*dqxxx_RT[1] - fy_RT()*dqxxx_RT[2] - N2_xx_RT_val*dqxxx_RT[3] - N2_xy_RT_val*dqxxx_RT[4] - N2_yy_RT_val*dqxxx_RT[5];
    
    dqxxy_RT[0] = qxxy_RT_val + e_RT()*dN3_112_RT.dI0 - fx_RT()*dqxxy_RT[1] - fy_RT()*dqxxy_RT[2] - N2_xx_RT_val*dqxxy_RT[3] - N2_xy_RT_val*dqxxy_RT[4] - N2_yy_RT_val*dqxxy_RT[5];
    
    dqxyy_RT[0] = qxyy_RT_val + e_RT()*dN3_122_RT.dI0 - fx_RT()*dqxyy_RT[1] - fy_RT()*dqxyy_RT[2] - N2_xx_RT_val*dqxyy_RT[3] - N2_xy_RT_val*dqxyy_RT[4] - N2_yy_RT_val*dqxyy_RT[5];
    
    dqyyy_RT[0] = qyyy_RT_val + e_RT()*dN3_222_RT.dI0 - fx_RT()*dqyyy_RT[1] - fy_RT()*dqyyy_RT[2] - N2_xx_RT_val*dqyyy_RT[3] - N2_xy_RT_val*dqyyy_RT[4] - N2_yy_RT_val*dqyyy_RT[5];    
//     } else {
//         cout << "qxxx_RT_val has not been computed yet !!!!!!" << endl;
//         exit(0);
//     }
    
//     cout << "dqxxx_RT[0] = " << dqxxx_RT[0] << "  " << "dqxxx_RT[1] = " << dqxxx_RT[1] << "  " << "dqxxx_RT[2] = " << dqxxx_RT[2] << "  " << "dqxxx_RT[3] = " << dqxxx_RT[3] << "  " << "dqxxx_RT[4] = " << dqxxx_RT[4] << "  " << "dqxxx_RT[5] = " << dqxxx_RT[5] << endl;
    
}

inline void Closure_RT:: Evaluate_d_q_Orig_Basis_dU0_Adept() {
    dqxxx_dU[0] = qxxx_val + W_array[0]*dqxxx_dU[0] - W_array[1]*dqxxx_dU[1] - W_array[2]*dqxxx_dU[2] - W_array[3]*dqxxx_dU[3] - W_array[4]*dqxxx_dU[4] - W_array[5]*dqxxx_dU[5];
    
    dqxxy_dU[0] = qxxy_val + W_array[0]*dqxxy_dU[0] - W_array[1]*dqxxy_dU[1] - W_array[2]*dqxxy_dU[2] - W_array[3]*dqxxy_dU[3] - W_array[4]*dqxxy_dU[4] - W_array[5]*dqxxy_dU[5];
    
    dqxyy_dU[0] = qxyy_val + W_array[0]*dqxyy_dU[0] - W_array[1]*dqxyy_dU[1] - W_array[2]*dqxyy_dU[2] - W_array[3]*dqxyy_dU[3] - W_array[4]*dqxyy_dU[4] - W_array[5]*dqxyy_dU[5];
    
    dqyyy_dU[0] = qyyy_val + W_array[0]*dqyyy_dU[0] - W_array[1]*dqyyy_dU[1] - W_array[2]*dqyyy_dU[2] - W_array[3]*dqyyy_dU[3] - W_array[4]*dqyyy_dU[4] - W_array[5]*dqyyy_dU[5];    
    
    // cout << "dqxxx_dU[0] = " << dqxxx_dU[0] << "  " << "dqxxx_dU[1] = " << dqxxx_dU[1] << "  " << "dqxxx_dU[2] = " << dqxxx_dU[2] << "  " << "dqxxx_dU[3] = " << dqxxx_dU[3] << "  " << "dqxxx_dU[4] = " << dqxxx_dU[4] << "  " << "dqxxx_dU[5] = " << dqxxx_dU[5] << endl;
}

//***********************************************************************
// This routine ...
//***********************************************************************
// inline void Closure_RT :: Setup_Independent_Vars_Closing_Fluxes_N3(Independent_Vars_Closing_Fluxes_N3 &Vars_N3_RT, int v) {
//     switch(Absorption_Model) {
//         case MEDIUM2D_ABSORB_GRAY : 
//         case MEDIUM2D_ABSORB_FSCK :
//             Vars_N3_RT.I0_star = e_RT();
//             break;
//         case MEDIUM2D_ABSORB_SNBCK:
//             Vars_N3_RT.I0_star = I0_star(v);
//             break;
//     }
//     
//     if (Vars_N3_RT.I0_star < ZERO) {
//         cout << "e_RT() = " << e_RT() << "  " << "I0_star = " << Vars_N3_RT.I0_star << endl;
//         exit(0);
//     }
//     
//     Vars_N3_RT.N1_1 = fx_RT();
//     Vars_N3_RT.N1_2 = fy_RT();
//     Vars_N3_RT.N1_3 = 0.0;
//     
//     Cartesian_to_spherical_Coordinates(Vars_N3_RT.norm_f, Vars_N3_RT.x_SH, Vars_N3_RT.y_SH, Vars_N3_RT.z_SH, Vars_N3_RT.N1_1, Vars_N3_RT.N1_2, Vars_N3_RT.N1_3);
//     
//     Vars_N3_RT.gam1 = gam_1();
//     Vars_N3_RT.gam2 = gam_2();
// }

// inline void Closure_RT :: Cartesian_to_spherical_Coordinates(long double &norm_f, long double &x_SH, long double &y_SH, long double &z_SH, const long double &N1_1, const long double &N1_2, const long double &N1_3) {
//     norm_f = sqrt(pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2));
//     
//     x_SH = N1_1/norm_f;
//     z_SH = N1_3/norm_f;
//     
//     if (norm_f < 1.0e-8) {
//         x_SH = 0.0;
//         y_SH = 0.0;
//         z_SH = 1.0;
//     }
// }

#endif // _RADMOM2D_STATE_SECOND_ORDER_RT_INCLUDED
