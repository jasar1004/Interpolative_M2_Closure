#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#define _NG_MN_Model_3D_OPTIM_H_INCLUDED

// NG_MN_Model_3D_OPTIM.h
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

#ifndef _MPI_UTILITIES_H_INCLUDED
#include "./MPI_Utilities.h"
#endif // _MPI_UTILITIES_H_INCLUDED

#ifndef _M2_STATE_PARAMETERS_H_INCLUDED
#include "./M2_State_Parameters.h"
#endif // _M2_STATE_PARAMETERS_H_INCLUDED

#include "../../Packages/Permutations/Permutations_With_Order.h"

# define nmax_Lebedev 65

#define PATHVAR "NLOPT_MN_OPTIMIZATION_Path"

#define TOLER_CONSTRAINTS 1e-12
#define TOLER 1e-8

#define CLOSING_FLUX                        100
#define PARTIAL_MOMENTS                     101

#define GRAY                 900
#define NON_GRAY             901

#define BOSE_EINSTEIN        1000
#define HYPERBOLIC_LIMIT     1001
#define LOGARITHMIC_LIMIT    1002

#define PREDEFINED                3000
#define DISCRETE                  3001

#define VAR_E                          2000
#define VAR_N1                         2001
#define VAR_PHI                        2002
#define VAR_THETA                      2003
#define VAR_GAM1                       2004
#define VAR_GAM2                       2005

#define FULL_TRIANGLE                         4000
#define GAM1_EQ_GAM2_EQ_0_5_TRIANGLE          4001
#define GAM1_EQ_GAM2_TRIANGLE                 4002
#define GAM1_EQ_GAM2_EQ_GAM3_TRIANGLE         4003
#define GAM2_EQ_0_TRIANGLE                    4004
#define GAM3_EQ_0_TRIANGLE                    4005

#define DISCRETE_SET_ENERGY                   5000

#define UNIFORM_ALGEBRAIC_MAPPING_N1          6000
#define UNIFORM_ALGEBRAIC_MAPPING_GAM1        6001

#define FACE_A                                7000
#define FACE_B                                7001
#define FACE_C                                7002

#define PLANE_GAM1_GAM2                       8000
#define PLANE_GAM2_GAM3                       8001
#define PLANE_GAM1_GAM3                       8002

#define FORWARD_FINITE_DIFFERENCE             9000
#define BACKWARD_FINITE_DIFFERENCE            9001
#define CENTRAL_FINITE_DIFFERENCE             9002

#define SPECTRUM_ENERGY_FINITE_DIFFERENCE     10000
#define RADIUS_SPHERE_FINITE_DIFFERENCE       10001
#define EDGE_TRIANGLE_FINITE_DIFFERENCE       10002
#define MEDIAN_TRIANGLE_FINITE_DIFFERENCE     10003

#define FLUX_ONE_D                            11

#define EDGE_ONE                              11000
#define EDGE_TWO                              11001

using namespace std;
using namespace nlopt;
// g++ -I/usr/local/include M1_NLOPT_TEST.cc -L/usr/local/lib -lnlopt -lm -o tuutorial

#include "../../Packages/Finite_Difference/differ.hpp"

extern int DISPLAY_ID;
extern int PRIMARY_ID;

extern long double L_vals[11], E_vals[11];

extern int Lebed_Rule_Set[6];
extern int N_quad_points_Circle_Set[4];
extern long double r_l[4]; //declare extern in here and declare whithout extern in .cc file to make it a global variable

extern long double finite_diff_h_N1;
extern long double finite_diff_h_gam;
extern long double tol_grad;

typedef int (*max_ent_obj_grad_hess)(long double *F_obj_grad_hess, const int &NFUN, const int &Index_Angle, void *fdata);

struct record_Ncoeffs {                                                                                                   
    int N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_Points_gam2;
};

struct record_Derivatives_Triangle {                                                        
    long double d_N3_111, d_N3_122, d_N3_123;                                                   
    long double d2_N3_111, d2_N3_122, d2_N3_123;         
    long double d3_N3_111_dgam_dN1_1, d3_N3_122_dgam_dN1_1, d3_N3_123_dgam_dN1_1;
    long double d4_N3_111_dnorm_f_dgam_dN1_1, d4_N3_122_dnorm_f_dgam_dN1_1;
    long double d3_N3_111_dnorm_f_dgam, d3_N3_122_dnorm_f_dgam;
    
    long double d2_N3_111_dgam_dN1_1, d2_N3_122_dgam_dN1_1;
    long double d2_N3_111_dnorm_f_dgam, d2_N3_122_dnorm_f_dgam; 
    long double d3_N3_111_dnorm_f_dgam_dN1_1, d3_N3_122_dnorm_f_dgam_dN1_1;
    
    long double d3_N3_123_dN1_2_dgam, d3_N3_123_dN1_3_dgam;
    long double d4_N3_123_dN1_1_dN1_2_dgam;
    long double d4_N3_123_dN1_1_dN1_3_dgam;
    long double d4_N3_123_dN1_2_dN1_3_dgam;
    long double d5_N3_123_dN1_1_dN1_2_dN1_3_dgam;
};

struct Finite_Diff_Data {
    long double h;
    int order;
    int prec;
};

struct record_N3 {             
    long double ratio_I0;
    long double I0, N1, N1_1, N1_2, N1_3, gam1, gam2, N3_111, N3_122, N3_123;
    long double x_SH, y_SH, z_SH;
    
    // N3_111
    long double dN3_111_dnorm_f, dN3_111_dgam1, dN3_111_dgam2;
    long double dN3_111_dN1_1, dN3_111_dN1_2, dN3_111_dN1_3;
    long double d2_N3_111_dnorm_f_dN1_1;
    long double d2_N3_111_dgam1_dN1_1;
    long double d2_N3_111_dnorm_f_dgam1;
    long double d3_N3_111_dnorm_f_dgam1_dN1_1;
    
    // N3_122
    long double dN3_122_dnorm_f, dN3_122_dgam1, dN3_122_dgam2;
    long double dN3_122_dN1_1, dN3_122_dN1_2, dN3_122_dN1_3;
    long double d2_N3_122_dnorm_f_dN1_1;
    long double d2_N3_122_dnorm_f_dgam1;
    long double d2_N3_122_dnorm_f_dgam2;
    long double d2_N3_122_dgam1_dN1_1;
    long double d2_N3_122_dgam2_dN1_1;
    long double d2_N3_122_dgam1_dgam2;
    long double d3_N3_122_dnorm_f_dgam1_dN1_1;
    long double d3_N3_122_dnorm_f_dgam2_dN1_1;
    long double d3_N3_122_dnorm_f_dgam1_dgam2;
    long double d3_N3_122_dgam1_dgam2_dN1_1;
    long double d4_N3_122_dnorm_f_dgam1_dgam2_dN1_1; // ==> 41
    
    // N3_123
    long double dN3_123_dN1_1, dN3_123_dN1_2, dN3_123_dN1_3, dN3_123_dgam1, dN3_123_dgam2, dN3_123_dgam3;
    long double d2_N3_123_dN1_1_dN1_2, d2_N3_123_dN1_1_dN1_3, d2_N3_123_dN1_2_dN1_3;
    long double d2_N3_123_dgam1_dgam2, d2_N3_123_dgam1_dgam3, d2_N3_123_dgam2_dgam3; // ==> 53
    
    long double d2_N3_123_dgam1_dN1_1, d2_N3_123_dgam2_dN1_1, d2_N3_123_dN1_1_dgam3;
    long double d2_N3_123_dN1_2_dgam1, d2_N3_123_dN1_2_dgam2, d2_N3_123_dN1_2_dgam3;
    long double d2_N3_123_dN1_3_dgam1, d2_N3_123_dN1_3_dgam2, d2_N3_123_dN1_3_dgam3; // ==> 62
    
    long double d3_N3_123_dN1_1_dN1_2_dN1_3;
    
    long double d3_N3_123_dN1_1_dN1_2_dgam1, d3_N3_123_dN1_1_dN1_2_dgam2, d3_N3_123_dN1_1_dN1_2_dgam3;
    long double d3_N3_123_dN1_1_dN1_3_dgam1, d3_N3_123_dN1_1_dN1_3_dgam2, d3_N3_123_dN1_1_dN1_3_dgam3;
    long double d3_N3_123_dN1_2_dN1_3_dgam1, d3_N3_123_dN1_2_dN1_3_dgam2, d3_N3_123_dN1_2_dN1_3_dgam3; // ==> 72
    
    long double d3_N3_123_dgam1_dgam2_dN1_1, d3_N3_123_dgam1_dN1_1_dgam3, d3_N3_123_dgam2_dN1_1_dgam3;
    long double d3_N3_123_dN1_2_dgam1_dgam2, d3_N3_123_dN1_2_dgam1_dgam3, d3_N3_123_dN1_2_dgam2_dgam3;
    long double d3_N3_123_dN1_3_dgam1_dgam2, d3_N3_123_dN1_3_dgam1_dgam3, d3_N3_123_dN1_3_dgam2_dgam3; // ==> 81
    
    long double d4_N3_123_dN1_1_dN1_2_dgam1_dgam2, d4_N3_123_dN1_1_dN1_2_dgam1_dgam3, d4_N3_123_dN1_1_dN1_2_dgam2_dgam3;
    long double d4_N3_123_dN1_1_dN1_3_dgam1_dgam2, d4_N3_123_dN1_1_dN1_3_dgam1_dgam3, d4_N3_123_dN1_1_dN1_3_dgam2_dgam3;
    long double d4_N3_123_dN1_2_dN1_3_dgam1_dgam2, d4_N3_123_dN1_2_dN1_3_dgam1_dgam3, d4_N3_123_dN1_2_dN1_3_dgam2_dgam3; // ==> 90
    
    long double d4_N3_123_dN1_1_dN1_2_dN1_3_dgam1, d4_N3_123_dN1_1_dN1_2_dN1_3_dgam2, d4_N3_123_dN1_1_dN1_2_dN1_3_dgam3;
    long double d5_N3_123_dN1_1_dN1_2_dN1_3_dgam1_dgam2, d5_N3_123_dN1_1_dN1_2_dN1_3_dgam1_dgam3, d5_N3_123_dN1_1_dN1_2_dN1_3_dgam2_dgam3; // ==> 96
    
    long double N3_112, N3_222; // ==> 98
};

struct record_Partial_Moments {             
    long double ratio_I0;
    long double I0, N1_1, N1_2, N1_3, gam1, gam2;
    long double I0_plus, N1_1_plus, N1_2_plus, N2_11_plus, N2_12_plus, N2_22_plus;
    long double N3_111_plus, N3_122_plus, N3_123_plus;
};

struct record_diff_Moments_dLag_Mult {
    long double dIn_dx0, dIn_dx1, dIn_dx2, dIn_dx3, dIn_dx4, dIn_dx5, dIn_dx6, dIn_dx7, dIn_dx8;
};

struct record_Lag_Mult {                                            
    long double I0, N1_1, N1_2, N1_3, gam1, gam2;
    long double x0, x1, x2, x3, x4, x5, x6, x7, x8;
    
    long double &operator[](int &index) {
        static long double temp_val;
        
        switch(index) {
            case 0:
                temp_val = x0;
                break;
            case 1:
                temp_val = x1;
                break;
            case 2:
                temp_val = x2;
                break;
            case 3:
                temp_val = x3;
                break;
            case 4:
                temp_val = x4;
                break;
            case 5:
                temp_val = x5;
                break;
            case 6:
                temp_val = x6;
                break;
            case 7:
                temp_val = x7;
                break;
            case 8:
                temp_val = x8;
                break;
        }
        return temp_val;
    }
    
    void reset() {
        x0 = 0.0; x1 = 0.0; x2 = 0.0; x3 = 0.0; x4 = 0.0; x5 = 0.0; x6 = 0.0; x7 = 0.0; x8 = 0.0;
    }
};

// Function to write a record at a specific position
template<class record_Vals>
ostream& write(ostream& ios, const record_Vals &rec) {
    ios.clear(); // clear any errors
    ios.write(reinterpret_cast<const char*>(&rec), sizeof(rec)); //write the record to file
    return ios; // return the stream (for easy error detection/chaining)
}

// Function to write a record at a specific position
template<class record_Vals>
iostream& read(iostream& ios, record_Vals &rec, int &id) {
 ios.clear(); // clear any errors
 ios.seekg(sizeof(record_Vals) * id); // move to record's position
 ios.read(reinterpret_cast<char *>(&rec), sizeof(rec)); //read record from file
 return ios; // return the stream (for easy error detection/chaining)
}

// Function to write a record at a specific position
template<class record_Vals_1, class record_Vals_2>
iostream& read(iostream& ios, record_Vals_1 &rec, int &id_1, int &id_2) {
 ios.clear(); // clear any errors
 ios.seekg(sizeof(record_Vals_1) * id_1 + sizeof(record_Vals_2) * id_2); // move to record's position
 ios.read(reinterpret_cast<char *>(&rec), sizeof(rec)); //read record from file
 return ios; // return the stream (for easy error detection/chaining)
}

struct Mobius_Scale_Parameters {
    int N_pts_Mob_Scale;
    int Length_Scale_Dist_Type;
    int N_pts_SH, Order_SH, N_pts_f, N_pts_gam1;
    long double *Coefficients_Mobius_Scale_Fit = NULL;
    
    // Constructor
    Mobius_Scale_Parameters(const int &Num_points_f, const int &Order_Sphere_Harmo, const int &N_Coeffs_SH, const int &Num_points_gam1, const long double *Coeffs_Mobius_Scale) {
        Order_SH = Order_Sphere_Harmo;
        N_pts_SH = N_Coeffs_SH;
        N_pts_f = Num_points_f;
        N_pts_gam1 = Num_points_gam1;
        N_pts_Mob_Scale = N_pts_SH*N_pts_f*(N_pts_gam1*(N_pts_gam1 + 1))/2;
        allocate();
        Set_Coefficients_Mob_Scale_Fit(Coeffs_Mobius_Scale);
    }
    
    // Destructor
    ~Mobius_Scale_Parameters() {
        deallocate();
    }
    
    void Set_Coefficients_Mob_Scale_Fit(const long double *Coeffs_Mobius_Scale) {
        for (int i = 0; i < N_pts_Mob_Scale; i++) {
            Coefficients_Mobius_Scale_Fit[i] = Coeffs_Mobius_Scale[i];
        }
    }
    
    void allocate() {
        deallocate();
        Coefficients_Mobius_Scale_Fit = new long double[N_pts_Mob_Scale];
    }
    
    void deallocate() {
        if (Coefficients_Mobius_Scale_Fit != NULL) {
            delete Coefficients_Mobius_Scale_Fit; Coefficients_Mobius_Scale_Fit = NULL;
        }
    }
    
    long double Evaluate_Length_Scale(const long double &norm_f, const long double &x, const long double &y, const long double &z, const long double &gam1, const long double &gam2) {
        long double Length_Scale, L_f, L_m, L_gam1, L_gam2;
        long double norm_f_2 = norm_f*norm_f;
        int index_Coeffs = N_pts_Mob_Scale - 1;
        long double x2, y2, z2;
        x2 = x*x;
        y2 = y*y;
        z2 = z*z;
        
        for (int l = Order_SH; l >= 0; l-=2) {
            for (int m = Order_SH - l; m >= 0; m-=2) {
                for (int i_fit_f = N_pts_f - 1; i_fit_f >= 0; i_fit_f--) {
                    for (int i_fit_gam1 = N_pts_gam1 - 1; i_fit_gam1 >= 0; i_fit_gam1--) {
                        for (int i_fit_gam2 = N_pts_gam1 - i_fit_gam1 - 1; i_fit_gam2 >= 0; i_fit_gam2--) {
                            if (i_fit_gam2 == N_pts_gam1 - i_fit_gam1 - 1) {
                                L_gam2 = Coefficients_Mobius_Scale_Fit[index_Coeffs];
                            } else {
                                L_gam2 = Coefficients_Mobius_Scale_Fit[index_Coeffs] + L_gam2*gam2;
                            }
                            index_Coeffs--;
                        }
                        if (i_fit_gam1 == N_pts_gam1 - 1) {
                            L_gam1 = L_gam2;
                        } else {
                            L_gam1 = L_gam2 + L_gam1*gam1;
                        }
                    }
                    if (i_fit_f == N_pts_f - 1) {
                        L_f = L_gam1;
                    } else {
                        L_f = L_gam1 + L_f*norm_f_2;
                    }
                }
                if (m == Order_SH - l) {
                    L_m = L_f;
                } else {
                    // L_m = L_f + L_m*x2;
                    L_m = L_f + L_m*y2;
                }
            }
            if (l == Order_SH) {
                Length_Scale = L_m;
            } else {
                // Length_Scale = L_m + Length_Scale*z2;
                Length_Scale = L_m + Length_Scale*x2;
            }
        }
        
        // Length_Scale = Coefficients_Mobius_Scale_Fit[0];
        
        if (Length_Scale_Dist_Type == LENGTH_SCALE_DIST_FIT) {
            Length_Scale = exp(Length_Scale);
        } else if (Length_Scale_Dist_Type != LENGTH_SCALE_DIST_UNIF) {
            cout << "Length scale distribution type not specified!!!!!!!!!!!!!!!" << endl;
            exit(0);
        }
        
        return Length_Scale;
    }
};

struct Finite_Diff_Parameters {
    long double x0;
    long double N3_111_x0, N3_122_x0, N3_123_x0;
    long double dN3_111_dN1, d2N3_111_dN1, d3N3_111_dN1, d4N3_111_dN1;
    long double dN3_122_dN1, d2N3_122_dN1, d3N3_122_dN1, d4N3_122_dN1;
    long double dN3_123_dN1, d2N3_123_dN1, d3N3_123_dN1, d4N3_123_dN1;
    
    long double gam1_0, gam2_0;
    long double delta_gam_Triangle;
    
    long double dN3_111_dgam, d2N3_111_dgam, d3N3_111_dgam, d4N3_111_dgam;
    long double dN3_122_dgam, d2N3_122_dgam, d3N3_122_dgam, d4N3_122_dgam;
    long double dN3_123_dgam, d2N3_123_dgam, d3N3_123_dgam, d4N3_123_dgam;
    
    int flag_finite_diff_setup = 0;
    
    void Taylor_Series_N3_ijk(long double &N3_111_Approx, long double &N3_122_Approx, long double &N3_123_Approx, const long double &norm_f) {
        long double delta_x;
        delta_x = norm_f - x0;
        
        N3_111_Approx = N3_111_x0 + delta_x*dN3_111_dN1 + pow(delta_x,2)*d2N3_111_dN1/2.0 + pow(delta_x,3)*d3N3_111_dN1/6.0 + pow(delta_x,4)*d4N3_111_dN1/24.0;
        
        N3_122_Approx = N3_122_x0 + delta_x*dN3_122_dN1 + pow(delta_x,2)*d2N3_122_dN1/2.0 + pow(delta_x,3)*d3N3_122_dN1/6.0 + pow(delta_x,4)*d4N3_122_dN1/24.0;
        
        N3_123_Approx = N3_123_x0 + delta_x*dN3_123_dN1 + pow(delta_x,2)*d2N3_123_dN1/2.0 + pow(delta_x,3)*d3N3_123_dN1/6.0 + pow(delta_x,4)*d4N3_123_dN1/24.0;
    }
    
    void Taylor_Series_N3_ijk_Triangle(long double &N3_111_Approx, long double &N3_122_Approx, long double &N3_123_Approx) {
        long double delta_x;
        delta_x = delta_gam_Triangle;
        
        N3_111_Approx = N3_111_x0 + delta_x*dN3_111_dgam + pow(delta_x,2)*d2N3_111_dgam/2.0 + pow(delta_x,3)*d3N3_111_dgam/6.0 + pow(delta_x,4)*d4N3_111_dgam/24.0;
        
        N3_122_Approx = N3_122_x0 + delta_x*dN3_122_dgam + pow(delta_x,2)*d2N3_122_dgam/2.0 + pow(delta_x,3)*d3N3_122_dgam/6.0 + pow(delta_x,4)*d4N3_122_dgam/24.0;
        
        N3_123_Approx = N3_123_x0 + delta_x*dN3_123_dgam + pow(delta_x,2)*d2N3_123_dgam/2.0 + pow(delta_x,3)*d3N3_123_dgam/6.0 + pow(delta_x,4)*d4N3_123_dgam/24.0;
    }
};


struct my_constraint_data{
    long double a = 0.0, b = 0.0, c = 0.0, f = 0.0, g = 0.0, h = 0.0, i = 0.0, j = 0.0, k = 0.0;
    
    long double operator[](int &index) {
        long double temp_val;
        
        switch(index) {
            case 0:
                temp_val = a;
                break;
            case 1:
                temp_val = b;
                break;
            case 2:
                temp_val = c;
                break;
            case 3:
                temp_val = f;
                break;
            case 4:
                temp_val = g;
                break;
            case 5:
                temp_val = h;
                break;
            case 6:
                temp_val = i;
                break;
            case 7:
                temp_val = j;
                break;
            case 8:
                temp_val = k;
                break;
            default:
                cout << "Invalid value for index !!!!!!!!!!!" << endl;
                exit(0);
                break;
        }
        return temp_val;
    }
    
//     long double &operator[](int &index) {
//         long double temp_val;
//         
// //         if(index >= NVARS || index < 0){
// //             cout << "Invalid values for index" << endl;
// //         }
//         switch(index) {
//             case 0:
//                 temp_val = a;
//                 return a;
//                 break;
//             case 1:
//                 temp_val = b;
//                 return b;
//                 break;
//             case 2:
//                 temp_val = c;
//                 return c;
//                 break;
//             case 3:
//                 temp_val = f;
//                 return f;
//                 break;
//             case 4:
//                 temp_val = g;
//                 return g;
//                 break;
//             case 5:
//                 temp_val = h;
//                 return h;
//                 break;
//             case 6:
//                 temp_val = i;
//                 return i;
//                 break;
//             case 7:
//                 temp_val = j;
//                 return j;
//                 break;
//             case 8:
//                 temp_val = k;
//                 return k;
//                 break;
//         }
// //         return temp_val;
//     }
};

///////////////////////////////////////////////////////
// Functions to be defined in .cc files
///////////////////////////////////////////////////////
long double roundval( const long double &val );

long double roundn( const long double &val, const int &n );

void Modified_Gram_Schmidt_Factorization(const int &n, const long double *A, long double *Q, long double *R);

long double Inverse_Mobius_Transformation(const long double &ratio, const long double &Length_Scale);
 
long double Mobius_Transformation(long double &Moment, const long double &Length_Scale);

void Rotate_Moments(M2_State_Param &M2_State, const long double &phi);

long double Rotation_Matrix_X_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle);

long double Rotation_Matrix_Y_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle);

long double Rotation_Matrix_Z_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle);

int Matrix_to_Vector_Indexing(const int &i, const int &j, const int &NVARS);

void Rotate_Frame_N3_Basis(long double &Omega1, long double &Omega2, long double &Omega3, const long double &theta, const long double &phi);

void Rotate_Frame(M2_State_Param &M2_State, const long double &theta, const long double &phi);

void transpose_matrix(const int &nl, const int &nc, long double *A_transpose, const long double *A);

void Backward_Substitution(const int &nl, const int &nc, long double *A, long double *Q, const long double *R);
 
void Compute_H_ONE_Matrix(long double *Q_data, const long double *x, const long double *Sk, M2_State_Param &M2_State);

void Compute_H_TWO_Matrix(long double *H_TWO, const long double *x, const long double *Sk, M2_State_Param &M2_State);
 
void Compute_A_ONE_Matrix(long double *Q_data, const long double *x, const long double *Sk, M2_State_Param &M2_State);

void Compute_A_TWO_Matrix(long double *Q_data, long double *Inv_Hessian, const int &j, const int &l, M2_State_Param &M2_State);
 
void Compute_A_TWO_Matrix(long double *A_TWO, long double *A_ONE, M2_State_Param &M2_State);

void Compute_A_THREE_Matrix(long double *Q_data, long double *Inv_Hessian, const int &j, const int &l, M2_State_Param &M2_State);

long double Calculate_N3_Second_Derivatives(long double *d2N3dlambda, long double *Inv_Hessian, const int &j, const int &l, M2_State_Param &M2_State);
 
long double Calculate_N3_Third_Derivatives(long double *dN3dlambda, long double *Inv_Hessian, const int &j, const int &l, M2_State_Param &M2_State);

void Triangle_to_Square_Mapping_Nodes(long double &v_x, long double &v_y, const int &i, const int &j, const int &m, const int &Node_Dist_gam1);

void Partial_Lebedev_Quadrature_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata);

void Lebedev_Quadrature_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata);

void Lebedev_Quadrature_Matrix(long double *Vals, const int &NFUN, const long double *Matrix, const int &Quad_Rule);

void Set_Regime(M2_State_Param *M2_State);

long double define_Moment_Set(const long double &E, const long double &fx_test, const long double &fy_test, const long double &fz_test, const long double &pxx_test, const long double &pyy_test, const int &index);

long double define_Moment_Set_Not_RT(const long double &E, const long double &fx_test, const long double &fy_test, const long double &fz_test, const long double &pxx_test, const long double &pxy_test, const long double &pyy_test, const int &index);

void Check_Domain_Type(int &Domain_Type, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, M2_State_Param &M2_State);

void Set_Initial_Guess_And_Tol(M2_State_Param &M2_State);

long double myconstraint(unsigned n, const long double *x, const long double *Sk, long double *grad, void *data);

long double sum_Lag_Mom(const long double *x, const int &NVARS, const long double* MOMENTS);

void add_constraints(nlopt_opt &opt, my_constraint_data *data, M2_State_Param &M2_State);
 
void new_basis_Moments(const int &n, long double *MOMENTS_NEW_BASIS, const long double *MOMENTS, const long double *Sk_cur);

void Regularize_Moment(const long double &r_l, M2_State_Param *M2_State);

void Numerical_Cubature(unsigned fdim, integrand f, unsigned dim, unsigned maxEval, long double *val, void *fdata);

void setup_spherical_harmonics_data (const int &N_Points_Phi, const int &N_Points_Theta, long double *x_Lebed, long double *y_Lebed, long double *z_Lebed, long double *w_quad);

long double generate_monomials_basis(const long double &Omega1, const long double &Omega2, const long double &Omega3, const int &index);
// long double generate_monomials_basis(const long double &mu, const long double &Phi, const int &index);

long double generate_higher_order_monomials_basis(const long double &Omega1, const long double &Omega2, const long double &Omega3, const int &index);

void generate_polynomials_Lebedev(long double* poly, const int &Id_angle, const M2_State_Param &M2_State);

void generate_polynomials(long double* poly, const int &Id_angle, const M2_State_Param &M2_State);

// void generate_polynomials(long double* poly, const long double &mu, const long double &phi, const M2_State_Param &M2_State);

long double sum_Lagrange(const long double *x, const long double *poly_Sk, const int &NVARS);

int N3_M2_3D_DIRECTIONS(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata);
 
bool Check_Realizability_Third_Order_Moment(const long double &N1_1, const long double &N2_11, const long double &N2_22, const long double &N3_111, const long double &N3_122, const long double &N3_123);

void N3_M2_3D_RT(long double *N3, const long double *x, const long double *Sk, M2_State_Param &M2_State);

void Check_Gradient(const long double *x, const long double *Sk, M2_State_Param &M2_State, const long double &grad_SLSQP);

int Partial_Moments_DIRECTIONS(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata);

void Partial_Moments_M2_GRAY_3D_RT(long double *PARTIAL_MOMS, const long double *x, const long double *Sk, M2_State_Param &M2_State);

void Partial_Moments_M2_GRAY_3D(long double *PARTIAL_MOMS, const int &NVARS, const int &NFUN, const long double *x, const long double *Sk, const long double &x_Lebed, const long double &y_Lebed, const long double &z_Lebed) ;
 
long double Uniform_Distribution(const int &i, const int &Np, const long double &val_min, const long double &val_max);

int Check_Realizability(M2_State_Param *M2_State);

void Set_Moments(M2_State_Param *M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, Mobius_Scale_Parameters *Mobius_Scale_Params);

void Set_Iso_Moments(M2_State_Param &M2_State);

int H_ONE_Matrix_Entries(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata);

int Higher_Order_Moments_Derivatives(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata);

int H_TWO_Matrix_Entries(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata);

int Higher_Order_Moments_Second_Derivatives(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata);

int H_THREE_Matrix_Entries(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata);

int Higher_Order_Moments_Third_Derivatives(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata);

void Setup_finite_diff_Triangle_Median(M2_State_Param &M2_State, const long double &h);

void finite_diff_triangle_median(M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere, const long double &h, const int &order, const int &prec, record_Derivatives_Triangle &diff_median);

void finite_diff_triangle_edges(M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere, const long double &h, const int &order, const int &prec, record_Derivatives_Triangle &diff_edge);

void Set_Closure_Free_Streaming(record_N3 *rec_N3_local, const int &id_count, M2_State_Param &M2_State);

void Set_Closure_Triangle_Vertices(record_N3 *rec_N3_local, const int &id_count, M2_State_Param &M2_State);

void Compute_Third_Order_Closing_Fluxes(record_N3 *rec_N3_local, const int &index_count, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, M2_State_Param &M2_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, long double *Sk_final);

void Compute_Third_Order_Closing_Fluxes_Derivatives(record_N3 *rec_N3_local, const int &index_count, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, M2_State_Param &M2_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, long double *Sk_final);

long double dIn_dInm1(const int &i, const int &j, const long double &I0, const long double &N1_1, const long double &N1_2, const long double &N1_3, const long double &gam1, const long double &gam2);

void Calculate_Higher_Order_Moments_Derivatives(record_N3 &rec_N3, const long double *x, const long double *Sk, M2_State_Param &M2_State);

void Setup_Taylor_Series_Coefficients(Finite_Diff_Parameters &Finite_Diff_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed);

void Setup_finite_diff_Parameters_Triangle_Taylor_Series(M2_State_Param &M2_State, Finite_Diff_Parameters &Finite_Diff_N3);

void Setup_Taylor_Series_Coefficients_Triangle(Finite_Diff_Parameters &Finite_Diff_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed);

void Calculate_Higher_Order_Moments_Derivatives_Boundary_Norm_f(record_N3 &rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed);

void Calculate_Higher_Order_Moments_Derivatives_Boundary_gam1_gam2_gam3(record_N3 &rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere = false);

void Calculate_Higher_Order_Moments_Derivatives_Vertex_Median(record_Derivatives_Triangle &rec_dN3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere, const Finite_Diff_Data &rec_Finite_Diff);

void Calculate_Higher_Order_Moments_Derivatives_Vertex_Edges(record_Derivatives_Triangle &rec_dN3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere, const int &Edge_Type, const Finite_Diff_Data &rec_Finite_Diff);

void Calculate_Higher_Order_Moments_Derivatives_Triangle_Vertices(record_N3 &rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere = false);

void Q_func_Gram_Schmidt(unsigned n, const long double *x, const long double *Sk, long double *Q_data, void *my_func_data, void *f_data_MN_State);

int obj_function(long double *F_obj, const int &NFUN, const int &Id_angle, void *fdata);
 
int gradient(long double *grad, const int &NFUN, const int &Id_angle, void *fdata);
 
long double myfunc(unsigned n, const long double *x, const long double *Sk, long double *grad, void *my_func_data);

void NLOPT_Optim_Algo(record_N3 *rec_N3_local, const int &id_count, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed);

int Optimization_Algorithm_Array(record_N3 *rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *num_points, const MPI_proc_params *MPI_proc_parameters, const int &N_pts_Mob_Scale, const long double *Coefficients_Mobius_Scale_Fit);

int Optimization_Algorithm(M2_State_Param &M2_State, MN_Var_Num_Points *num_points, fstream &in_N3_out, const int &N_pts_Mob_Scale, const long double *Coefficients_Mobius_Scale_Fit);

////////////////////////////////////////////////////////////
// Routine for computations over the triangle 
////////////////////////////////////////////////////////////
long double a_median(const M2_State_Param &M2_State);

long double b_median(const M2_State_Param &M2_State);

void Triangle_Median_Edge_Intersection_Points(M2_State_Param &M2_State);

void Median_Line_Unit_Vector(long double &u1, long double &u2, M2_State_Param &M2_State);

void Median_Line_Unit_Vector_Plane_Gam1_Gam3(long double &u1, long double &u2, M2_State_Param &M2_State);

#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED
