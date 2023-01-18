#ifndef _M2_STATE_PARAMETERS_H_INCLUDED
#define _M2_STATE_PARAMETERS_H_INCLUDED

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
#include "/home/g/groth/jasar10/CFFC/nlopt/include/nlopt.hpp"
// #include "../../../../../CFFC/nlopt/include/nlopt.hpp"
#include "../../Packages/cubature_1_0_3/cubature.h"

#include "./circle_rule.hpp"

#ifndef _SPHERE_LEBEDEV_RULE_HPP_INCLUDED
#include "../../Packages/Lebedev_Quadrature/sphere_lebedev_rule.hpp"
#endif // _SPHERE_LEBEDEV_RULE_HPP_INCLUDED

#ifndef _ONE_DIMENSION_QUADRATURE_HPP_INCLUDED
#include "../../Packages/One_Dimensional_Quadratures/quadrule.hpp"
#endif // _ONE_DIMENSION_QUADRATURE_HPP_INCLUDED

#ifndef _CHEBYSHEV_HPP_INCLUDED
#include "../../Packages/Chebyshev/chebyshev.hpp"
#endif // _CHEBYSHEV_HPP_INCLUDED

#define PATHVAR "NLOPT_MN_OPTIMIZATION_Path"

#define ONE_DIMENSIONAL      1
#define TWO_DIMENSIONAL      2
#define THREE_DIMENSIONAL    3

#define BOUNDARY_FREE_STREAMING        3000
#define GAM2_EQUAL_GAM3                3001
#define BOUNDARY_GAM1                  3002
#define BOUNDARY_GAM2                  3003
#define BOUNDARY_GAM3                  3004
#define GAM1_GAM2_GAM3                 3005
#define BOUNDARY_GAM1_GAM2_EQ_0        3006
#define BOUNDARY_GAM1_GAM3_EQ_0        3007
#define BOUNDARY_GAM2_GAM3_EQ_0        3008

#define MOMENTS_TEST_NORM_F                  4000
#define MOMENTS_TEST_SPHERICAL_HARMONIC      4001
#define TEST_SINGLE_SET                      4002

#define LENGTH_SCALE_DIST_FIT               10000000
#define LENGTH_SCALE_DIST_UNIF              10000001

#define PI 3.14159265358979323846

using namespace std;
using namespace nlopt;

struct MN_Var_Num_Points{
    int E = 0, f = 0, theta = 0, phi = 0, gam1 = 0, gam2 = 0;
    int Order_SH, N_pts_SH;
    int Maximum_Entropy_Solution_Regime;
    int Length_Scale_Dist_Type = 0;
    int N_Points_Triangle_gam1_gam2 = 0;
};

struct Record_Moments_Tests {
    long double ratio_I0 = -1.0, I0 = 0.0, N1 = 0.0, N1_1 = 0.0, N1_2 = 0.0, N1_3 = 0.0;
    long double gamma_1 = 0.0, gamma_2 = 0.0;
    long double x_SH = 0.0, y_SH = 0.0, z_SH = 0.0;
    int flag_Moments_Test = false;
    int flag_Moments_Test_Type = MOMENTS_TEST_SPHERICAL_HARMONIC;
    int Problem_Type;
    
    void Override_Moments_Tests(long double &I0_val, long double &ratio_I0_val, long double &norm_f_val, long double &x_SH_val, long double &y_SH_val, long double &z_SH_val, long double &gam1_val, long double &gam2_val, int &Problem_Type_val) {
        Problem_Type_val = Problem_Type;
        switch (flag_Moments_Test_Type) {
            case MOMENTS_TEST_NORM_F:
                ratio_I0_val = ratio_I0;
                I0_val = I0;
                x_SH_val = x_SH;
                y_SH_val = y_SH;
                z_SH_val = z_SH;
                gam1_val = gamma_1;
                gam2_val = gamma_2;
                break;
            case MOMENTS_TEST_SPHERICAL_HARMONIC:
                ratio_I0_val = ratio_I0;
                I0_val = I0;
                norm_f_val = N1;
                gam1_val = gamma_1;
                gam2_val = gamma_2;
                break;
            case TEST_SINGLE_SET:
                ratio_I0_val = ratio_I0;
                I0_val = I0;
                norm_f_val = N1;
                x_SH_val = x_SH;
                y_SH_val = y_SH;
                z_SH_val = z_SH;
                gam1_val = gamma_1;
                gam2_val = gamma_2;
                break;
            default:
                cout << "flag_Moments_Test_Type not specified !!!!!!!!!!!!!!" << endl;
                exit(0);
                break;
        }
    }
    
    void Copy(Record_Moments_Tests &Rec_Moms_Test) {
        ratio_I0 = Rec_Moms_Test.ratio_I0;
        I0 = Rec_Moms_Test.I0;
        N1 = Rec_Moms_Test.N1;
        N1_1 = Rec_Moms_Test.N1_1;
        N1_2 = Rec_Moms_Test.N1_2;
        N1_3 = Rec_Moms_Test.N1_3;
        x_SH = Rec_Moms_Test.x_SH;
        y_SH = Rec_Moms_Test.y_SH;
        z_SH = Rec_Moms_Test.z_SH;
        gamma_1 = Rec_Moms_Test.gamma_1;
        gamma_2 = Rec_Moms_Test.gamma_2;
        
        flag_Moments_Test = Rec_Moms_Test.flag_Moments_Test;
        flag_Moments_Test_Type = Rec_Moms_Test.flag_Moments_Test_Type;
        
        Problem_Type = Rec_Moms_Test.Problem_Type;
    }
};

struct M2_State_Param {
    int num_proc_E = -10;
    int num_proc_f = -10;
    int num_proc_Phi = -10;
    int num_proc_Theta = -10;
    int num_proc_Triangle_gam1_gam2 = -10;
    
    mutable Record_Moments_Tests Rec_Moms_Test;
    static int Max_Ent_Solution_Type;
    long double I0 = 0.0, N1 = 0.0, N1_1 = 0.0, N1_2 = 0.0, N1_3 = 0.0;
    long double gamma_1 = 0.0, gamma_2 = 0.0, gamma_3 = 0.0;
    long double x_SH = 0.0, y_SH = 0.0, z_SH = 0.0;
    long double N2_11, N2_12, N2_22;
    long double ratio_I0;
    int face_triangle = 0;
    int Triangle_Plane = 0;
    int id_finite_diff_ratio_I0 = 0;
    int id_finite_diff_norm_f = 0;
    int id_finite_diff_gams = 0;
    long double h_ratio_I0 = 0.0;
    long double h_norm_f = 0.0;
    long double h_gam1 = 0.0;
    long double h_gam2 = 0.0;
    long double h_gam = 0.0;
    long double a_median = 0.0;
    long double b_median = 0.0;
    long double ratio_I0_knot = 0.0;
    long double norm_f_knot = 0.0;
    long double gam1_knot = 0.0;
    long double gam2_knot = 0.0;
    long double gam3_knot = 0.0;
    long double gam1_bound_min = 0.0;
    long double gam2_bound_min = 0.0;
    long double gam1_bound_max = 0.0;
    long double gam2_bound_max = 0.0;
    long double cos_Median_Line = 1.0;
    long double sin_Median_Line = 0.0;
    long double x_finite_diff_ratio_I0 = 0.0;
    long double x_finite_diff_norm_f = 0.0;
    long double x_finite_diff_gam1 = 0.0;
    long double x_finite_diff_gam2 = 0.0;
    int Finite_Diff_Spectrum_Boundary_Type;
    int finite_diff_type = 0;
    int finite_diff_domain_Spectrum = 0;
    int finite_diff_domain_Sphere = 0;
    int finite_diff_domain_Triangle = 0;
    bool finite_diff_Sphere = false;
    bool finite_diff_Triangle = false;
    int NVARS = 9;
    int Lebedev_rule = 10;
    int Order_Quad_mu = 1;
    int *Index = NULL;
    long double *MOM_VEC = NULL;
    long double *MOM_VEC_iso = NULL;
    long double *Sk = NULL;
    long double *Sk_cur = NULL;
    long double *ak = NULL;
    long double *x = NULL;
    long double *tol_x = NULL;
    int index_p = 0;
    int index_a = 0;
    int index_Moment = 0;
    int index_Lag_i = 0;
    int index_Lag_j = 0;
    int index_Lag_k = 0;
    int Index_Higher_Order_Moments = 0;
    int flag_Realizability = 0;
    static int Regime; 
    int Problem_Type = 0, Dimension = 0, Domain_Type = 0;
    long double Length_Scale_Mobius;
    long double theta_N1 = 0.0, phi_N1 = 0.0;
    int Triangle_Domain;
    int display = false;
    int proc_display = 0;
    int id_proc = 0, num_proc = 1;
    bool flag_Taylor_Series_Expansion_Sphere = false;
    bool flag_Taylor_Series_Expansion_Triangle = false;
    
    bool flag_finite_diff_Sphere = false;
    bool flag_finite_diff_Triangle = false;
    
    int Node_Dist_E, Node_Dist_f, Node_Dist_Phi_Theta, Node_Dist_gam1;
    
    // Weight and abscissas for Lebedev quadrature
    long double *phi_quad = NULL;
    long double *w_phi_quad = NULL;
    long double *mu_quad = NULL;
    long double *w_mu_quad = NULL;
    
    long double *Omega1 = NULL, *Omega2 = NULL, *Omega3 = NULL;
    long double *weight_total = NULL;
    
    // Parameters for quadrature refinement
    static const int max_refinement_level = 20;
    long double bounds_quad_mu_refin[2*max_refinement_level+3];
    long double bounds_quad_phi_refin[4*max_refinement_level+3];
    int refinement_level = 0;
    long double mu_peak[2], phi_peak[2];
    int N_peaks_mu, N_peaks_phi;
    int N_intervals_mu_quad = 2, N_intervals_phi_quad = 4;
    
    void Find_Peak_Locations();
    void Cartesian_to_Spherical_Coordinates(long double &theta, long double &phi, const long double &Omega1, const long double &Omega2, const long double &Omega3);
    void set_bounds_quad_refin_mu(const int &id_refinement_level);
    void set_bounds_quad_refin_phi(const int &id_refinement_level);
    void set_bounds_quad_refin(const int &id_refinement_level);
    
    long double linear_interpolation(const long double &a_new, const long double &b_new, const long double &val_orig, const long double &a_orig, const long double &b_orig);
    
    long double diff_linear_interpolation(const long double &a_new, const long double &b_new, const long double &a_orig, const long double &b_orig);
    
    int N_dirs();
    
    int Id_Angle(const int &Id_angle_mu, const int &Id_angle_Phi);
    
    long double interp1d(const long double &x, const long double &a_orig, const long double &b_orig, const long double &a, const long double &b);
    
    void get_angles_N1();
    
    void compute_Omegas(const long double &mu_start, const long double &mu_end, const long double &phi_start, const long double &phi_end);
   
    void Set_Lebed_Quad(void);
    
    void Allocate();
    
    void Allocate_Quad();

    void Deallocate();
    
    void Deallocate_Quad();
    
    void Set_Num_Vars(void);
    
    static int Get_Num_Var(const int &Dimension, const int &Domain_Type);
    
    void set_x(const long double *x_iter) {
	    for (int i = 0; i < NVARS; i++) {
		    x[i] = x_iter[i];
	    }
    }

    void set_Sk(const long double *Sk_iter) {
	    for (int i = 0; i < NVARS*NVARS; i++) {
		Sk[i] = Sk_iter[i];
	    }
    }

    void set_Sk_cur(const long double *Sk_cur_iter) {
	 for (int i = 0; i < NVARS*NVARS; i++) {
	     Sk_cur[i] = Sk_cur_iter[i];
	 }
    } 

    void set_ak(const long double *ak_cur) {
	 for (int i = 0; i < NVARS*NVARS; i++) {
	     ak[i] = ak_cur[i];
	 }
    }
    
    void Set_Indices(void);
    
    void Reset_Num_Vars(const int &Dom_Type) {
        Deallocate();
        Domain_Type = Dom_Type;
        Set_Num_Vars();
        Allocate();
        Set_Indices();
    }
    
     M2_State_Param(const int &NUM_VARS) { // Constructor
        NVARS = NUM_VARS;
        Dimension = THREE_DIMENSIONAL;
        Domain_Type = GAM1_GAM2_GAM3;
//         Lebedev_rule = Lebed_rule;
//         Order_Quad_mu = order_table ( Lebedev_rule );
        Allocate();
        Set_Indices();
//         Set_Lebed_Quad();
    }
    
    M2_State_Param(const int &Dim, const int &Dom_Type,
                   const int &Prob_Type) { // Constructor
        Dimension = Dim;
        Domain_Type = Dom_Type;
        Problem_Type = Prob_Type;
//         Lebedev_rule = Lebed_rule;
//         Order_Quad_mu = order_table ( Lebedev_rule );
        Set_Num_Vars();
        Allocate();
        Set_Indices();
//         Set_Lebed_Quad();
    }
    
    void copy(const M2_State_Param &M2_State);
    
    M2_State_Param(const M2_State_Param &M2_State);
    
    M2_State_Param &operator=(const M2_State_Param &M2_State) {
        copy(M2_State);
        Allocate();
        Set_Indices();
        set_bounds_quad_refin(refinement_level);
	    return (*this);
    }
    
   ~M2_State_Param() { Deallocate();} // Destructor
   
   void Setup_Quadrature(const int &id_refinement_level);
};

inline void M2_State_Param :: Setup_Quadrature(const int &id_refinement_level) {
    // Order_Quad_mu = pow(2,id_refinement_level)*5;
    
    Allocate_Quad();
    compute_Omegas(-1.0, 1.0, 0.0, 2.0*PI);
    
    set_bounds_quad_refin(id_refinement_level);
}

inline M2_State_Param :: M2_State_Param(const M2_State_Param &M2_State) {
    copy(M2_State);
    Allocate();
    Set_Indices();
    set_bounds_quad_refin(refinement_level);
}

inline void M2_State_Param :: copy(const M2_State_Param &M2_State) {
    M2_State.Rec_Moms_Test.Copy(Rec_Moms_Test);
    
    num_proc_E = M2_State.num_proc_E;
    num_proc_f = M2_State.num_proc_f;
    num_proc_Phi = M2_State.num_proc_Phi;
    num_proc_Theta = M2_State.num_proc_Theta;
    num_proc_Triangle_gam1_gam2 = M2_State.num_proc_Triangle_gam1_gam2;
    
    NVARS = M2_State.NVARS;
    Dimension = M2_State.Dimension;
    Domain_Type = M2_State.Domain_Type;
    Problem_Type = M2_State.Problem_Type;
    Regime = M2_State.Regime;
    
    ratio_I0 = M2_State.ratio_I0;
    I0 = M2_State.I0;
    N1 = M2_State.N1;
    N1_1 = M2_State.N1_1;
    N1_2 = M2_State.N1_2;
    N1_3 = M2_State.N1_3;
    
    x_SH = M2_State.x_SH;
    y_SH = M2_State.y_SH;
    z_SH = M2_State.z_SH;
    
    gamma_1 = M2_State.gamma_1;
    gamma_2 = M2_State.gamma_2;
    gamma_3 = M2_State.gamma_3;
    N2_11 = M2_State.N2_11;
    N2_12 = M2_State.N2_12;
    N2_22 = M2_State.N2_22;
    
    face_triangle = M2_State.face_triangle;
    Triangle_Plane = M2_State.Triangle_Plane;
    
    id_finite_diff_ratio_I0 = M2_State.id_finite_diff_ratio_I0;
    id_finite_diff_norm_f = M2_State.id_finite_diff_norm_f;
    id_finite_diff_gams = M2_State.id_finite_diff_gams;
    h_ratio_I0 = M2_State.h_ratio_I0;
    h_norm_f = M2_State.h_norm_f;
    h_gam1 = M2_State.h_gam1;
    h_gam2 = M2_State.h_gam2;
    h_gam = M2_State.h_gam;
    ratio_I0_knot = M2_State.ratio_I0_knot;
    norm_f_knot = M2_State.norm_f_knot;
    a_median = M2_State.a_median;
    b_median = M2_State.b_median;
    gam1_knot = M2_State.gam1_knot;
    gam2_knot = M2_State.gam2_knot;
    gam3_knot = M2_State.gam3_knot;
    gam1_bound_min = M2_State.gam1_bound_min;
    gam2_bound_min = M2_State.gam2_bound_min;
    gam1_bound_max = M2_State.gam1_bound_max;
    gam2_bound_max = M2_State.gam2_bound_max;
    cos_Median_Line = M2_State.cos_Median_Line;
    sin_Median_Line = M2_State.sin_Median_Line;
    
    x_finite_diff_ratio_I0 = M2_State.x_finite_diff_ratio_I0;
    x_finite_diff_norm_f = M2_State.x_finite_diff_norm_f;
    x_finite_diff_gam1 = M2_State.x_finite_diff_gam1;
    x_finite_diff_gam2 = M2_State.x_finite_diff_gam2;
    Finite_Diff_Spectrum_Boundary_Type = M2_State.Finite_Diff_Spectrum_Boundary_Type;
    finite_diff_type = M2_State.finite_diff_type;
    finite_diff_domain_Spectrum = M2_State.finite_diff_domain_Spectrum;
    finite_diff_domain_Sphere = M2_State.finite_diff_domain_Sphere;
    finite_diff_domain_Triangle = M2_State.finite_diff_domain_Triangle;
    finite_diff_Sphere = M2_State.finite_diff_Sphere;
    finite_diff_Triangle = M2_State.finite_diff_Triangle;
    Lebedev_rule = M2_State.Lebedev_rule;
    Order_Quad_mu = M2_State.Order_Quad_mu;
    
    index_p = M2_State.index_p;
    index_a = M2_State.index_a;
    index_Moment = M2_State.index_Moment;
    index_Lag_i = M2_State.index_Lag_i;
    index_Lag_j = M2_State.index_Lag_j;
    index_Lag_k = M2_State.index_Lag_k;
    Index_Higher_Order_Moments = M2_State.Index_Higher_Order_Moments;
    flag_Realizability = M2_State.flag_Realizability;
    Length_Scale_Mobius = M2_State.Length_Scale_Mobius;
    theta_N1 = M2_State.theta_N1;
    phi_N1 = M2_State.phi_N1;
    Triangle_Domain = M2_State.Triangle_Domain;
    display = M2_State.display;
    proc_display = M2_State.proc_display;
    id_proc = M2_State.id_proc;
    num_proc = M2_State.num_proc;
    
    flag_Taylor_Series_Expansion_Sphere = M2_State.flag_Taylor_Series_Expansion_Sphere;
    flag_Taylor_Series_Expansion_Triangle = M2_State.flag_Taylor_Series_Expansion_Triangle;
    flag_finite_diff_Sphere = M2_State.flag_finite_diff_Sphere;
    flag_finite_diff_Triangle = M2_State.flag_finite_diff_Triangle;
    
    Node_Dist_E = M2_State.Node_Dist_E;
    Node_Dist_f = M2_State.Node_Dist_f;
    Node_Dist_Phi_Theta = M2_State.Node_Dist_Phi_Theta;
    Node_Dist_gam1 = M2_State.Node_Dist_gam1;
    
    phi_peak[0] = M2_State.phi_peak[0];
    phi_peak[1] = M2_State.phi_peak[1];
    mu_peak[0] = M2_State.mu_peak[0];
    mu_peak[1] = M2_State.mu_peak[1];
    refinement_level = M2_State.refinement_level;
}

inline long double M2_State_Param :: linear_interpolation(const long double &a_new, const long double &b_new, const long double &val_orig, const long double &a_orig, const long double &b_orig) {
    long double interp;
    
    interp = a_new + (b_new - a_new)*(val_orig - a_orig)/(b_orig - a_orig);
    
    return interp;
}

inline long double M2_State_Param :: diff_linear_interpolation(const long double &a_new, const long double &b_new, const long double &a_orig, const long double &b_orig) {
    long double dinterp;
    
    dinterp = (b_new - a_new)/(b_orig - a_orig);
        
    return dinterp;
}

int inline M2_State_Param :: N_dirs() {
    int num_dirs;
    switch (Domain_Type) {
        case BOUNDARY_GAM1:
        case BOUNDARY_GAM2:
        case BOUNDARY_GAM3:
            num_dirs = 2*Order_Quad_mu;
            break;
        case GAM1_GAM2_GAM3:
        case BOUNDARY_GAM1_GAM2_EQ_0:
        case BOUNDARY_GAM1_GAM3_EQ_0:
        case BOUNDARY_GAM2_GAM3_EQ_0:
        case BOUNDARY_FREE_STREAMING:
            num_dirs = 2*Order_Quad_mu*Order_Quad_mu;
            break;
        default:
            cout << "Domain Type not defined" << endl;
            exit(0);
            break;
    }
    return num_dirs;
}

inline int M2_State_Param :: Id_Angle(const int &Id_angle_mu, const int &Id_angle_Phi) {
    int index;
    index = Id_angle_Phi*Order_Quad_mu + Id_angle_mu;
    return index;
}

inline long double M2_State_Param :: interp1d(const long double &x, const long double &a_orig, const long double &b_orig, const long double &a, const long double &b) {
    long double x_new;
    x_new = a + (b - a)*(x - a_orig)/(b_orig - a_orig);
    return x_new;
}

inline void M2_State_Param :: get_angles_N1() {
    long double norm_f;
    norm_f = pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2);
    norm_f = sqrt(norm_f);
    theta_N1 = N1_1/norm_f;
    theta_N1 = acos(theta_N1);
    phi_N1 = N1_2/sqrt(pow(N1_2, 2) + pow(N1_3, 2));
    phi_N1 = acos(phi_N1);
}

//**************************************************************************************
// This routine computes the locations of the potential peaks for the M2 distribution
// and sorts them in increasing order
//**************************************************************************************
inline void M2_State_Param :: Find_Peak_Locations() {
    long double Omega1_peak, Omega2_peak, Omega3_peak;
    long double theta_peak_val, phi_peak_val;
    long double gamma_3 = 1.0 - gamma_1 - gamma_2;
    
    if (gamma_1 < 1.0/3.0) {
        N_peaks_mu = 1;
        N_peaks_phi = 2;
        if (gamma_2 < 1.0/3.0) {
            // In this case gam3 ==> 1
            // The distribution is zero everywhere except along the direction 
            // s = [\omega_1, \omega_2, \omega_3] = [N1_1, N1_2, \pm \sqrt(1 - N1_1^2 - N1_2^2)]
            Omega1_peak = N1_1;
            Omega2_peak = N1_2;
            Omega3_peak = sqrt(1.0 - pow(N1_1, 2) - pow(N1_2, 2));
        } else if (gamma_3 < 1.0/3.0) {
            // In this case gam2 ==> 1
            // The distribution is zero everywhere except along the direction 
            // s = [\omega_1, \omega_2, \omega_3] = [N1_1, \pm \sqrt(1 - N1_1^2 - N1_3^2), N1_3]
            Omega1_peak = N1_1;
            Omega2_peak = sqrt(1.0 - pow(N1_1, 2) - pow(N1_3, 2));
            Omega3_peak = N1_3;
        } else {
            // In this case gam1 ==> 0
            // The distribution is zero everywhere except on the surface \omega_1 = N1_1
            Omega1_peak = N1_1;
            if (gamma_2 > gamma_3) {
                // Then gam2 ==> 1
                // s = [N1_1, \pm \sqrt(1 - N1_1^2 - N1_3^2), N1_3]
                Omega2_peak = sqrt(1.0 - pow(N1_1, 2) - pow(N1_3, 2));
                Omega3_peak = N1_3;
            } else {
                // Then gam3 ==> 1
                // s = [N1_1, N1_2, \pm \sqrt(1 - N1_1^2 - N1_2^2)]
                Omega2_peak = N1_2;
                Omega3_peak = sqrt(1.0 - pow(N1_1, 2) - pow(N1_2, 2));
            }
        }
    } else if (gamma_2 < 1.0/3.0) {
        if (gamma_3 < 1.0/3.0) {
            // In this case gam1 ==> 1
            // The distribution is zero everywhere except along the direction 
            // s = [\omega_1, \omega_2, \omega_3] = [\pm \sqrt(1 - N1_2^2 - N1_3^2), N1_2, N1_3]
            N_peaks_mu = 2;
            N_peaks_phi = 1;
        
            Omega1_peak = sqrt(1.0 - pow(N1_2, 2) - pow(N1_3, 2));
            Omega2_peak = N1_2;
            Omega3_peak = N1_3;
        } else {
            // In this case gam2 ==> 0
            // The distribution is zero everywhere except on the surface \omega_2 = N1_2
            if (gamma_1 > gamma_3) {
                // In this case gam1 ==> 1
                // The distribution is zero everywhere except along the direction 
                // s = [\pm \sqrt(1 - N1_2^2 - N1_3^2), N1_2, N1_3]
                N_peaks_mu = 2;
                N_peaks_phi = 1;
                
                Omega1_peak = sqrt(1.0 - pow(N1_2, 2) - pow(N1_3, 2));
                Omega2_peak = N1_2;
                Omega3_peak = N1_3;
            } else {
                // Then gam3 ==> 1
                // s = [N1_1, N1_2, \pm \sqrt(1 - N1_1^2 - N1_2^2)]
                N_peaks_mu = 1;
                N_peaks_phi = 2;
            
                Omega1_peak = N1_1;
                Omega2_peak = N1_2;
                Omega3_peak = sqrt(1.0 - pow(N1_1, 2) - pow(N1_2, 2));
            }
        }
    } else { // Then gamma_3 < 1.0/3.0
        // In this case gam3 ==> 0
        // The distribution is zero everywhere except on the surface \omega_3 = N1_3
        if (gamma_1 > gamma_2) {
            // Then gam1 ==> 1
            // s = [\pm \sqrt(1 - N1_2^2 - N1_3^2), N1_2, N1_3]
            N_peaks_mu = 2;
            N_peaks_phi = 1;
            
            Omega1_peak = sqrt(1.0 - pow(N1_2, 2) - pow(N1_3, 2));
            Omega2_peak = N1_2;
            Omega3_peak = N1_3;
        } else {
            // Then gam2 ==>1
            // s = [N1_1, \pm \sqrt(1 - N1_1^2 - N1_3^2), N1_3]
            N_peaks_mu = 1;
            N_peaks_phi = 2;
                
            Omega1_peak = N1_1;
            Omega2_peak = sqrt(1.0 - pow(N1_1, 2) - pow(N1_3, 2));
            Omega3_peak = N1_3;
        }
    }
    
    Cartesian_to_Spherical_Coordinates(theta_peak_val, phi_peak_val, Omega1_peak, Omega2_peak, Omega3_peak);
    
    phi_peak[0] = phi_peak_val;
    mu_peak[0] = cos(theta_peak_val);
            
    if (gamma_1 < 1.0/3.0) {
        mu_peak[1] = mu_peak[0];
        if (gamma_2 < 1.0/3.0) {
            // In this case gam3 ==> 1
            // The distribution is zero everywhere except along the direction 
            // s = [\omega_1, \omega_2, \omega_3] = [N1_1, N1_2, \pm \sqrt(1 - N1_1^2 - N1_2^2)]
            phi_peak[1] = 2.0*PI - phi_peak[0];
        } else if (gamma_3 < 1.0/3.0) {
            // In this case gam2 ==> 1
            // The distribution is zero everywhere except along the direction 
            // s = [\omega_1, \omega_2, \omega_3] = [N1_1, \pm \sqrt(1 - N1_1^2 - N1_3^2), N1_3]
            phi_peak[1] = PI - phi_peak[0];
        } else {
            // In this case gam1 ==> 0
            // The distribution is zero everywhere except on the surface \omega_1 = N1_1
            if (gamma_2 > gamma_3) {
                // Then gam2 ==> 1
                // s = [N1_1, \pm \sqrt(1 - N1_1^2 - N1_3^2), N1_3]
                phi_peak[1] = PI - phi_peak[0];
            } else {
                // Then gam3 ==> 1
                // s = [N1_1, N1_2, \pm \sqrt(1 - N1_1^2 - N1_2^2)]
                phi_peak[1] = 2.0*PI - phi_peak[0];
            }
        }
    } else if (gamma_2 < 1.0/3.0) {
        if (gamma_3 < 1.0/3.0) {
            // In this case gam1 ==> 1
            // The distribution is zero everywhere except along the direction 
            // s = [\omega_1, \omega_2, \omega_3] = [\pm \sqrt(1 - N1_2^2 - N1_3^2), N1_2, N1_3]
            mu_peak[1] = -mu_peak[0];
            phi_peak[1] = phi_peak[0];
        } else {
            // In this case gam2 ==> 0
            // The distribution is zero everywhere except on the surface \omega_2 = N1_2
            if (gamma_1 > gamma_3) {
                // In this case gam1 ==> 1
                // The distribution is zero everywhere except along the direction 
                // s = [\pm \sqrt(1 - N1_2^2 - N1_3^2), N1_2, N1_3]
                mu_peak[1] = -mu_peak[0];
                phi_peak[1] = phi_peak[0];
            } else {
                // Then gam3 ==> 1
                // s = [N1_1, N1_2, \pm \sqrt(1 - N1_1^2 - N1_2^2)]
                mu_peak[1] = mu_peak[0];
                phi_peak[1] = 2.0*PI - phi_peak[0];
            }
        }
    } else { // Then gamma_3 < 1.0/3.0
        // In this case gam3 ==> 0
        // The distribution is zero everywhere except on the surface \omega_3 = N1_3
        if (gamma_1 > gamma_2) {
            // Then gam1 ==> 1
            // s = [\pm \sqrt(1 - N1_2^2 - N1_3^2), N1_2, N1_3]
            mu_peak[1] = -mu_peak[0];
            phi_peak[1] = phi_peak[0];
        } else {
            // Then gam2 ==>1
            // s = [N1_1, \pm \sqrt(1 - N1_1^2 - N1_3^2), N1_3]
            mu_peak[1] = mu_peak[0];
            phi_peak[1] = PI - phi_peak[0];
        }
    }
    
    for (int i = 0; i < 2; i++) {
        // Make sure phi is in [0.0, 2.0*PI]
        if (phi_peak[i] < 0.0) {
            phi_peak[i] = 2.0*PI + phi_peak[i];
        } else if (phi_peak[i] > 2.0*PI) {
            phi_peak[i] = phi_peak[i] - 2.0*PI;
        }
    }
    
    // Sort peaks in increasing order
    long double phi_peak_max, phi_peak_min;
    long double mu_peak_max, mu_peak_min;
    
    if (N_peaks_mu == 2) {
        mu_peak_min = min(mu_peak[0], mu_peak[1]);
        mu_peak_max = max(mu_peak[0], mu_peak[1]);
        
        mu_peak[0] = mu_peak_min;
        mu_peak[1] = mu_peak_max;
    }
    
    if (N_peaks_phi == 2) {
        phi_peak_min = min(phi_peak[0], phi_peak[1]);
        phi_peak_max = max(phi_peak[0], phi_peak[1]);
        
        phi_peak[0] = phi_peak_min;
        phi_peak[1] = phi_peak_max;
    }
    
//     cout << "N1_1 = " << N1_1 << "  " << "N1_2 = " << N1_2 << "  " << "N1_3 = " << N1_3 << "  " << "gamma_1 = " << gamma_1 << "  " << "gamma_2 = " << gamma_2 << endl;
//     cout << "mu_peak[0] = " << mu_peak[0] << "  " << "mu_peak[1] = " << mu_peak[1] << endl;
//     cout << "phi_peak[0] = " << phi_peak[0] << "  " << "phi_peak[1] = " << phi_peak[1] << endl;
}

inline void M2_State_Param :: Cartesian_to_Spherical_Coordinates(long double &theta, long double &phi, const long double &Omega1, const long double &Omega2, const long double &Omega3) {
    long double norm_Omegas, norm_Omegas_partial;
    long double cos_phi, sin_phi;
    
    norm_Omegas = pow(Omega1, 2) + pow(Omega2, 2) + pow(Omega3, 2);
    norm_Omegas = sqrt(norm_Omegas);
    
    norm_Omegas_partial = pow(Omega2, 2) + pow(Omega3, 2);
    norm_Omegas_partial = sqrt(norm_Omegas_partial);
    
    theta = Omega1/norm_Omegas;
    theta = acos(theta);
    
    if (norm_Omegas_partial < 1.0e-8) {
        cos_phi = 1.0;
        sin_phi = 0.0;
    } else {
        cos_phi = Omega2/norm_Omegas_partial;
        sin_phi = Omega3/norm_Omegas_partial;
    }
    
    phi = acos(cos_phi);
    
    if (sin_phi < 0.0) {
       phi = 2.0*PI - phi; 
    }
}

inline void M2_State_Param :: set_bounds_quad_refin_mu(const int &id_refinement_level) {
    long double bound_max, bound_min, bound_intermediate;
    long double bound_peak_1, bound_peak_2;
    long double temp_val;
    long double mu_min, mu_max;
    int i_start = 0;
    bool flag_refine_block;
    bool flag_peak;
    
    // We first start by splitting the interval [-1, 1]
    bound_min = -1.0;
    bound_intermediate = 0.0;
    bound_max = 1.0;
//     bound_min = -1.0;
//     bound_peak_1 = mu_peak[0];
//     bound_peak_2 = mu_peak[1];
//     bound_max = 1.0;
    
    long double bounds_quad_mu_refin_coarse[N_intervals_mu_quad+1];
    
    if (refinement_level == 0) {
        bounds_quad_mu_refin[0] = bound_min;
        bounds_quad_mu_refin[1] = bound_intermediate;
        bounds_quad_mu_refin[2] = bound_max;
        N_intervals_mu_quad = 2;
        
        for (int id_peak = 0; id_peak < N_peaks_mu; id_peak++) {
            for (int i = 0; i < N_intervals_mu_quad + 1; i++) {
                if (fabs(mu_peak[id_peak] - bounds_quad_mu_refin[i]) < 1.0e-12) {
                    mu_peak[id_peak] = bounds_quad_mu_refin[i];
                }
            }
        }
        
//         if (N_peaks_mu == 1) {
//             bounds_quad_mu_refin[0] = bound_min;
//             bounds_quad_mu_refin[1] = bound_peak_1;
//             bounds_quad_mu_refin[2] = bound_max;
//             N_intervals_mu_quad = 2;
//         } else if (N_peaks_mu == 2) {
//             bounds_quad_mu_refin[0] = bound_min;
//             bounds_quad_mu_refin[1] = bound_peak_1;
//             bounds_quad_mu_refin[2] = bound_peak_2;
//             bounds_quad_mu_refin[3] = bound_max;
//             N_intervals_mu_quad = 3;
//         } else {
//             cout << "Invalid value for N_peaks_mu !!!!" << endl;
//             exit(0);
//         }
    } else {
        for (int i = 0; i < N_intervals_mu_quad + 1; i++) {
            bounds_quad_mu_refin_coarse[i] = bounds_quad_mu_refin[i];
        }
        
        for (int i = 0; i < N_intervals_mu_quad; i++) {
            mu_min = bounds_quad_mu_refin_coarse[i];
            mu_max = bounds_quad_mu_refin_coarse[i+1];
            
            for (int id_peak = 0; id_peak < N_peaks_mu; id_peak++) {
                if (mu_peak[id_peak] >= mu_min && mu_peak[id_peak] <= mu_max) {
                    // In this case one of the peaks is within the interval of interest
                    flag_refine_block = true;
                    break;
                } else {
                    // In this case none of the peaks is within the interval of interest
                    flag_refine_block = false;
                }
            }
            if (flag_refine_block) {
                // Then the interval [mu_min, mu_max] must be refined
                temp_val = (mu_min + mu_max)/2.0;
                bounds_quad_mu_refin[i_start] = mu_min; 
                bounds_quad_mu_refin[i_start+1] = temp_val;
                bounds_quad_mu_refin[i_start+2] = mu_max;
                i_start = i_start + 2;
            } else {
                // Then the interval [mu_min, mu_max] must not be refined
                bounds_quad_mu_refin[i_start] = mu_min;
                bounds_quad_mu_refin[i_start+1] = mu_max;
                i_start = i_start + 1;
            }
        }
        N_intervals_mu_quad = i_start;
    }
    
    long double bounds_quad_mu_refin_fine[N_intervals_mu_quad+1];
    long double mu_peak_val;
    if (refinement_level == id_refinement_level) {
        // Update block to properly isolate peak
        for (int i = 0; i < N_intervals_mu_quad + 1; i++) {
            bounds_quad_mu_refin_fine[i] = bounds_quad_mu_refin[i];
        }
        
        i_start = 0;
        for (int i = 0; i < N_intervals_mu_quad; i++) {
            mu_min = bounds_quad_mu_refin_fine[i];
            mu_max = bounds_quad_mu_refin_fine[i+1];
            for (int id_peak = 0; id_peak < N_peaks_mu; id_peak++) {
                if (mu_peak[id_peak] > mu_min && mu_peak[id_peak] < mu_max) {
                    flag_peak = true;
                    mu_peak_val = mu_peak[id_peak];
                    break;
                } else {
                    flag_peak = false;
                }
            }
            if (flag_peak) {
                bounds_quad_mu_refin[i_start] = mu_min; 
                bounds_quad_mu_refin[i_start+1] = mu_peak_val;
                bounds_quad_mu_refin[i_start+2] = mu_max;
                i_start = i_start + 2;
            } else {
                bounds_quad_mu_refin[i_start] = mu_min;
                bounds_quad_mu_refin[i_start+1] = mu_max;
                i_start = i_start + 1;
            }
        }
        N_intervals_mu_quad = i_start;
    }
}

inline void M2_State_Param :: set_bounds_quad_refin_phi(const int &id_refinement_level) {
    long double bound_max, bound_min, bound_intermediate_1, bound_intermediate_2, bound_intermediate_3;
    long double bound_peak_1, bound_peak_2;
    long double temp_val;
    long double phi_min, phi_max;
    int i_start = 0;
    bool flag_refine_block;
    bool flag_peak;
    
    // We first start by splitting the interval [0, 2 PI]
    bound_min = 0.0;
    bound_intermediate_1 = PI/2.0;
    bound_intermediate_2 = PI;
    bound_intermediate_3 = 3.0*PI/2.0;
    bound_max = 2.0*PI;
    
//     bound_min = 0.0;
//     bound_peak_1 = phi_peak[0];
//     bound_peak_2 = phi_peak[1];
//     bound_max = 2.0*PI;
    
    long double bounds_quad_phi_refin_coarse[N_intervals_phi_quad+1];
    
    if (refinement_level == 0) {
        bounds_quad_phi_refin[0] = bound_min;
        bounds_quad_phi_refin[1] = bound_intermediate_1;
        bounds_quad_phi_refin[2] = bound_intermediate_2;
        bounds_quad_phi_refin[3] = bound_intermediate_3;
        bounds_quad_phi_refin[4] = bound_max;
        N_intervals_phi_quad = 4;
        
        for (int id_peak = 0; id_peak < N_peaks_phi; id_peak++) {
            for (int i = 0; i < N_intervals_phi_quad + 1; i++) {
                if (fabs(phi_peak[id_peak] - bounds_quad_phi_refin[i]) < 1.0e-12) {
                    phi_peak[id_peak] = bounds_quad_phi_refin[i];
                }
            }
        }
        
//         if (N_peaks_phi == 1) {
//             bounds_quad_phi_refin[0] = bound_min;
//             bounds_quad_phi_refin[1] = bound_peak_1;
//             bounds_quad_phi_refin[2] = bound_max;
//             N_intervals_phi_quad = 2;
//         } else if (N_peaks_phi == 2) {
//             bounds_quad_phi_refin[0] = bound_min;
//             bounds_quad_phi_refin[1] = bound_peak_1;
//             bounds_quad_phi_refin[2] = bound_peak_2;
//             bounds_quad_phi_refin[3] = bound_max;
//             N_intervals_phi_quad = 3;
//         } else {
//             cout << "Invalid value for N_peaks_phi !!!!" << endl;
//             exit(0);
//         }
    } else {
        for (int i = 0; i < N_intervals_phi_quad + 1; i++) {
            bounds_quad_phi_refin_coarse[i] = bounds_quad_phi_refin[i];
        }
        
        for (int i = 0; i < N_intervals_phi_quad; i++) {
            phi_min = bounds_quad_phi_refin_coarse[i];
            phi_max = bounds_quad_phi_refin_coarse[i+1];
            for (int id_peak = 0; id_peak < N_peaks_phi; id_peak++) {
                if (phi_peak[id_peak] == 0.0 && (phi_min == 0.0 || phi_max == 2.0*PI)) {
                    flag_refine_block = true;
                    break;
                } else if (phi_peak[id_peak] >= phi_min && phi_peak[id_peak] <= phi_max) {
                    // In this case one of the peaks is within the interval of interest
                    flag_refine_block = true;
                    break;
                } else {
                    // In this case none of the peaks is within the interval of interest
                    flag_refine_block = false;
                }
            }
            if (flag_refine_block) {
                // Then the interval [mu_min, mu_max] must be refined
                temp_val = (phi_min + phi_max)/2.0;
                bounds_quad_phi_refin[i_start] = phi_min; 
                bounds_quad_phi_refin[i_start+1] = temp_val;
                bounds_quad_phi_refin[i_start+2] = phi_max;
                i_start = i_start + 2;
            } else {
                // Then the interval [mu_min, mu_max] must not be refined
                bounds_quad_phi_refin[i_start] = phi_min;
                bounds_quad_phi_refin[i_start+1] = phi_max;
                i_start = i_start + 1;
            }
        }
        N_intervals_phi_quad = i_start;
    }
    
    
    long double bounds_quad_phi_refin_fine[N_intervals_phi_quad+1];
    long double phi_peak_val;
    if (refinement_level == id_refinement_level) {
        // Update block to properly isolate peak
        for (int i = 0; i < N_intervals_phi_quad + 1; i++) {
            bounds_quad_phi_refin_fine[i] = bounds_quad_phi_refin[i];
        }
        
        i_start = 0;
        for (int i = 0; i < N_intervals_phi_quad; i++) {
            phi_min = bounds_quad_phi_refin_fine[i];
            phi_max = bounds_quad_phi_refin_fine[i+1];
            for (int id_peak = 0; id_peak < N_peaks_phi; id_peak++) {
                if (phi_peak[id_peak] > phi_min && phi_peak[id_peak] < phi_max) {
                    flag_peak = true;
                    phi_peak_val = phi_peak[id_peak];
                    break;
                } else {
                    flag_peak = false;
                }
            }
            if (flag_peak) {
                bounds_quad_phi_refin[i_start] = phi_min; 
                bounds_quad_phi_refin[i_start+1] = phi_peak_val;
                bounds_quad_phi_refin[i_start+2] = phi_max;
                i_start = i_start + 2;
            } else {
                bounds_quad_phi_refin[i_start] = phi_min;
                bounds_quad_phi_refin[i_start+1] = phi_max;
                i_start = i_start + 1;
            }
        }
        N_intervals_phi_quad = i_start;
    }
}

inline void M2_State_Param :: set_bounds_quad_refin(const int &id_refinement_level) {
    Find_Peak_Locations();
    
    switch (Domain_Type) {
        case GAM1_GAM2_GAM3:
        case BOUNDARY_GAM1_GAM2_EQ_0:
        case BOUNDARY_GAM1_GAM3_EQ_0:
        case BOUNDARY_GAM2_GAM3_EQ_0:
        case BOUNDARY_FREE_STREAMING:
            for (int i = 0; i <= id_refinement_level; i++) {
                refinement_level = i;
                set_bounds_quad_refin_mu(id_refinement_level);
                set_bounds_quad_refin_phi(id_refinement_level);
            }
            break;
        case BOUNDARY_GAM1:
        case BOUNDARY_GAM2:
        case BOUNDARY_GAM3:
            for (int i = 0; i <= id_refinement_level; i++) {
                refinement_level = i;
                set_bounds_quad_refin_phi(id_refinement_level);
            }
            break;
        default:
            cout << "Domain type not specified for quadrature refinement" << endl;
            exit(0);
            break;
    }
    
//     long double phi_start, phi_end;
//     long double mu_start, mu_end;
//     if (refinement_level > -1) {
//         cout << "refinement_level = " << refinement_level << "  "  << "id_refinement_level = " << id_refinement_level << endl;
//         cout << "N_peaks_mu = " << N_peaks_mu << "  " << "N_peaks_phi = " << N_peaks_phi << endl;
//         for (int i_phi_peak = 0; i_phi_peak < N_peaks_phi; i_phi_peak++) {
//             cout << "i_phi_peak = " << i_phi_peak << "  " << "phi_peak = " << phi_peak[i_phi_peak] << endl;
//         }
//         
//         for (int i_mu_peak = 0; i_mu_peak < N_peaks_mu; i_mu_peak++) {
//             cout << "i_mu_peak = " << i_mu_peak << "  " << "mu_peak = " << mu_peak[i_mu_peak] << endl;
//         }
//         
//         for (int i_phi_quad = 0; i_phi_quad < N_intervals_phi_quad; i_phi_quad++) {
//             phi_start = bounds_quad_phi_refin[i_phi_quad];
//             phi_end = bounds_quad_phi_refin[i_phi_quad+1];
//             cout << "i_phi_quad = " << i_phi_quad << "  " << "N_intervals_phi_quad = " << N_intervals_phi_quad << "  " << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << endl;
//         }
//         
//         for (int i_mu_quad = 0; i_mu_quad < N_intervals_mu_quad; i_mu_quad++) {
//             mu_start = bounds_quad_mu_refin[i_mu_quad];
//             mu_end = bounds_quad_mu_refin[i_mu_quad+1];
//             cout << "i_mu_quad = " << i_mu_quad << "  " << "N_intervals_mu_quad = " << N_intervals_mu_quad << "  " << "mu_start = " << mu_start << "  " << "mu_end = " << mu_end << "  " << endl;
//         }
//     }
}

inline void M2_State_Param :: compute_Omegas(const long double &mu_start, const long double &mu_end, const long double &phi_start, const long double &phi_end) {
    int index;
    long double mu_temp, phi_temp;
    long double w_phi_refin, w_mu_refin;
    switch (Domain_Type) {    
        case BOUNDARY_GAM1:
            for (int Id_Phi = 0; Id_Phi < 2*Order_Quad_mu; Id_Phi++) {
                mu_temp = N1_1;
                // phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], 0.0, 2.0*PI);
                phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], -1.0, 1.0);
                Omega1[Id_Phi] = mu_temp;
                Omega2[Id_Phi] = sqrt(1.0 - pow(mu_temp,2))*cos(phi_temp);
                Omega3[Id_Phi] = sqrt(1.0 - pow(mu_temp,2))*sin(phi_temp);
                
                // w_phi_refin = diff_linear_interpolation(phi_start, phi_end, 0.0, 2.0*PI);
                w_phi_refin = diff_linear_interpolation(phi_start, phi_end, -1.0, 1.0);
                
                if (w_phi_refin < 0.0) {
                    cout << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << "w_phi_refin = " << w_phi_refin << endl;
                    exit(0);
                }
                
                weight_total[Id_Phi] = /*2.0 * PI **/ w_phi_refin * w_phi_quad[Id_Phi];
            }
            break;    
        case BOUNDARY_GAM2:
            for (int Id_Phi = 0; Id_Phi < 2*Order_Quad_mu; Id_Phi++) {
                mu_temp = N1_2;
                // phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], 0.0, 2.0*PI);
                phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], -1.0, 1.0);
                Omega2[Id_Phi] = mu_temp;
                Omega3[Id_Phi] = sqrt(1.0 - pow(mu_temp,2))*cos(phi_temp);
                Omega1[Id_Phi] = sqrt(1.0 - pow(mu_temp,2))*sin(phi_temp);
                
                // w_phi_refin = diff_linear_interpolation(phi_start, phi_end, 0.0, 2.0*PI);
                w_phi_refin = diff_linear_interpolation(phi_start, phi_end, -1.0, 1.0);
                
                if (w_phi_refin < 0.0) {
                    cout << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << "w_phi_refin = " << w_phi_refin << endl;
                    exit(0);
                }
                
                weight_total[Id_Phi] = /*2.0 * PI * */ w_phi_refin * w_phi_quad[Id_Phi];
            }
            break;    
        case BOUNDARY_GAM3:
            for (int Id_Phi = 0; Id_Phi < 2*Order_Quad_mu; Id_Phi++) {
                mu_temp = N1_3;
                // phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], 0.0, 2.0*PI);
                phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], -1.0, 1.0);
                Omega3[Id_Phi] = mu_temp;
                Omega1[Id_Phi] = sqrt(1.0 - pow(mu_temp,2))*cos(phi_temp);
                Omega2[Id_Phi] = sqrt(1.0 - pow(mu_temp,2))*sin(phi_temp);
                
                // w_phi_refin = diff_linear_interpolation(phi_start, phi_end, 0.0, 2.0*PI);
                w_phi_refin = diff_linear_interpolation(phi_start, phi_end, -1.0, 1.0);
                
                if (w_phi_refin < 0.0) {
                    cout << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << "w_phi_refin = " << w_phi_refin << endl;
                    exit(0);
                }
                
                weight_total[Id_Phi] = /*2.0 * PI **/ w_phi_refin * w_phi_quad[Id_Phi];
            }
            break;  
        default:
            for (int Id_Phi = 0; Id_Phi < 2*Order_Quad_mu; Id_Phi++) {
                for (int Id_mu = 0; Id_mu < Order_Quad_mu; Id_mu++) {
                    index = Id_Angle(Id_mu, Id_Phi);
                    
                    // phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], 0.0, 2.0*PI);
                    phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], -1.0, 1.0);
                    mu_temp = linear_interpolation(mu_start, mu_end, mu_quad[Id_mu], -1.0, 1.0);
                    
                    // w_phi_refin = diff_linear_interpolation(phi_start, phi_end, 0.0, 2.0*PI);
                    w_phi_refin = diff_linear_interpolation(phi_start, phi_end, -1.0, 1.0);
                    
                    if (w_phi_refin < 0.0) {
                        cout << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << "w_phi_refin = " << w_phi_refin << endl;
                        exit(0);
                    }
                    
                    w_mu_refin = diff_linear_interpolation(mu_start, mu_end, -1.0, 1.0);
                    
                    if (w_mu_refin < 0.0) {
                        cout << "w_mu_refin = " << w_mu_refin << endl;
                        exit(0);
                    }
                    
                    // cout << "w_mu_refin = " << w_mu_refin << "  " << "w_phi_refin = " << w_phi_refin << "  " << "phi_temp = " << phi_temp << "  " << "mu_temp = " << mu_temp << "  "  << "mu_quad[Id_mu] = " << mu_quad[Id_mu] << endl;
                    
                    Omega1[index] = mu_temp;
                    Omega2[index] = sqrt(1.0 - pow(mu_temp,2))*cos(phi_temp);
                    Omega3[index] = sqrt(1.0 - pow(mu_temp,2))*sin(phi_temp);
                    
                    weight_total[index] = /*2.0 * PI **/ w_phi_refin * w_phi_quad[Id_Phi];
                    weight_total[index] *= w_mu_refin * w_mu_quad[Id_mu];
                }
            }
            break;
    };
}

inline void M2_State_Param :: Allocate() {
    Deallocate();
    
    if (Index == NULL) {
        Index = new int [NVARS];
    }
    if (MOM_VEC == NULL) {
        MOM_VEC = new long double [NVARS];
    }
    if (MOM_VEC_iso == NULL) {
        MOM_VEC_iso = new long double [NVARS];
    }
    if (Sk == NULL) {
        Sk = new long double [NVARS*NVARS];
    }
    if (Sk_cur == NULL) {
        Sk_cur = new long double [NVARS*NVARS];
    }
    if (ak == NULL) {
        ak = new long double [NVARS*NVARS];
    }
    if (x == NULL) {
        x = new long double [NVARS];
    }
    if (tol_x == NULL) {
        tol_x = new long double [NVARS];
    }
    
    Allocate_Quad();
}

inline void M2_State_Param :: Allocate_Quad() {
    Deallocate_Quad();
    
    // Weight and abscissas for Lebedev quadrature
    if (Omega1 == NULL) {
        Omega1 = new long double[N_dirs()];
    }
    if (Omega2 == NULL) {
        Omega2 = new long double[N_dirs()];
    }
    if (Omega3 == NULL) {
        Omega3 = new long double[N_dirs()];
    }
    if (weight_total == NULL) {
        weight_total = new long double[N_dirs()];
    }
    
    if (phi_quad == NULL) {
        phi_quad = new long double[2*Order_Quad_mu];
    }
    if (w_phi_quad == NULL) {
        w_phi_quad = new long double[2*Order_Quad_mu];
    }
    
    if (mu_quad == NULL) {
        mu_quad = new long double[Order_Quad_mu];
    }
    if (w_mu_quad == NULL) {
        w_mu_quad = new long double[Order_Quad_mu];
    }

    // circle_rule ( 2*Order_Quad_mu, w_phi_quad, phi_quad );
    lobatto_compute ( 2*Order_Quad_mu, phi_quad, w_phi_quad );
    lobatto_compute ( Order_Quad_mu, mu_quad, w_mu_quad );
    
    // half_circle_rule( Order_Quad_mu, w_mu_quad, mu_quad );
    // legendre_set ( Order_Quad_mu, mu_quad, w_mu_quad );
}

inline void M2_State_Param :: Deallocate() {
    if (Index != NULL) {
        delete[] Index; Index = NULL;
    }
    if (MOM_VEC != NULL) {
        delete[] MOM_VEC; MOM_VEC = NULL;
    }
    if (MOM_VEC_iso != NULL) {
        delete[] MOM_VEC_iso; MOM_VEC_iso = NULL;
    }
    if (Sk != NULL) {
        delete[] Sk; Sk = NULL;
    }
    if (Sk_cur != NULL) {
        delete[] Sk_cur; Sk_cur = NULL;
    }
    if (ak != NULL) {
        delete[] ak; ak = NULL;
    }
    if (x != NULL) {
        delete[] x; x = NULL;
    }
    if (tol_x != NULL) {
        delete[] tol_x; tol_x = NULL;
    }
}

inline void M2_State_Param :: Deallocate_Quad() {
    if (Omega1 != NULL) {
        delete[] Omega1; Omega1 = NULL;
    }
    if (Omega2 != NULL) {
        delete[] Omega2; Omega2 = NULL;
    }
    if (Omega3 != NULL) {
        delete[] Omega3; Omega3 = NULL;
    }
    if (weight_total != NULL) {
        delete[] weight_total; weight_total = NULL;
    }
    if (phi_quad != NULL) {
        delete[] phi_quad; phi_quad = NULL;
    }
    if (w_phi_quad != NULL) {
        delete[] w_phi_quad; w_phi_quad = NULL;
    }
    if (mu_quad != NULL) {
        delete[] mu_quad; mu_quad = NULL;
    }
    if (w_mu_quad != NULL) {
        delete[] w_mu_quad; w_mu_quad = NULL;
    }
}

inline void M2_State_Param :: Set_Num_Vars(void) {
    switch (Dimension) {
        case ONE_DIMENSIONAL:
            NVARS = 3;
            break;
        case TWO_DIMENSIONAL:
            switch (Domain_Type){
                case BOUNDARY_GAM1:
                    NVARS = 3;
                    break;
                case BOUNDARY_GAM2:
                    NVARS = 3;
                    break;
                case BOUNDARY_GAM3:
                    NVARS = 3;
                    break;
                case GAM1_GAM2_GAM3:
                case BOUNDARY_GAM1_GAM2_EQ_0:
                case BOUNDARY_GAM1_GAM3_EQ_0:
                case BOUNDARY_GAM2_GAM3_EQ_0:
                case BOUNDARY_FREE_STREAMING:
                    NVARS = 6;
                    break;
                default:
                    NVARS = 6;
                    break;
            };
            break;
        case THREE_DIMENSIONAL:
            switch (Domain_Type){
                case BOUNDARY_GAM1:
                    NVARS = 5;
                    break;
                case BOUNDARY_GAM2:
                    NVARS = 5; 
                    break;
                case BOUNDARY_GAM3:
                    NVARS = 5;  
                    break;
                case GAM1_GAM2_GAM3:
                case BOUNDARY_GAM1_GAM2_EQ_0:
                case BOUNDARY_GAM1_GAM3_EQ_0:
                case BOUNDARY_GAM2_GAM3_EQ_0:
                case BOUNDARY_FREE_STREAMING:
                    NVARS = 9;
                    break;
                default:
                    NVARS = 9;
                    break;
            };
            break;
        default:
            NVARS = 9;
            break;
    };
}

inline int M2_State_Param :: Get_Num_Var(const int &Dimension, const int &Domain_Type) {
    int nvars;
    switch (Dimension) {
        case ONE_DIMENSIONAL:
            nvars = 3;
            break;
        case TWO_DIMENSIONAL:
            switch (Domain_Type){
                case BOUNDARY_GAM1:
                    nvars = 3;
                    break;
                case BOUNDARY_GAM2:
                    nvars = 3;
                    break;
                case BOUNDARY_GAM3:
                    nvars = 3;
                    break;
                case GAM1_GAM2_GAM3:
                case BOUNDARY_GAM1_GAM2_EQ_0:
                case BOUNDARY_GAM1_GAM3_EQ_0:
                case BOUNDARY_GAM2_GAM3_EQ_0:
                case BOUNDARY_FREE_STREAMING:
                    nvars = 6;
                    break;
                default:
                    nvars = 6;
                    break;
            };
            break;
        case THREE_DIMENSIONAL:
            switch (Domain_Type){
                case BOUNDARY_GAM1:
                    nvars = 5;
                    break;
                case BOUNDARY_GAM2:
                    nvars = 5; 
                    break;
                case BOUNDARY_GAM3:
                    nvars = 5;  
                    break;
                case GAM1_GAM2_GAM3:
                case BOUNDARY_GAM1_GAM2_EQ_0:
                case BOUNDARY_GAM1_GAM3_EQ_0:
                case BOUNDARY_GAM2_GAM3_EQ_0:
                case BOUNDARY_FREE_STREAMING:
                    nvars = 9;
                    break;
                default:
                    nvars = 9;
                    break;
            };
            break;
        default:
            nvars = 9;
            break;
    };
    return nvars;
}

inline void M2_State_Param :: Set_Indices(void) {
    //     if gam1 = 0, then moment that depend on Omega1 are not included
    //     Similarly, if gam2 = 0, then moment that depend on Omega2 are not included
    //     and, if gam3 = 0, then moment that depend on Omega3 are not included
    
    for (int i = 0; i < NVARS; i++) {
        Index[i] = i;
    }
    
    switch (Dimension) {
        case ONE_DIMENSIONAL:
            Index[0] = 0;
            Index[1] = 1;
            Index[2] = 4;
            break;
        case TWO_DIMENSIONAL:
            switch (Domain_Type){
                case BOUNDARY_GAM1:
                    Index[0] = 0;
                    Index[1] = 2; 
                    Index[2] = 7;
                    break;
                case BOUNDARY_GAM2:
                    Index[0] = 0;
                    Index[1] = 1;
                    Index[2] = 4;
                    break;
                case BOUNDARY_GAM3:
                    Index[0] = 0;
                    Index[1] = 1;
                    Index[2] = 4;
                    break;
                case GAM1_GAM2_GAM3:
                case BOUNDARY_GAM1_GAM2_EQ_0:
                case BOUNDARY_GAM1_GAM3_EQ_0:
                case BOUNDARY_GAM2_GAM3_EQ_0:
                case BOUNDARY_FREE_STREAMING:
                    Index[0] = 0;
                    Index[1] = 1;
                    Index[2] = 2;
                    Index[3] = 4; 
                    Index[4] = 5; 
                    Index[5] = 7; 
                    break;
                default:
                    for (int i = 0; i < NVARS; i++)
                        Index[i] = i;
                    break;
            };
            break;
        case THREE_DIMENSIONAL:
            switch (Domain_Type){
                case BOUNDARY_GAM1:
                    Index[0] = 0;
                    Index[1] = 2;
                    Index[2] = 3; 
                    Index[3] = 7; 
                    Index[4] = 8;    
                    break;
                case BOUNDARY_GAM2:
                    Index[0] = 0;
                    Index[1] = 1;
                    Index[2] = 3; 
                    Index[3] = 4; 
                    Index[4] = 6;    
                    break;
                case BOUNDARY_GAM3:
                    Index[0] = 0;
                    Index[1] = 1;
                    Index[2] = 2; 
                    Index[3] = 4; 
                    Index[4] = 5;    
                    break;
                case GAM1_GAM2_GAM3:
                case BOUNDARY_FREE_STREAMING:
                case BOUNDARY_GAM1_GAM2_EQ_0:
                case BOUNDARY_GAM1_GAM3_EQ_0:
                case BOUNDARY_GAM2_GAM3_EQ_0:
                    Index[0] = 0;
                    Index[1] = 1;
                    Index[2] = 2;
                    Index[3] = 3; 
                    Index[4] = 4; 
                    Index[5] = 5; 
                    Index[6] = 6; 
                    Index[7] = 7; 
                    Index[8] = 8;
                    break;
                default:
                    cout << "******************** No Domain Type Found !!!!!!! *******************" << endl;
                    exit(0);
                    break;
            };
            break;
        default:
            cout << "******************** No Dimension Type Found !!!!!!! *******************" << endl;
            exit(0);
            break;
    };
}

#endif // _M2_STATE_PARAMETERS_H_INCLUDED
