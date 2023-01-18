#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#include "NG_MN_Model_3D_OPTIM.h"
#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED

//**************************************************************************************
// This routines computes the third-order closing fluxes in the case of the free-
// streaming limit
//**************************************************************************************
void Set_Closure_Free_Streaming(record_N3 *rec_N3_local, const int &id_count, M2_State_Param &M2_State) {
    rec_N3_local[id_count].ratio_I0 = M2_State.ratio_I0;
    rec_N3_local[id_count].I0 = M2_State.I0;
    rec_N3_local[id_count].N1 = M2_State.N1;
    rec_N3_local[id_count].N1_1 = M2_State.N1_1;
    rec_N3_local[id_count].N1_2 = M2_State.N1_2;
    rec_N3_local[id_count].N1_3 = M2_State.N1_3;
    rec_N3_local[id_count].gam1 = M2_State.gamma_1;
    rec_N3_local[id_count].gam2 = M2_State.gamma_2;
    
    rec_N3_local[id_count].x_SH = M2_State.x_SH;
    rec_N3_local[id_count].y_SH = M2_State.y_SH;
    rec_N3_local[id_count].z_SH = M2_State.z_SH;
    
    rec_N3_local[id_count].N3_111 = pow(M2_State.N1_1, 3);
    rec_N3_local[id_count].N3_122 = M2_State.N1_1*pow(M2_State.N1_2, 2);
    rec_N3_local[id_count].N3_123 = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
    
    rec_N3_local[id_count].dN3_111_dN1_1 = 3.0*pow(M2_State.N1_1, 2);
    rec_N3_local[id_count].dN3_122_dN1_1 = pow(M2_State.N1_2, 2);
}

//**************************************************************************************
// This routines computes the third-order closing fluxes for the M2 closure at one of the
// vertices of the triangle (P1 P2 P3)
//**************************************************************************************
void Set_Closure_Triangle_Vertices(record_N3 *rec_N3_local, const int &id_count, M2_State_Param &M2_State) {
    long double norm_f_2;
    rec_N3_local[id_count].ratio_I0 = M2_State.ratio_I0;
    rec_N3_local[id_count].I0 = M2_State.I0;
    rec_N3_local[id_count].N1 = M2_State.N1;
    rec_N3_local[id_count].N1_1 = M2_State.N1_1;
    rec_N3_local[id_count].N1_2 = M2_State.N1_2;
    rec_N3_local[id_count].N1_3 = M2_State.N1_3;
    rec_N3_local[id_count].gam1 = M2_State.gamma_1;
    rec_N3_local[id_count].gam2 = M2_State.gamma_2;
    
    rec_N3_local[id_count].x_SH = M2_State.x_SH;
    rec_N3_local[id_count].y_SH = M2_State.y_SH;
    rec_N3_local[id_count].z_SH = M2_State.z_SH;
    
    norm_f_2 = pow(M2_State.N1_1, 2) + pow(M2_State.N1_2, 2) + pow(M2_State.N1_3, 2);
    
    rec_N3_local[id_count].N3_111 = M2_State.N1_1*(pow(M2_State.N1_1, 2) + M2_State.gamma_1*(1.0 - norm_f_2));
    rec_N3_local[id_count].N3_122 = M2_State.N1_1*(pow(M2_State.N1_2, 2) + M2_State.gamma_2*(1.0 - norm_f_2));
    rec_N3_local[id_count].N3_123 = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
    
    rec_N3_local[id_count].dN3_111_dN1_1 = 3.0*pow(M2_State.N1_1, 2) + M2_State.gamma_1*(1.0 - norm_f_2 - 2.0*pow(M2_State.N1_1, 2));
    // rec_N3_local[id_count].dN3_111_dN1_2 = -2.0*M2_State.gamma_1*M2_State.N1_1*M2_State.N1_2;
    // rec_N3_local[id_count].dN3_111_dN1_3 = -2.0*M2_State.gamma_1*M2_State.N1_1*M2_State.N1_3;
    
    rec_N3_local[id_count].dN3_122_dN1_1 = pow(M2_State.N1_2, 2) + M2_State.gamma_2*(1.0 - norm_f_2 - 2.0*pow(M2_State.N1_1, 2));
    // rec_N3_local[id_count].dN3_122_dN1_2 = 2.0*M2_State.N1_1*M2_State.N1_2*(1.0 - M2_State.gamma_2);
    // rec_N3_local[id_count].dN3_122_dN1_3 = -2.0*M2_State.gamma_2*M2_State.N1_1*M2_State.N1_3;
            
    rec_N3_local[id_count].dN3_123_dN1_1 = M2_State.N1_2*M2_State.N1_3;
    rec_N3_local[id_count].dN3_123_dN1_2 = M2_State.N1_1*M2_State.N1_3;
    rec_N3_local[id_count].dN3_123_dN1_3 = M2_State.N1_1*M2_State.N1_2;
}

//**************************************************************************************
// This routines computes the third-order closing fluxes based on the regime of radiation
// encountered
//**************************************************************************************
void Compute_Third_Order_Closing_Fluxes(record_N3 *rec_N3_local, const int &index_count, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, M2_State_Param &M2_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, long double *Sk_final) {
    long double N3[3];
    Finite_Diff_Parameters Finite_Diff_N3;
    
    if (M2_State.Domain_Type == BOUNDARY_FREE_STREAMING) {
        Set_Closure_Free_Streaming(rec_N3_local, index_count, M2_State);
    } else if ((M2_State.Domain_Type == BOUNDARY_GAM1_GAM2_EQ_0) || 
               (M2_State.Domain_Type == BOUNDARY_GAM1_GAM3_EQ_0) || 
               (M2_State.Domain_Type == BOUNDARY_GAM2_GAM3_EQ_0)) {
        // Corresponds to one of the vertices of the triangel P1 P2 P3
        Set_Closure_Triangle_Vertices(rec_N3_local, index_count, M2_State);
    } else if (M2_State.flag_Taylor_Series_Expansion_Sphere && M2_State.flag_finite_diff_Sphere == false) {
        rec_N3_local[index_count].ratio_I0 = M2_State.ratio_I0;
        rec_N3_local[index_count].I0 = M2_State.I0;
        rec_N3_local[index_count].N1 = M2_State.N1;
        rec_N3_local[index_count].N1_1 = M2_State.N1_1;
        rec_N3_local[index_count].N1_2 = M2_State.N1_2;
        rec_N3_local[index_count].N1_3 = M2_State.N1_3;
        rec_N3_local[index_count].gam1 = M2_State.gamma_1;
        rec_N3_local[index_count].gam2 = M2_State.gamma_2;
        
        rec_N3_local[index_count].x_SH = M2_State.x_SH;
        rec_N3_local[index_count].y_SH = M2_State.y_SH;
        rec_N3_local[index_count].z_SH = M2_State.z_SH;
        
        // In this case we perform Taylor series expansion in the vicinity of the free-streaming 
        // limit to compute the Eddington factor
        Setup_Taylor_Series_Coefficients(Finite_Diff_N3, M2_State, Var_index, num_points, M2_State.id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
        
        Finite_Diff_N3.Taylor_Series_N3_ijk(rec_N3_local[index_count].N3_111, rec_N3_local[index_count].N3_122, rec_N3_local[index_count].N3_123, M2_State.N1);
    } else if (M2_State.flag_Taylor_Series_Expansion_Triangle && M2_State.flag_finite_diff_Triangle == false) {
        rec_N3_local[index_count].ratio_I0 = M2_State.ratio_I0;
        rec_N3_local[index_count].I0 = M2_State.I0;
        rec_N3_local[index_count].N1 = M2_State.N1;
        rec_N3_local[index_count].N1_1 = M2_State.N1_1;
        rec_N3_local[index_count].N1_2 = M2_State.N1_2;
        rec_N3_local[index_count].N1_3 = M2_State.N1_3;
        rec_N3_local[index_count].gam1 = M2_State.gamma_1;
        rec_N3_local[index_count].gam2 = M2_State.gamma_2;
        
        rec_N3_local[index_count].x_SH = M2_State.x_SH;
        rec_N3_local[index_count].y_SH = M2_State.y_SH;
        rec_N3_local[index_count].z_SH = M2_State.z_SH;
        
        // In this case we perform Taylor series expansion in the vicinity of the free-streaming 
        // limit to compute the Eddington factor
        Setup_Taylor_Series_Coefficients_Triangle(Finite_Diff_N3, M2_State, Var_index, num_points, M2_State.id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
        
        Finite_Diff_N3.Taylor_Series_N3_ijk_Triangle(rec_N3_local[index_count].N3_111, rec_N3_local[index_count].N3_122, rec_N3_local[index_count].N3_123);
    } else {
        N3_M2_3D_RT(N3, x, Sk_final, M2_State);
        
        if (M2_State.Regime == HYPERBOLIC_LIMIT) {
            M2_State.MOM_VEC[0] = 0.0;
        } else if (M2_State.Regime == LOGARITHMIC_LIMIT) {
            M2_State.MOM_VEC[0] = 1.0e32;
        }
        
        rec_N3_local[index_count].ratio_I0 = M2_State.ratio_I0;
        rec_N3_local[index_count].I0 = M2_State.I0;
        rec_N3_local[index_count].N1 = M2_State.N1;
        rec_N3_local[index_count].N1_1 = M2_State.N1_1;
        rec_N3_local[index_count].N1_2 = M2_State.N1_2;
        rec_N3_local[index_count].N1_3 = M2_State.N1_3;
        rec_N3_local[index_count].gam1 = M2_State.gamma_1;
        rec_N3_local[index_count].gam2 = M2_State.gamma_2;
        
        rec_N3_local[index_count].x_SH = M2_State.x_SH;
        rec_N3_local[index_count].y_SH = M2_State.y_SH;
        rec_N3_local[index_count].z_SH = M2_State.z_SH;
        
        rec_N3_local[index_count].N3_111 = N3[0];
        rec_N3_local[index_count].N3_112 = N3[1];
        rec_N3_local[index_count].N3_122 = N3[2];
        rec_N3_local[index_count].N3_123 = N3[3];
        rec_N3_local[index_count].N3_222 = N3[4];
        
        if (M2_State.display) {
            if (M2_State.id_proc == M2_State.proc_display) {
                cout << "N3_111 = " << N3[0] << "  "  << "N3_112 = " << N3[1] << "  "  << "N3_122 = " << N3[2] << "  "  << "N3_123 = " << N3[3] << "  "  << "N3_122 = " << N3[4] << endl;
                cout << "\n" << endl;
            }
        }
        
    }
}

//**************************************************************************************
// This routines computes the derivatives of the third-order closing fluxes based on the 
// regime of radiation encountered
//**************************************************************************************
void Compute_Third_Order_Closing_Fluxes_Derivatives(record_N3 *rec_N3_local, const int &index_count, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, M2_State_Param &M2_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, long double *Sk_final) {
    
    if (M2_State.Domain_Type == BOUNDARY_FREE_STREAMING) {
        // Free-streaming limit
        if (!M2_State.flag_finite_diff_Sphere) {
            // flag_finite_diff_Sphere = true means that a finite differencing procedure for the free-streaming limit has already been initiated, so no need to start another one
            // flag_finite_diff_Sphere = false means that there is no finite differencing procedure currently being carried out. 
            Calculate_Higher_Order_Moments_Derivatives_Boundary_Norm_f(rec_N3_local[index_count], M2_State, Var_index, num_points, M2_State.id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
        }
        cout << "exit for now!!!!!!!!" << endl;
        exit(0);
    } else if ((M2_State.Domain_Type == BOUNDARY_GAM1_GAM2_EQ_0) || 
               (M2_State.Domain_Type == BOUNDARY_GAM1_GAM3_EQ_0) || 
               (M2_State.Domain_Type == BOUNDARY_GAM2_GAM3_EQ_0)) {
                       
        if (!M2_State.flag_finite_diff_Triangle) {
            // flag_finite_diff_Triangle = true means that a finite differencing procedure on one of the vertices of the triangle P1 P2 P3 has already been initiated, so no need to start another one
            // flag_finite_diff_Triangle = false means that there is no finite differencing procedure currently being carried out on any of the vertices of the triangle P1 P2 P3.
            Calculate_Higher_Order_Moments_Derivatives_Triangle_Vertices(rec_N3_local[index_count], M2_State, Var_index, num_points, M2_State.id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, M2_State.flag_finite_diff_Sphere);
        }
        cout << "exit for now!!!!!!!!" << endl;
        exit(0);
    } else if ((M2_State.Domain_Type == BOUNDARY_GAM1) || 
               (M2_State.Domain_Type == BOUNDARY_GAM2) || 
               (M2_State.Domain_Type == BOUNDARY_GAM3)) {
            // This means the current set of angular moments up to second-order of
            // interest is either one of the domains BOUNDARY_GAM1, BOUNDARY_GAM2,
            // or BOUNDARY_GAM3
                   
            if (!M2_State.flag_finite_diff_Triangle) {
                // flag_finite_diff_Triangle = true means that a finite differencing procedure on one of the edges of the triangle P1 P2 P3 has already been initiated, so no need to start another one
                // flag_finite_diff_Triangle = false means that there is no finite differencing procedure currently being carried out on any of the vertices of the triangle P1 P2 P3. 
                // In this case, if Domain_Type is BOUNDARY_GAM1 or BOUNDARY_GAM2 or BOUNDARY_GAM3, then 
                // we use a finite difference procedure to compute the derivatives at either one of the edges
                Calculate_Higher_Order_Moments_Derivatives_Boundary_gam1_gam2_gam3(rec_N3_local[index_count], M2_State, Var_index, num_points, M2_State.id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, M2_State.flag_finite_diff_Sphere);
            }
        cout << "exit for now!!!!!!!!" << endl;
        exit(0);
    } else {
        // This means the current set of angular moments up to second-order of
        // interest is in the domain GAM1_GAM2_GAM3
        Calculate_Higher_Order_Moments_Derivatives(rec_N3_local[index_count], x, Sk_final, M2_State);
    }
}

//**************************************************************************************
// This routines performs finite differencing in either the Hyperbolic limit or the 
// logarithmic limit in order to compute derivatives of the third-order closing fluxes 
// with respect to the (exponential) mapping of the radiative energy density
//**************************************************************************************
void Calculate_Higher_Order_Moments_Derivatives_Boundary_Spectrum_Energy(record_N3 &rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed) {
    // Perform finite difference approximation to obtain derivative
    M2_State_Param M2_State_Temp(M2_State);
    record_N3 *rec_N3_local_diff;
    long double h;
    int order, prec, n;
    long double *c, *x;
    long double d_N3_111_dnorm_f, d_N3_122_dnorm_f, d_N3_123_dnorm_f;
    long double d2_N3_111_dnorm_f_dN1_1, d2_N3_122_dnorm_f_dN1_1;
    order = 1;
    prec = 2;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    rec_N3_local_diff = new record_N3[n];
    h = finite_diff_h_N1;
    
    // Compute values for performing finite difference
    M2_State_Temp.finite_diff_domain_Spectrum = SPECTRUM_ENERGY_FINITE_DIFFERENCE;
    M2_State_Temp.h_ratio_I0 = h;
    
    switch (M2_State_Temp.Finite_Diff_Spectrum_Boundary_Type) {
        case HYPERBOLIC_LIMIT:
            M2_State_Temp.ratio_I0_knot = -1.0;
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_ratio_I0, order, prec, c, x );
            break;
        case LOGARITHMIC_LIMIT:
            M2_State_Temp.ratio_I0_knot = 1.0;
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_ratio_I0, order, prec, c, x );
            break;
        default:
            cout << "Finite_Diff_Spectrum_Boundary_Type not specified !!!!" << endl;
            exit(0);
            break;
    }
    
    M2_State_Temp.flag_finite_diff_Sphere = false;
    M2_State_Temp.flag_finite_diff_Triangle = false;
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M2_State_Temp.id_finite_diff_ratio_I0 = id_diff;
        M2_State_Temp.x_finite_diff_ratio_I0 = x[id_diff];
        
        NLOPT_Optim_Algo(rec_N3_local_diff, id_diff, M2_State_Temp, Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
    }
    
    for (int i = 0; i < n; i++) {
        
    }
    
    if (M2_State_Temp.display) {
        // cout << "dChi2_drI0 = " << rec_Chi2.dChi2_drI0 << "   " << "I0 = " << M1_State.I0 << "   " << "ratio_I0 = " << M1_State.ratio_I0 << endl;
    }
    
    delete[] c;
    delete[] x;
    delete[] rec_N3_local_diff;
}

//**************************************************************************************
// This routines computes the coefficients required for performing a taylor series expansion 
// in the vicinity of the free-streaming limit along the radius of the unit ball defined
// by the realizable space for the first-order angular moments in order to compute 
// the third-order closing fluxes near such limit
//**************************************************************************************
void Setup_Taylor_Series_Coefficients(Finite_Diff_Parameters &Finite_Diff_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed) {
    // Perform finite difference approximation to obtain derivative
    M2_State_Param M2_State_Temp(M2_State);
    record_N3 *rec_N3_local_diff;
    long double h;
    int order, prec, n, nmax;
    long double *c, *x;
    long double d_N3_111, d_N3_122, d_N3_123;
    long double d_N3_111_dnorm_f, d_N3_122_dnorm_f, d_N3_123_dnorm_f;
    long double d2_N3_111_dnorm_f_dN1_1, d2_N3_122_dnorm_f_dN1_1;
    int Domain_Type;
    order = 4;
    prec = 2;
    n = order + prec;
    nmax = n;
    c = new long double[n];
    x = new long double[n];
    rec_N3_local_diff = new record_N3[n];
    h = finite_diff_h_N1;
    
    long double N1_1_proj, N1_2_proj, N1_3_proj;
    // Projections of the first-order angular moments on the surface of the sphere 
    N1_1_proj = M2_State.N1_1/M2_State.N1;
    N1_2_proj = M2_State.N1_2/M2_State.N1;
    N1_3_proj = M2_State.N1_3/M2_State.N1;
    
    Finite_Diff_N3.x0 = 1.0;
    Finite_Diff_N3.N3_111_x0 = pow(N1_1_proj,3);
    Finite_Diff_N3.N3_122_x0 = N1_1_proj*pow(N1_2_proj, 2);
    Finite_Diff_N3.N3_123_x0 = N1_1_proj*N1_2_proj*N1_3_proj;
    
    // Compute values for performing finite difference
    M2_State_Temp.finite_diff_domain_Sphere = RADIUS_SPHERE_FINITE_DIFFERENCE;
    M2_State_Temp.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
    M2_State_Temp.h_norm_f = h;
    
    // backward finite difference approximation to fourth derivative
    differ_backward ( M2_State_Temp.h_norm_f, order, prec, c, x );
    
    M2_State_Temp.norm_f_knot = 1.0;
    
    M2_State_Temp.flag_finite_diff_Sphere = true;
    M2_State_Temp.flag_finite_diff_Triangle = false;
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M2_State_Temp.id_finite_diff_norm_f = id_diff;
        M2_State_Temp.x_finite_diff_norm_f = x[id_diff];
        
        Check_Domain_Type(Domain_Type, Var_index, num_points, M2_State_Temp);
        if ((Domain_Type != BOUNDARY_GAM1_GAM2_EQ_0) || (Domain_Type != BOUNDARY_GAM1_GAM3_EQ_0) || (Domain_Type != BOUNDARY_GAM2_GAM3_EQ_0)) {
            if (M2_State_Temp.Domain_Type != Domain_Type) {
                M2_State_Temp.Reset_Num_Vars(Domain_Type);
            }
        }
        
        NLOPT_Optim_Algo(rec_N3_local_diff, id_diff, M2_State_Temp, Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
    }
    
    d_N3_111 = 0.0;
    d_N3_122 = 0.0;
    d_N3_123 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_111;
        d_N3_122 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_122;
        d_N3_123 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_123;
    }
    
    Finite_Diff_N3.d4N3_111_dN1 = d_N3_111;
    Finite_Diff_N3.d4N3_122_dN1 = d_N3_122;
    Finite_Diff_N3.d4N3_123_dN1 = d_N3_123;
    
    // backward finite difference approximation to third derivative
    order = 3;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    
    d_N3_111 = 0.0;
    d_N3_122 = 0.0;
    d_N3_123 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_111;
        d_N3_122 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_122;
        d_N3_123 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_123;
    }
    
    Finite_Diff_N3.d3N3_111_dN1 = d_N3_111;
    Finite_Diff_N3.d3N3_122_dN1 = d_N3_122;
    Finite_Diff_N3.d3N3_123_dN1 = d_N3_123;
    
    // backward finite difference approximation to second derivative
    order = 2;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    
    d_N3_111 = 0.0;
    d_N3_122 = 0.0;
    d_N3_123 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_111;
        d_N3_122 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_122;
        d_N3_123 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_123;
    }
    
    Finite_Diff_N3.d2N3_111_dN1 = d_N3_111;
    Finite_Diff_N3.d2N3_122_dN1 = d_N3_122;
    Finite_Diff_N3.d2N3_123_dN1 = d_N3_123;
    
    // backward finite difference approximation to first derivative
    order = 1;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    
    d_N3_111 = 0.0;
    d_N3_122 = 0.0;
    d_N3_123 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_111;
        d_N3_122 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_122;
        d_N3_123 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_123;
    }
    
    Finite_Diff_N3.dN3_111_dN1 = d_N3_111;
    Finite_Diff_N3.dN3_122_dN1 = d_N3_122;
    Finite_Diff_N3.dN3_123_dN1 = d_N3_123;
    
    if (M2_State_Temp.display) {
        if (M2_State_Temp.id_proc == M2_State.proc_display) {
            for (int id_diff = 0; id_diff < nmax; id_diff++) {
                cout << "id_diff = " << id_diff << "   " << "x = " << x[id_diff] << "   " << "N3_111 = " << rec_N3_local_diff[id_diff].N3_111 << "   " << "N3_122 = " << rec_N3_local_diff[id_diff].N3_122 << "   " << "N3_123 = " << rec_N3_local_diff[id_diff].N3_123 << endl;
            }
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] rec_N3_local_diff;
}

void Setup_finite_diff_Parameters_Triangle_Taylor_Series(M2_State_Param &M2_State, Finite_Diff_Parameters &Finite_Diff_N3) {
    
    switch (M2_State.Domain_Type) {
        case GAM1_GAM2_GAM3:
            M2_State.finite_diff_domain_Triangle = MEDIAN_TRIANGLE_FINITE_DIFFERENCE;
            cout << "Fix this !!!!!!!!!!!!!!!!!!" << endl;
            // Setup_finite_diff_Triangle_Median(M2_State);
            
            Finite_Diff_N3.delta_gam_Triangle = pow(M2_State.gamma_1 - M2_State.gam1_knot, 2) + pow(M2_State.gamma_2 - M2_State.gam2_knot, 2);
            Finite_Diff_N3.delta_gam_Triangle /= pow(M2_State.gam1_bound_min - M2_State.gam1_bound_max, 2) + pow(M2_State.gam2_bound_min - M2_State.gam2_bound_max, 2);
            
            switch(M2_State.finite_diff_type) {
                case FORWARD_FINITE_DIFFERENCE:
                    break;
                case BACKWARD_FINITE_DIFFERENCE:
                    Finite_Diff_N3.delta_gam_Triangle = -Finite_Diff_N3.delta_gam_Triangle;
                    break;
                default:
                    cout << "Finite difference type not specified" << endl;
                    exit(0);
                    break;
            }
            break;
        case BOUNDARY_GAM1:
            M2_State.finite_diff_domain_Triangle = EDGE_TRIANGLE_FINITE_DIFFERENCE;
            // gam1 = 0 in this case
            if (M2_State.gamma_2 < finite_diff_h_gam) {
                M2_State.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
                M2_State.gam1_knot = 0.0;
                M2_State.gam2_knot = 0.0;
                
                // Point P3
                Finite_Diff_N3.N3_111_x0 = pow(M2_State.N1_1,3);
                Finite_Diff_N3.N3_122_x0 = M2_State.N1_1*pow(M2_State.N1_2,2);
                Finite_Diff_N3.N3_123_x0 = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
                
                Finite_Diff_N3.delta_gam_Triangle = M2_State.gamma_2 - M2_State.gam2_knot;
            } else if ((1.0 - M2_State.gamma_2) < finite_diff_h_gam) {
                M2_State.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
                M2_State.gam1_knot = 0.0;
                M2_State.gam2_knot = 1.0;
                
                // Point P2
                Finite_Diff_N3.N3_111_x0 = pow(M2_State.N1_1,3);
                Finite_Diff_N3.N3_122_x0 = M2_State.N1_1*(pow(M2_State.N1_2,2) + (1.0 - pow(M2_State.N1,2)));
                Finite_Diff_N3.N3_123_x0 = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
                
                Finite_Diff_N3.delta_gam_Triangle = M2_State.gamma_2 - M2_State.gam2_knot;
            } else {
                cout << "Inconsistency in Taylor series finite differencing over the triangle" << endl;
                exit(0);
            }
            break;
        case BOUNDARY_GAM2:
            M2_State.finite_diff_domain_Triangle = EDGE_TRIANGLE_FINITE_DIFFERENCE;
            // gam2 = 0 in this case
            if (M2_State.gamma_1 < finite_diff_h_gam) {
                M2_State.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
                M2_State.gam1_knot = 0.0;
                M2_State.gam2_knot = 0.0;
                
                // Point P3
                Finite_Diff_N3.N3_111_x0 = pow(M2_State.N1_1,3);
                Finite_Diff_N3.N3_122_x0 = M2_State.N1_1*pow(M2_State.N1_2,2);
                Finite_Diff_N3.N3_123_x0 = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
                
                Finite_Diff_N3.delta_gam_Triangle = M2_State.gamma_1 - M2_State.gam1_knot;
            } else if ((1.0 - M2_State.gamma_1) < finite_diff_h_gam) {
                M2_State.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
                M2_State.gam1_knot = 1.0;
                M2_State.gam2_knot = 0.0;
                
                // Point P1
                Finite_Diff_N3.N3_111_x0 = M2_State.N1_1*(pow(M2_State.N1_1,2) + (1.0 - pow(M2_State.N1,2)));
                Finite_Diff_N3.N3_122_x0 = M2_State.N1_1*pow(M2_State.N1_2,2);
                Finite_Diff_N3.N3_123_x0 = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
                
                Finite_Diff_N3.delta_gam_Triangle = M2_State.gamma_1 - M2_State.gam1_knot;
            } else {
                cout << "Inconsistency in Taylor series finite differencing over the triangle" << endl;
                exit(0);
            }
            break;
        case BOUNDARY_GAM3:
            M2_State.finite_diff_domain_Triangle = EDGE_TRIANGLE_FINITE_DIFFERENCE;
            // gam3 = 1 - gam1 - gam2 = 0 in this case
            if (M2_State.gamma_1 < finite_diff_h_gam) {
                M2_State.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
                M2_State.gam1_knot = 0.0;
                M2_State.gam2_knot = 1.0;
                
                // Point P2
                Finite_Diff_N3.N3_111_x0 = pow(M2_State.N1_1,3);
                Finite_Diff_N3.N3_122_x0 = M2_State.N1_1*(pow(M2_State.N1_2,2) + (1.0 - pow(M2_State.N1,2)));
                Finite_Diff_N3.N3_123_x0 = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
                
                Finite_Diff_N3.delta_gam_Triangle = pow(M2_State.gamma_1 - M2_State.gam1_knot, 2) + pow(M2_State.gamma_2 - M2_State.gam2_knot, 2);
                Finite_Diff_N3.delta_gam_Triangle /= sqrt(2.0);
            } else if ((1.0 - M2_State.gamma_1) < finite_diff_h_gam) {
                M2_State.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
                M2_State.gam1_knot = 1.0;
                M2_State.gam2_knot = 0.0;
                
                // Point P1
                Finite_Diff_N3.N3_111_x0 = M2_State.N1_1*(pow(M2_State.N1_1,2) + (1.0 - pow(M2_State.N1,2)));
                Finite_Diff_N3.N3_122_x0 = M2_State.N1_1*pow(M2_State.N1_2,2);
                Finite_Diff_N3.N3_123_x0 = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
                
                Finite_Diff_N3.delta_gam_Triangle = pow(M2_State.gamma_1 - M2_State.gam1_knot, 2) + pow(M2_State.gamma_2 - M2_State.gam2_knot, 2);
                Finite_Diff_N3.delta_gam_Triangle /= -sqrt(2.0);
                // The minus sign accounts for the fact that we are going backward on the line segment
            } else {
                cout << "Inconsistency in Taylor series finite differencing over the triangle" << endl;
                exit(0);
            }
            break;
        default:
            cout << "Incorrect domain type for Taylor series finite differencing over the triangle" << endl;
            exit(0);
            break;
    }
}

//**************************************************************************************
// This routines computes the coefficients required for performing a taylor series expansion 
// in the vicinity of the free-streaming limit along the radius of the unit ball defined
// by the realizable space for the first-order angular moments in order to compute 
// the third-order closing fluxes near such limit
//**************************************************************************************
void Setup_Taylor_Series_Coefficients_Triangle(Finite_Diff_Parameters &Finite_Diff_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed) {
    // In the case where the set of angular moments up to second order of interest lies 
    // within the triangle (P1 P2 P3) but is sufficiently close to either one of the 
    // edges of the triangle (P1 P2 P3), or one of the vertices of the triangle (P1 P2 P3), 
    // we perform a Taylor series expansion along the median line passing through such a 
    // point and the barycenter of the triangle
    // On the other hand, if the point of interest lies along one of the edges of the 
    // triangle (P1 P2 P3), and is sufficiently close to one of the vertices, we then 
    // apply a taylor series expansion along such a edge to obtain an approximation of
    // the maximum entropy solutions at that point
    M2_State_Param M2_State_Temp(M2_State);
    record_N3 *rec_N3_local_diff;
    long double h;
    int order, prec, n, nmax;
    long double *c, *x;
    long double d_N3_111, d_N3_122, d_N3_123;
    int Domain_Type;
    order = 4;
    prec = 2;
    n = order + prec;
    nmax = n;
    c = new long double[n];
    x = new long double[n];
    rec_N3_local_diff = new record_N3[n];
    h = finite_diff_h_gam;
    
    // Setup parameters for performing finite difference
    M2_State_Temp.h_gam = h;
    Setup_finite_diff_Parameters_Triangle_Taylor_Series(M2_State_Temp, Finite_Diff_N3);
    
    Finite_Diff_N3.gam1_0 = M2_State_Temp.gam1_knot;
    Finite_Diff_N3.gam2_0 = M2_State_Temp.gam2_knot;
    
    switch(M2_State_Temp.finite_diff_type) {
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            exit(0);
            break;
    }
    
    M2_State_Temp.norm_f_knot = 1.0;
    
    M2_State_Temp.flag_finite_diff_Sphere = false;
    M2_State_Temp.flag_finite_diff_Triangle = true;
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M2_State_Temp.id_finite_diff_norm_f = id_diff;
        M2_State_Temp.x_finite_diff_norm_f = x[id_diff];
        
        Check_Domain_Type(Domain_Type, Var_index, num_points, M2_State_Temp);
        if ((Domain_Type != BOUNDARY_GAM1_GAM2_EQ_0) || (Domain_Type != BOUNDARY_GAM1_GAM3_EQ_0) || (Domain_Type != BOUNDARY_GAM2_GAM3_EQ_0)) {
            if (M2_State_Temp.Domain_Type != Domain_Type) {
                M2_State_Temp.Reset_Num_Vars(Domain_Type);
            }
        }
        
        NLOPT_Optim_Algo(rec_N3_local_diff, id_diff, M2_State_Temp, Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
        
        if (M2_State_Temp.finite_diff_domain_Triangle == MEDIAN_TRIANGLE_FINITE_DIFFERENCE) {
            switch(M2_State_Temp.finite_diff_type) {
                case FORWARD_FINITE_DIFFERENCE:
                    if (id_diff == 0) {
                        Finite_Diff_N3.N3_111_x0 = rec_N3_local_diff[id_diff].N3_111;
                        Finite_Diff_N3.N3_122_x0 = rec_N3_local_diff[id_diff].N3_122;
                        Finite_Diff_N3.N3_123_x0 = rec_N3_local_diff[id_diff].N3_123;
                    }
                    break;
                case BACKWARD_FINITE_DIFFERENCE:
                    if (id_diff == n - 1) {
                        Finite_Diff_N3.N3_111_x0 = rec_N3_local_diff[id_diff].N3_111;
                        Finite_Diff_N3.N3_122_x0 = rec_N3_local_diff[id_diff].N3_123;
                        Finite_Diff_N3.N3_123_x0 = rec_N3_local_diff[id_diff].N3_123;
                    }
                    break;
                default:
                    cout << "Finite difference type not specified" << endl;
                    exit(0);
                    break;
            }
        }
    }
    
    d_N3_111 = 0.0;
    d_N3_122 = 0.0;
    d_N3_123 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_111;
        d_N3_122 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_122;
        d_N3_123 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_123;
    }
    
    Finite_Diff_N3.d4N3_111_dgam = d_N3_111;
    Finite_Diff_N3.d4N3_122_dgam = d_N3_122;
    Finite_Diff_N3.d4N3_123_dgam = d_N3_123;
    
    // backward finite difference approximation to third derivative
    order = 3;
    n = order + prec;
    switch(M2_State_Temp.finite_diff_type) {
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            exit(0);
            break;
    }
    
    d_N3_111 = 0.0;
    d_N3_122 = 0.0;
    d_N3_123 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_111;
        d_N3_122 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_122;
        d_N3_123 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_123;
    }
    
    Finite_Diff_N3.d3N3_111_dgam = d_N3_111;
    Finite_Diff_N3.d3N3_122_dgam = d_N3_122;
    Finite_Diff_N3.d3N3_123_dgam = d_N3_123;
    
    // backward finite difference approximation to second derivative
    order = 2;
    n = order + prec;
    switch(M2_State_Temp.finite_diff_type) {
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            exit(0);
            break;
    }
    
    d_N3_111 = 0.0;
    d_N3_122 = 0.0;
    d_N3_123 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_111;
        d_N3_122 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_122;
        d_N3_123 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_123;
    }
    
    Finite_Diff_N3.d2N3_111_dgam = d_N3_111;
    Finite_Diff_N3.d2N3_122_dgam = d_N3_122;
    Finite_Diff_N3.d2N3_123_dgam = d_N3_123;
    
    // backward finite difference approximation to first derivative
    order = 1;
    n = order + prec;
    switch(M2_State_Temp.finite_diff_type) {
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            exit(0);
            break;
    }
    
    d_N3_111 = 0.0;
    d_N3_122 = 0.0;
    d_N3_123 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_111;
        d_N3_122 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_122;
        d_N3_123 += c[n - i - 1] * rec_N3_local_diff[nmax - i - 1].N3_123;
    }
    
    Finite_Diff_N3.dN3_111_dgam = d_N3_111;
    Finite_Diff_N3.dN3_122_dgam = d_N3_122;
    Finite_Diff_N3.dN3_123_dgam = d_N3_123;
    
    if (M2_State_Temp.display) {
        if (M2_State_Temp.id_proc == M2_State.proc_display) {
            for (int id_diff = 0; id_diff < nmax; id_diff++) {
                cout << "id_diff = " << id_diff << "   " << "x = " << x[id_diff] << "   " << "N3_111 = " << rec_N3_local_diff[id_diff].N3_111 << "   " << "N3_122 = " << rec_N3_local_diff[id_diff].N3_122 << "   " << "N3_123 = " << rec_N3_local_diff[id_diff].N3_123 << endl;
            }
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] rec_N3_local_diff;
}

//**************************************************************************************
// This routines performs finite differencing along the radius of the unit ball defined
// by the realizable space for the first-order angular moments in order to compute 
// derivatives of the third-order closing fluxes with respect to the norm of the first-
// order angular moments
//**************************************************************************************
void Calculate_Higher_Order_Moments_Derivatives_Boundary_Norm_f(record_N3 &rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed) {
    // Perform finite difference approximation to obtain derivative
    M2_State_Param M2_State_Temp(M2_State);
    record_N3 *rec_N3_local_diff;
    long double h;
    int order, prec, n;
    long double *c, *x;
    long double d_N3_111_dnorm_f, d_N3_122_dnorm_f, d_N3_123_dnorm_f;
    long double d2_N3_111_dnorm_f_dN1_1, d2_N3_122_dnorm_f_dN1_1;
    int Domain_Type;
    order = 1;
    prec = 2;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    rec_N3_local_diff = new record_N3[n];
    h = finite_diff_h_N1;
    
    // Compute values for performing finite difference
    M2_State_Temp.finite_diff_domain_Sphere = RADIUS_SPHERE_FINITE_DIFFERENCE;
    M2_State_Temp.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
    M2_State_Temp.h_norm_f = h;
    
    // backward finite difference approximation to first derivative
    differ_backward ( M2_State_Temp.h_norm_f, order, prec, c, x );
    
    M2_State_Temp.norm_f_knot = 1.0;
    
    M2_State_Temp.flag_finite_diff_Sphere = true;
    M2_State_Temp.flag_finite_diff_Triangle = false;
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M2_State_Temp.id_finite_diff_norm_f = id_diff;
        M2_State_Temp.x_finite_diff_norm_f = x[id_diff];
        
        Check_Domain_Type(Domain_Type, Var_index, num_points, M2_State_Temp);
        if ((Domain_Type != BOUNDARY_GAM1_GAM2_EQ_0) || (Domain_Type != BOUNDARY_GAM1_GAM3_EQ_0) || (Domain_Type != BOUNDARY_GAM2_GAM3_EQ_0)) {
            if (M2_State_Temp.Domain_Type != Domain_Type) {
                M2_State_Temp.Reset_Num_Vars(Domain_Type);
            }
        }
        
        NLOPT_Optim_Algo(rec_N3_local_diff, id_diff, M2_State_Temp, Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
    }
    
    // We must have dN3_dgam = 0.0 in the free streaming limit ???
    
    d_N3_111_dnorm_f = 0.0;
    d_N3_122_dnorm_f = 0.0;
    
    d2_N3_111_dnorm_f_dN1_1 = 0.0;
    d2_N3_122_dnorm_f_dN1_1 = 0.0;
    
    for (int i = 0; i < n; i++) {
        d_N3_111_dnorm_f += c[i] * rec_N3_local_diff[i].N3_111;
        d_N3_122_dnorm_f += c[i] * rec_N3_local_diff[i].N3_122;
        
        d2_N3_111_dnorm_f_dN1_1 += c[i] * rec_N3_local_diff[i].dN3_111_dN1_1;
        d2_N3_122_dnorm_f_dN1_1 += c[i] * rec_N3_local_diff[i].dN3_122_dN1_1;
    }
    
    rec_N3.dN3_111_dnorm_f = d_N3_111_dnorm_f;
    rec_N3.dN3_122_dnorm_f = d_N3_122_dnorm_f;
    
    rec_N3.d2_N3_111_dnorm_f_dN1_1 = d2_N3_111_dnorm_f_dN1_1;
    rec_N3.d2_N3_122_dnorm_f_dN1_1 = d2_N3_122_dnorm_f_dN1_1;
    
    // In this case: N3_111 = (N1_1)^3
    //               N3_122 = N1_1*(N1_2)^2
    rec_N3.dN3_111_dN1_1 = 3.0*pow(rec_N3.N1_1, 2);
    rec_N3.dN3_122_dN1_1 = pow(rec_N3.N1_2, 2);
            
    if (M2_State.display) {
        if (M2_State.id_proc == M2_State.proc_display) {
            cout << "d_N3_111_dnorm_f = " << rec_N3.dN3_111_dnorm_f << "  " << "d_N3_122_dnorm_f = " << rec_N3.dN3_122_dnorm_f << endl;
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] rec_N3_local_diff;
}

//**************************************************************************************
// This routines performs finite differencing along the line passing through any given
// point on the edges of the triangle (P1 P2 P3) and the barycenter of the triangle in 
// order to compute derivatives of the third-order closing fluxes with respect to the 
// eigenvalues of the covariance matrix
//**************************************************************************************
void finite_diff_triangle_median(M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere, const long double &h, const int &order, const int &prec, record_Derivatives_Triangle &diff_median) {
    int Domain_Type;
    record_N3 *rec_N3_local_diff;
    int n;
    long double *c, *x;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    rec_N3_local_diff = new record_N3[n];
    
    // Compute values for performing finite difference along the line passing
    // through the point on one of the edges of the triangle and the barycenter
    // of the triangle
    Setup_finite_diff_Triangle_Median(M2_State, h);
    M2_State_Param M2_State_Temp(M2_State);
    
    M2_State_Temp.finite_diff_domain_Triangle = MEDIAN_TRIANGLE_FINITE_DIFFERENCE;
    
    switch(M2_State_Temp.finite_diff_type) {
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            exit(0);
            break;
    }
    
    M2_State_Temp.flag_finite_diff_Triangle = true;
    
//     for (int id_diff = 0; id_diff < n; id_diff++) {
//         cout << "BBBBBBBBBBBBBBBBBBBBBB" << "id_diff = " << id_diff << "  " << "x = " << x[id_diff] << "  " << "h = " << M2_State_Temp.h_gam << "  " << "h orig = " << h << endl;
//     }
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M2_State_Temp.id_finite_diff_gams = id_diff;
        M2_State_Temp.x_finite_diff_gam1 = x[id_diff]*M2_State_Temp.cos_Median_Line;
        M2_State_Temp.x_finite_diff_gam2 = x[id_diff]*M2_State_Temp.sin_Median_Line;
        
        // cout << "BBBBBBBBBBBBBBBBBBBBBB" << "id_diff = " << id_diff << "  " << "x = " << x[id_diff] << "  " << "gam1_knot = " << M2_State_Temp.gam1_knot << "  " << "gam2_knot = " << M2_State_Temp.gam2_knot << "  " << "x_gam1 = " << M2_State_Temp.x_finite_diff_gam1 << "  " << "x_gam2 = " << M2_State_Temp.x_finite_diff_gam2 << endl;
        
        Check_Domain_Type(Domain_Type, Var_index, num_points, M2_State_Temp);
        if ((Domain_Type != BOUNDARY_GAM1_GAM2_EQ_0) && (Domain_Type != BOUNDARY_GAM1_GAM3_EQ_0) && (Domain_Type != BOUNDARY_GAM2_GAM3_EQ_0)) {
            if (M2_State_Temp.Domain_Type != Domain_Type) {
                M2_State_Temp.Reset_Num_Vars(Domain_Type);
            }
        }
        
        NLOPT_Optim_Algo(rec_N3_local_diff, id_diff, M2_State_Temp, Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
    }
    
    diff_median.d_N3_111 = 0.0;
    diff_median.d_N3_122 = 0.0;
    diff_median.d_N3_123 = 0.0;
    
    diff_median.d2_N3_111_dgam_dN1_1 = 0.0;
    diff_median.d2_N3_122_dgam_dN1_1 = 0.0;
    
    diff_median.d2_N3_111_dnorm_f_dgam = 0.0;
    diff_median.d2_N3_122_dnorm_f_dgam = 0.0;
    diff_median.d3_N3_111_dnorm_f_dgam_dN1_1 = 0.0;
    diff_median.d3_N3_122_dnorm_f_dgam_dN1_1 = 0.0;
    
    for (int i = 0; i < n; i++) {
//         cout << "i = " << i << "  " << "h = " << M2_State_Temp.h_gam << "  " << "c = " << c[i] << "  " << "N3_122 = " << rec_N3_local_diff[i].N3_122 << endl;
        diff_median.d_N3_111 += c[i] * rec_N3_local_diff[i].N3_111;
        diff_median.d_N3_122 += c[i] * rec_N3_local_diff[i].N3_122;
        diff_median.d_N3_123 += c[i] * rec_N3_local_diff[i].N3_123;
        
        diff_median.d2_N3_111_dgam_dN1_1 += c[i] * rec_N3_local_diff[i].dN3_111_dN1_1;
        diff_median.d2_N3_122_dgam_dN1_1 += c[i] * rec_N3_local_diff[i].dN3_122_dN1_1;
        
        if (flag_finite_diff_Sphere) {
            diff_median.d2_N3_111_dnorm_f_dgam += c[i] * rec_N3_local_diff[i].dN3_111_dnorm_f;
            diff_median.d2_N3_122_dnorm_f_dgam += c[i] * rec_N3_local_diff[i].dN3_122_dnorm_f;
            diff_median.d3_N3_111_dnorm_f_dgam_dN1_1 += c[i] * rec_N3_local_diff[i].d2_N3_111_dnorm_f_dN1_1;
            diff_median.d3_N3_122_dnorm_f_dgam_dN1_1 += c[i] * rec_N3_local_diff[i].d2_N3_122_dnorm_f_dN1_1;
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] rec_N3_local_diff;
}

//**************************************************************************************
// This routines performs finite differencing along any given edge of the triangle 
// (P1 P2 P3) in order to compute derivatives of the third-order closing fluxes with 
// respect to the eigenvalues of the covariance matrix
//**************************************************************************************
void finite_diff_triangle_edges(M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere, const long double &h, const int &order, const int &prec, record_Derivatives_Triangle &diff_edge) {
    int Domain_Type;
    M2_State_Param M2_State_Temp(M2_State);
    record_N3 *rec_N3_local_diff;
    int n;
    long double *c, *x;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    rec_N3_local_diff = new record_N3[n];
    
    // Now compute values for performing finite difference the edge of the triangle
    // that is of interest
    M2_State_Temp.finite_diff_domain_Triangle = EDGE_TRIANGLE_FINITE_DIFFERENCE;
    M2_State_Temp.gam1_knot = M2_State.gamma_1;
    M2_State_Temp.gam2_knot = M2_State.gamma_2;
    
    switch (M2_State_Temp.face_triangle) {
        case FACE_A:
            M2_State_Temp.h_gam1 = 0.0;
            M2_State_Temp.h_gam2 = h;
            
            if (M2_State_Temp.gam2_knot - ceil(n/2.0)*M2_State_Temp.h_gam2 < 0.0) {
                M2_State_Temp.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
            } else if (M2_State_Temp.gam2_knot + ceil(n/2.0)*M2_State_Temp.h_gam2 > 1.0) {
                M2_State_Temp.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
            } else {
                M2_State_Temp.finite_diff_type = CENTRAL_FINITE_DIFFERENCE;
            }
            break;
        case FACE_B:
            M2_State_Temp.h_gam1 = h;
            M2_State_Temp.h_gam2 = 0.0;
            
            if (M2_State_Temp.gam1_knot - ceil(n/2.0)*M2_State_Temp.h_gam1 < 0.0) {
                M2_State_Temp.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
            } else if (M2_State_Temp.gam1_knot + ceil(n/2.0)*M2_State_Temp.h_gam1 > 1.0) {
                M2_State_Temp.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
            } else {
                M2_State_Temp.finite_diff_type = CENTRAL_FINITE_DIFFERENCE;
            }
            break;
        case FACE_C:
            // h cos theta or sin
            M2_State_Temp.h_gam1 = h; //*sqrt(2.0)/2.0;
            M2_State_Temp.h_gam2 = -h; //*sqrt(2.0)/2.0;
            
            if (M2_State_Temp.gam1_knot - ceil(n/2.0)*M2_State_Temp.h_gam1 < 0.0) {
                M2_State_Temp.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
            } else if (M2_State_Temp.gam1_knot + ceil(n/2.0)*M2_State_Temp.h_gam1 > 1.0) {
                M2_State_Temp.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
            } else {
                M2_State_Temp.finite_diff_type = CENTRAL_FINITE_DIFFERENCE;
            }
            break;
        default:
            cout << "Face Type for triangle finite difference not specified" << endl;
            exit(0);
            break;
    };
    
    M2_State_Temp.h_gam = sqrt(pow(M2_State_Temp.h_gam1, 2) + pow(M2_State_Temp.h_gam2, 2));
    
    M2_State_Temp.flag_finite_diff_Triangle = true;
    
    switch(M2_State_Temp.finite_diff_type) {
        case CENTRAL_FINITE_DIFFERENCE:
            // central finite difference approximation to first derivative
            differ_central ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            break;
    }
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M2_State_Temp.id_finite_diff_gams = id_diff;
        
        M2_State_Temp.x_finite_diff_gam1 = x[id_diff]*M2_State_Temp.h_gam1/h;
        M2_State_Temp.x_finite_diff_gam2 = x[id_diff]*M2_State_Temp.h_gam2/h;
        
        Check_Domain_Type(Domain_Type, Var_index, num_points, M2_State_Temp);
        if ((Domain_Type != BOUNDARY_GAM1_GAM2_EQ_0) || (Domain_Type != BOUNDARY_GAM1_GAM3_EQ_0) || (Domain_Type != BOUNDARY_GAM2_GAM3_EQ_0)) {
            if (M2_State_Temp.Domain_Type != Domain_Type) {
                M2_State_Temp.Reset_Num_Vars(Domain_Type);
            }
        }
        
        NLOPT_Optim_Algo(rec_N3_local_diff, id_diff, M2_State_Temp, Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
    }
    
    diff_edge.d_N3_111 = 0.0;
    diff_edge.d_N3_122 = 0.0;
    diff_edge.d_N3_123 = 0.0;
    
    diff_edge.d2_N3_111_dgam_dN1_1 = 0.0;
    diff_edge.d2_N3_122_dgam_dN1_1 = 0.0;
    
    diff_edge.d2_N3_111_dnorm_f_dgam = 0.0;
    diff_edge.d2_N3_122_dnorm_f_dgam = 0.0;
    diff_edge.d3_N3_111_dnorm_f_dgam_dN1_1 = 0.0;
    diff_edge.d3_N3_122_dnorm_f_dgam_dN1_1 = 0.0;
    
    for (int i = 0; i < n; i++) {
//         cout << "i = " << i << "  " << "x = " << x[i] << endl;
        diff_edge.d_N3_111 += c[i] * rec_N3_local_diff[i].N3_111;
        diff_edge.d_N3_122 += c[i] * rec_N3_local_diff[i].N3_122;
        diff_edge.d_N3_123 += c[i] * rec_N3_local_diff[i].N3_123;
        
        diff_edge.d2_N3_111_dgam_dN1_1 += c[i] * rec_N3_local_diff[i].dN3_111_dN1_1;
        diff_edge.d2_N3_122_dgam_dN1_1 += c[i] * rec_N3_local_diff[i].dN3_122_dN1_1;
        
        if (flag_finite_diff_Sphere) {
            diff_edge.d2_N3_111_dnorm_f_dgam += c[i] * rec_N3_local_diff[i].dN3_111_dnorm_f;
            diff_edge.d2_N3_122_dnorm_f_dgam += c[i] * rec_N3_local_diff[i].dN3_122_dnorm_f;
            diff_edge.d3_N3_111_dnorm_f_dgam_dN1_1 += c[i] * rec_N3_local_diff[i].d2_N3_111_dnorm_f_dN1_1;
            diff_edge.d3_N3_122_dnorm_f_dgam_dN1_1 += c[i] * rec_N3_local_diff[i].d2_N3_122_dnorm_f_dN1_1;
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] rec_N3_local_diff;
}

//*****************************************************************************************************
// This routines performs finite differencing both along the line passing through any given point on the 
// edges of the triangle (P1 P2 P3) and the barycenter of the triangle, as well as finite differencing 
// along the correspoding edge of the triangle so as to compute derivatives of the third-order closing 
// fluxes with respect to the eigenvalues of the third-order closing fluxes
//*****************************************************************************************************
void Calculate_Higher_Order_Moments_Derivatives_Boundary_gam1_gam2_gam3(record_N3 &rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere) {
    // Perform finite difference approximation to obtain derivative
    M2_State_Param M2_State_Temp(M2_State);
    long double u1, u2, v1, v2;
    long double h;
    int order, prec;
    record_Derivatives_Triangle diff_median, diff_edge;
    
    order = 1;
    prec = 2;
    h = finite_diff_h_gam;
    
    // Setup the reference values for gam1 and gam2 used in the finite differencing procedure
    M2_State_Temp.gam1_knot = M2_State.gamma_1;
    M2_State_Temp.gam2_knot = M2_State.gamma_2;
    
    // Check what boundary the current point of interest belongs
    switch (M2_State_Temp.Domain_Type) {
        case BOUNDARY_GAM1:
        case BOUNDARY_GAM1_GAM2_EQ_0:
        case BOUNDARY_GAM1_GAM3_EQ_0:
            M2_State_Temp.face_triangle = FACE_A;
            break;
        case BOUNDARY_GAM2:
        case BOUNDARY_GAM2_GAM3_EQ_0:
            M2_State_Temp.face_triangle = FACE_B;
            break;
        case BOUNDARY_GAM3:
            M2_State_Temp.face_triangle = FACE_C;
            break;
        default:
            cout << "Domain Type for triangle finite difference not specified" << endl;
            exit(0);
            break;
    };
    
    if (M2_State_Temp.display) {
        cout << "Performing finite differencing along the edges" << endl;
    }
    
    // Now perform finite differencing along the edge of interest
    finite_diff_triangle_edges(M2_State_Temp, Var_index, num_points, id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, flag_finite_diff_Sphere, h, order, prec, diff_edge);
    
    if (M2_State_Temp.display) {
        cout << "Performing finite differencing along the median" << endl;
    }
    
    // Now perform finite differencing along the line passing through the the point of the edge of interest and
    // and barycenter of the triangle
    finite_diff_triangle_median(M2_State_Temp, Var_index, num_points, id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, flag_finite_diff_Sphere, h, order, prec, diff_median);
    
    // If we know the directional derivative along two different directions
    // We can then solve the resulting system of two equations for the two
    // unknowns that are the derivatives along the Cartesian axes
    switch (M2_State_Temp.Domain_Type) {
        // Direction vector along the line defined by a gam1 + b gam2 + c = 0 (edge )
        // is u1 = {b, -a}
        case BOUNDARY_GAM1:
        case BOUNDARY_GAM1_GAM2_EQ_0:
        case BOUNDARY_GAM1_GAM3_EQ_0:
            // For FACE_A (a = 1, b = 0)
            u1 = 0.0;
            u2 = 1.0;
            break;
        case BOUNDARY_GAM2:
        case BOUNDARY_GAM2_GAM3_EQ_0:
            // For FACE_B (a = 0, b = 1)
            u1 = 1.0;
            u2 = 0.0;
            break;
        case BOUNDARY_GAM3:
            // For FACE_C (a = 1, b = 1)
            u1 = 1.0/sqrt(2.0);
            u2 = -1.0/sqrt(2.0);
            break;
        default:
            cout << "Domain Type for triangle finite difference not specified" << endl;
            exit(0);
            break;
    };
    // Direction vector along the line defined by a gam1 + b gam2 + c = 0 (median)
    // is u2 = {b, -a}
    // For median we have equation gam2 = a_median * gam1 + c ==> (a = -a_median, b = 1)
    Median_Line_Unit_Vector(v1, v2, M2_State_Temp);
    
    // u is the unitt vector along the edge of interest and v is the unit
    // vector along the line passing through the barycenter of the triangle
    // and the point on the edge of the triangle
    
    // We then have the following system of equations
    // Directional derivative along u
    // u1 * dN3ijk_dgam1 + u2 * dN3ijk_dgam2 = d_N3_ijk_edge
    // Directional derivative along v
    // v1 * dN3ijk_dgam1 + v2 * dN3ijk_dgam2 = d_N3_ijk_median
    // The solution is
    // dN3ijk_dgam1 = (v2 * d_N3_ijk_edge - u2 * d_N3_ijk_median)/(u1 * v2 - u2 * v1);
    // dN3ijk_dgam2 = (u1 * d_N3_ijk_median - v1 * d_N3_ijk_edge)/(u1 * v2 - u2 * v1);
    long double determinant;
    determinant = u1*v2 - u2*v1;
    
    rec_N3.dN3_111_dgam1 = (v2*diff_edge.d_N3_111 - u2*diff_median.d_N3_111)/determinant;
    rec_N3.dN3_122_dgam1 = (v2*diff_edge.d_N3_122 - u2*diff_median.d_N3_122)/determinant;
    rec_N3.dN3_123_dgam1 = (v2*diff_edge.d_N3_123 - u2*diff_median.d_N3_123)/determinant;
    
    rec_N3.d2_N3_111_dgam1_dN1_1 = (v2*diff_edge.d2_N3_111_dgam_dN1_1 - u2*diff_median.d2_N3_111_dgam_dN1_1)/determinant;
    rec_N3.d2_N3_122_dgam1_dN1_1 = (v2*diff_edge.d2_N3_122_dgam_dN1_1 - u2*diff_median.d2_N3_122_dgam_dN1_1)/determinant;
    
    rec_N3.d2_N3_111_dnorm_f_dgam1 = (v2*diff_edge.d2_N3_111_dnorm_f_dgam - u2*diff_median.d2_N3_111_dnorm_f_dgam)/determinant;
    rec_N3.d2_N3_122_dnorm_f_dgam1 = (v2*diff_edge.d2_N3_122_dnorm_f_dgam - u2*diff_median.d2_N3_122_dnorm_f_dgam)/determinant;
    
    rec_N3.d3_N3_111_dnorm_f_dgam1_dN1_1 = (v2*diff_edge.d3_N3_111_dnorm_f_dgam_dN1_1 - u2*diff_median.d3_N3_111_dnorm_f_dgam_dN1_1)/determinant;
    rec_N3.d3_N3_122_dnorm_f_dgam1_dN1_1 = (v2*diff_edge.d3_N3_122_dnorm_f_dgam_dN1_1 - u2*diff_median.d3_N3_122_dnorm_f_dgam_dN1_1)/determinant;
    
    rec_N3.dN3_111_dgam2 = (u1*diff_median.d_N3_111 - v1*diff_edge.d_N3_111)/determinant;
    rec_N3.dN3_122_dgam2 = (u1*diff_median.d_N3_122 - v1*diff_edge.d_N3_122)/determinant;
    rec_N3.dN3_123_dgam2 = (u1*diff_median.d_N3_123 - v1*diff_edge.d_N3_123)/determinant;
    
//     rec_N3.d2_N3_111_dgam2_dN1_1 = (u1*diff_median.d2_N3_111_dgam_dN1_1 - v1*diff_edge.d2_N3_111_dgam_dN1_1)/determinant;
    rec_N3.d2_N3_122_dgam2_dN1_1 = (u1*diff_median.d2_N3_122_dgam_dN1_1 - v1*diff_edge.d2_N3_122_dgam_dN1_1)/determinant;
    
//     rec_N3.d2_N3_111_dnorm_f_dgam2 = (u1*diff_median.d2_N3_111_dnorm_f_dgam - v1*diff_edge.d2_N3_111_dnorm_f_dgam)/determinant;
    rec_N3.d2_N3_122_dnorm_f_dgam2 = (u1*diff_median.d2_N3_122_dnorm_f_dgam - v1*diff_edge.d2_N3_122_dnorm_f_dgam)/determinant;
    
//     rec_N3.d3_N3_111_dnorm_f_dgam2_dN1_1 = (u1*diff_median.d3_N3_111_dnorm_f_dgam_dN1_1 - v1*diff_edge.d3_N3_111_dnorm_f_dgam_dN1_1)/determinant;
    rec_N3.d3_N3_122_dnorm_f_dgam2_dN1_1 = (u1*diff_median.d3_N3_122_dnorm_f_dgam_dN1_1 - v1*diff_edge.d3_N3_122_dnorm_f_dgam_dN1_1)/determinant;
            
    if (M2_State_Temp.display) {
        // cout << "median = " << diff_median.d_N3_111 << "  " << "edge = " << diff_edge.d_N3_111 << "  " << "median = " << diff_median.d_N3_122 << "  " << "edge = " << diff_edge.d_N3_122 << "  " << "dN3_111_dgam1 = " << rec_N3.dN3_111_dgam1 << "  " << "dN3_111_dgam2 = " << rec_N3.dN3_111_dgam2 << "  " << "dN3_122_dgam1 = " << rec_N3.dN3_122_dgam1 << "  " << "dN3_122_dgam2 = " << rec_N3.dN3_122_dgam2 << "  " << "dN3_123_dgam1 = " << rec_N3.dN3_123_dgam1 << "  " << "dN3_123_dgam2 = " << rec_N3.dN3_123_dgam2 << endl;
    }
}

void Setup_finite_diff_Triangle_Median(M2_State_Param &M2_State, const long double &h) {
    long double d_to_min, d_to_max;
    long double u1 ,u2;
    
    // Compute parameters for the equation of the median line
    M2_State.a_median = a_median(M2_State);
    M2_State.b_median = b_median(M2_State);
    
    // Find intersections of median line with edges of the triangle
    Triangle_Median_Edge_Intersection_Points(M2_State);
    
    // Compare the distances between the point of interest and the min and max values
    // to see which one it is closer to and consequently set the knot value 
    d_to_min = pow(M2_State.gamma_1 - M2_State.gam1_bound_min, 2) +  pow(M2_State.gamma_2 - M2_State.gam2_bound_min, 2);
    d_to_max = pow(M2_State.gamma_1 - M2_State.gam1_bound_max, 2) +  pow(M2_State.gamma_2 - M2_State.gam2_bound_max, 2);
    
    if (d_to_min < d_to_max) {
        // Then  we are close to the "min" value
        M2_State.gam1_knot = M2_State.gam1_bound_min;
        M2_State.gam2_knot = M2_State.gam2_bound_min;
    } else {
        // Then we are closer to the "max" value
        M2_State.gam1_knot = M2_State.gam1_bound_max;
        M2_State.gam2_knot = M2_State.gam2_bound_max;
    }
    M2_State.gam3_knot = 1.0 - M2_State.gam1_knot - M2_State.gam2_knot;
    
    Median_Line_Unit_Vector(u1, u2, M2_State);
    
    M2_State.cos_Median_Line = u1;
    M2_State.sin_Median_Line = u2;
    
    M2_State.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
    long double tan_Median_Line, tan_Median_Line_gam1_gam3;
    long double h_gam3;
    long double u1_gam1_gam3, u2_gam1_gam3;
    // If the equation of the median line is gam2 = a * gam1  + b
    // and we have gam3 = 1 - gam1 - gam2
    // We then also have: gam3 = -(1 + a) * gam1  + (1 - b)
    Median_Line_Unit_Vector_Plane_Gam1_Gam3(u1_gam1_gam3, u2_gam1_gam3, M2_State);
    
    tan_Median_Line = M2_State.sin_Median_Line/M2_State.cos_Median_Line;
    tan_Median_Line_gam1_gam3 = u2_gam1_gam3/u1_gam1_gam3;
    
    
    if (M2_State.gam1_knot < 1.0e-8) {
        // Then we are along FACE A ==> gam1 = 0
        M2_State.h_gam1 = h;
        M2_State.h_gam2 = fabs(M2_State.h_gam1*tan_Median_Line);
    } else if (M2_State.gam2_knot < 1.0e-8) {
        // Then we are along FACE B ==> gam2 = 0
        M2_State.h_gam2 = h;
        M2_State.h_gam1 = fabs(M2_State.h_gam2/tan_Median_Line);
    } else if (1.0 - M2_State.gam1_knot - M2_State.gam2_knot < 1.0e-8) {
        h_gam3 = h;
        M2_State.h_gam1 = fabs(h_gam3/tan_Median_Line_gam1_gam3);
        M2_State.h_gam2 = fabs(M2_State.h_gam1*tan_Median_Line);
    } else {
        cout << "Issues with determining h_gam1 and h_gam2 for median line !!!!!!!!!!!" << endl;
        exit(0);
    }
    M2_State.h_gam = sqrt(pow(M2_State.h_gam1, 2) + pow(M2_State.h_gam2, 2));
    
    // cout << "h_gam = " << M2_State.h_gam << "  " << "h_gam1 = " << M2_State.h_gam1 << "  " << "h_gam2 = " << M2_State.h_gam2 << "  " << "cos = " << M2_State.cos_Median_Line << "  " << "sin = " << M2_State.sin_Median_Line << "  " << "tan = " << tan_Median_Line << endl;
}

//*****************************************************************************************************
// This routines performs finite differencing along the line passing through any of the vertices of the
// triangle (P1 P2 P3) and the barycenter of the triangle
//*****************************************************************************************************
void Calculate_Higher_Order_Moments_Derivatives_Vertex_Median(record_Derivatives_Triangle &rec_dN3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere, const Finite_Diff_Data &rec_Finite_Diff) {
    int Domain_Type;
    record_N3 *rec_N3_local_diff;
    int n, nmax;
    long double *c, *x;
    int order, prec;
    order = rec_Finite_Diff.order;
    prec = rec_Finite_Diff.prec;
    n = order + prec;
    nmax = n;
    c = new long double[n];
    x = new long double[n];
    rec_N3_local_diff = new record_N3[n];
    
    if (order != 2) {
        cout << "Issue with finite difference along median line passing through vertex !!" << endl;
        exit(0);
    }
    
    // Compute values for performing finite difference along the line passing
    // through the point on one of the edges of the triangle and the barycenter
    // of the triangle
    Setup_finite_diff_Triangle_Median(M2_State, rec_Finite_Diff.h);
    M2_State_Param M2_State_Temp(M2_State);
    M2_State_Temp.finite_diff_domain_Triangle = MEDIAN_TRIANGLE_FINITE_DIFFERENCE;
    M2_State_Temp.flag_finite_diff_Triangle = true;
    
    switch(M2_State_Temp.finite_diff_type) {
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            exit(0);
            break;
    }
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M2_State_Temp.id_finite_diff_gams = id_diff;
        M2_State_Temp.x_finite_diff_gam1 = x[id_diff]*M2_State_Temp.cos_Median_Line;
        M2_State_Temp.x_finite_diff_gam2 = x[id_diff]*M2_State_Temp.sin_Median_Line;
        
        // cout << "AAAAAAAAAAAAAAAAAAAAA" << "id_diff = " << id_diff << "  " << "gam1_knot = " << M2_State_Temp.gam1_knot << "  " << "gam2_knot = " << M2_State_Temp.gam2_knot << "  " << "x_gam1 = " << M2_State_Temp.x_finite_diff_gam1 << "  " << "x_gam2 = " << M2_State_Temp.x_finite_diff_gam2 << endl;
        
        Check_Domain_Type(Domain_Type, Var_index, num_points, M2_State_Temp);
        if ((Domain_Type != BOUNDARY_GAM1_GAM2_EQ_0) && (Domain_Type != BOUNDARY_GAM1_GAM3_EQ_0) && (Domain_Type != BOUNDARY_GAM2_GAM3_EQ_0)) {
            if (M2_State_Temp.Domain_Type != Domain_Type) {
                M2_State_Temp.Reset_Num_Vars(Domain_Type);
            }
        }
        
        NLOPT_Optim_Algo(rec_N3_local_diff, id_diff, M2_State_Temp, Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
    }
    
    // First perform finite difference for second-order derivatives
    rec_dN3.d2_N3_111 = 0.0;
    rec_dN3.d2_N3_122 = 0.0;
    rec_dN3.d2_N3_123 = 0.0;
    
    rec_dN3.d3_N3_111_dgam_dN1_1 = 0.0;
    rec_dN3.d3_N3_122_dgam_dN1_1 = 0.0;
    rec_dN3.d3_N3_111_dnorm_f_dgam = 0.0;
    rec_dN3.d3_N3_122_dnorm_f_dgam = 0.0;
    rec_dN3.d4_N3_111_dnorm_f_dgam_dN1_1 = 0.0;
    rec_dN3.d4_N3_122_dnorm_f_dgam_dN1_1 = 0.0;
    
    int id_c, id_vals;
    for (int i = 0; i < n; i++) {
        id_c = n - i - 1;
        id_vals = nmax - i - 1;
        
        rec_dN3.d_N3_111 += c[id_c] * rec_N3_local_diff[id_vals].N3_111;
        rec_dN3.d_N3_122 += c[id_c] * rec_N3_local_diff[id_vals].N3_122;
        rec_dN3.d_N3_123 += c[id_c] * rec_N3_local_diff[id_vals].N3_123;
        
        rec_dN3.d2_N3_111_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].dN3_111_dN1_1;
        rec_dN3.d2_N3_122_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].dN3_122_dN1_1;
        
        if (flag_finite_diff_Sphere) {
            rec_dN3.d2_N3_111_dnorm_f_dgam += c[id_c] * rec_N3_local_diff[id_vals].dN3_111_dnorm_f;
            rec_dN3.d2_N3_122_dnorm_f_dgam += c[id_c] * rec_N3_local_diff[id_vals].dN3_122_dnorm_f;
            rec_dN3.d3_N3_111_dnorm_f_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].d2_N3_111_dnorm_f_dN1_1;
            rec_dN3.d3_N3_122_dnorm_f_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].d2_N3_122_dnorm_f_dN1_1;
        }
    }
    
    // Now perform finite difference for first-order derivatives
    order = 1;
    n = order + prec;
    switch(M2_State_Temp.finite_diff_type) {
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            exit(0);
            break;
    }
    
    rec_dN3.d_N3_111 = 0.0;
    rec_dN3.d_N3_122 = 0.0;
    rec_dN3.d_N3_123 = 0.0;
    
    rec_dN3.d2_N3_111_dgam_dN1_1 = 0.0;
    rec_dN3.d2_N3_122_dgam_dN1_1 = 0.0;
    
    rec_dN3.d2_N3_111_dnorm_f_dgam = 0.0;
    rec_dN3.d2_N3_122_dnorm_f_dgam = 0.0;
    rec_dN3.d3_N3_111_dnorm_f_dgam_dN1_1 = 0.0;
    rec_dN3.d3_N3_122_dnorm_f_dgam_dN1_1 = 0.0;
    
    for (int i = 0; i < n; i++) {
        id_c = n - i - 1;
        switch(M2_State_Temp.finite_diff_type) {
            case FORWARD_FINITE_DIFFERENCE:
                // forward finite difference approximation to first derivative
                id_vals = n - i - 1;
                break;
            case BACKWARD_FINITE_DIFFERENCE:
                // backward finite difference approximation to first derivative
                id_vals = nmax - i - 1;
                break;
            default:
                cout << "Finite difference type not specified" << endl;
                exit(0);
                break;
        }
        
        rec_dN3.d_N3_111 += c[id_c] * rec_N3_local_diff[id_vals].N3_111;
        rec_dN3.d_N3_122 += c[id_c] * rec_N3_local_diff[id_vals].N3_122;
        rec_dN3.d_N3_123 += c[id_c] * rec_N3_local_diff[id_vals].N3_123;
        
        rec_dN3.d2_N3_111_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].dN3_111_dN1_1;
        rec_dN3.d2_N3_122_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].dN3_122_dN1_1;
        
        if (flag_finite_diff_Sphere) {
            rec_dN3.d2_N3_111_dnorm_f_dgam += c[id_c] * rec_N3_local_diff[id_vals].dN3_111_dnorm_f;
            rec_dN3.d2_N3_122_dnorm_f_dgam += c[id_c] * rec_N3_local_diff[id_vals].dN3_122_dnorm_f;
            rec_dN3.d3_N3_111_dnorm_f_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].d2_N3_111_dnorm_f_dN1_1;
            rec_dN3.d3_N3_122_dnorm_f_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].d2_N3_122_dnorm_f_dN1_1;
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] rec_N3_local_diff;
}

void Calculate_Higher_Order_Moments_Derivatives_Vertex_Edges(record_Derivatives_Triangle &rec_dN3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere, const int &Edge_Type, const Finite_Diff_Data &rec_Finite_Diff) {
    M2_State_Param M2_State_Temp(M2_State);
    int Domain_Type;
    long double *c, *x;
    int order, prec;
    order = rec_Finite_Diff.order;
    prec = rec_Finite_Diff.prec;
    int n, nmax;
    n = order + prec;
    nmax = n;
    c = new long double[n];
    x = new long double[n];
    
    record_N3 *rec_N3_local_diff;
    rec_N3_local_diff = new record_N3[n];
    
    // Setup gamma knot
    M2_State_Temp.gam1_knot = M2_State.gamma_1;
    M2_State_Temp.gam2_knot = M2_State.gamma_2;
    
    M2_State_Temp.finite_diff_domain_Triangle = EDGE_TRIANGLE_FINITE_DIFFERENCE;
    M2_State_Temp.flag_finite_diff_Triangle = true;
    switch (M2_State.Domain_Type) {
        case BOUNDARY_GAM1_GAM2_EQ_0:
            switch (Edge_Type) {
                case EDGE_ONE:
                    M2_State_Temp.face_triangle = FACE_A; 
                    M2_State_Temp.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
                    M2_State_Temp.h_gam1 = 0.0;
                    M2_State_Temp.h_gam2 = rec_Finite_Diff.h;
                    M2_State_Temp.h_gam = sqrt(pow(M2_State_Temp.h_gam1, 2) + pow(M2_State_Temp.h_gam2, 2));
                    
                    // forward finite difference approximation to first derivative
                    differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
                    break;
                case EDGE_TWO:
                    M2_State_Temp.face_triangle = FACE_B;
                    M2_State_Temp.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
                    M2_State_Temp.h_gam1 = rec_Finite_Diff.h;
                    M2_State_Temp.h_gam2 = 0.0;
                    M2_State_Temp.h_gam = sqrt(pow(M2_State_Temp.h_gam1, 2) + pow(M2_State_Temp.h_gam2, 2));
                    
                    // forward finite difference approximation to first derivative
                    differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
                    break;
                default:
                    cout << "Edge type for finite differencing at edges of triangle not specified!!!!" << endl;
                    exit(0);
                    break;
            };
            break;
        case BOUNDARY_GAM1_GAM3_EQ_0:
            switch (Edge_Type) {
                case EDGE_ONE:
                    M2_State_Temp.face_triangle = FACE_A;
                    M2_State_Temp.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
                    M2_State_Temp.h_gam1 = 0.0;
                    M2_State_Temp.h_gam2 = -rec_Finite_Diff.h;
                    M2_State_Temp.h_gam = sqrt(pow(M2_State_Temp.h_gam1, 2) + pow(M2_State_Temp.h_gam2, 2));
                    
                    // backward finite difference approximation to first derivative
                    differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
                    M2_State_Temp.h_gam2 = fabs(M2_State_Temp.h_gam2);
                    break;
                case EDGE_TWO:
                    M2_State_Temp.face_triangle = FACE_C;
                    M2_State_Temp.finite_diff_type = FORWARD_FINITE_DIFFERENCE;
                    // If we go backwards along FACE_A then gam2 decreases
                    M2_State_Temp.h_gam1 = rec_Finite_Diff.h; //*sqrt(2.0)/2.0;
                    M2_State_Temp.h_gam2 = -rec_Finite_Diff.h; //*sqrt(2.0)/2.0;
                    M2_State_Temp.h_gam = sqrt(pow(M2_State_Temp.h_gam1, 2) + pow(M2_State_Temp.h_gam2, 2));
                    
                    // forward finite difference approximation to first derivative
                    differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
                    M2_State_Temp.h_gam1 = fabs(M2_State_Temp.h_gam1);
                    M2_State_Temp.h_gam2 = -fabs(M2_State_Temp.h_gam2);
                    break;
                default:
                    cout << "Edge type for finite differencing at edges of triangle not specified!!!!" << endl;
                    exit(0);
                    break;
            };
            break;
        case BOUNDARY_GAM2_GAM3_EQ_0:
            switch (Edge_Type) {
                case EDGE_ONE:
                    M2_State_Temp.face_triangle = FACE_B;
                    M2_State_Temp.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
                    M2_State_Temp.h_gam1 = -rec_Finite_Diff.h;
                    M2_State_Temp.h_gam2 = 0.0;
                    M2_State_Temp.h_gam = sqrt(pow(M2_State_Temp.h_gam1, 2) + pow(M2_State_Temp.h_gam2, 2));
                    
                    // backward finite difference approximation to first derivative
                    differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
                    M2_State_Temp.h_gam1 = fabs(M2_State_Temp.h_gam1);
                    break;
                case EDGE_TWO:
                    M2_State_Temp.finite_diff_type = BACKWARD_FINITE_DIFFERENCE;
                    M2_State_Temp.face_triangle = FACE_C;
                    // If we go backwards along FACE_C then gam1 decreases while gam2 increases
                    M2_State_Temp.h_gam1 = -rec_Finite_Diff.h; //*sqrt(2.0)/2.0;
                    M2_State_Temp.h_gam2 = rec_Finite_Diff.h; //*sqrt(2.0)/2.0;
                    M2_State_Temp.h_gam = sqrt(pow(M2_State_Temp.h_gam1, 2) + pow(M2_State_Temp.h_gam2, 2));
                    
                    // backward finite difference approximation to first derivative
                    differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
                    M2_State_Temp.h_gam1 = fabs(M2_State_Temp.h_gam1);
                    M2_State_Temp.h_gam2 = -fabs(M2_State_Temp.h_gam2);
                    break;
                default:
                    cout << "Edge type for finite differencing at edges of triangle not specified!!!!" << endl;
                    exit(0);
                    break;
            };
            break;
        default:
            cout << "Domain_Type for finite differencing at edges of triangle not specified!!!!" << endl;
            exit(0);
            break;
    };
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M2_State_Temp.id_finite_diff_gams = id_diff;
        
        M2_State_Temp.x_finite_diff_gam1 = x[id_diff]*M2_State_Temp.h_gam1/rec_Finite_Diff.h;
        M2_State_Temp.x_finite_diff_gam2 = x[id_diff]*M2_State_Temp.h_gam2/rec_Finite_Diff.h;
        
        Check_Domain_Type(Domain_Type, Var_index, num_points, M2_State_Temp);
        if ((Domain_Type != BOUNDARY_GAM1_GAM2_EQ_0) || (Domain_Type != BOUNDARY_GAM1_GAM3_EQ_0) || (Domain_Type != BOUNDARY_GAM2_GAM3_EQ_0)) {
            if (M2_State_Temp.Domain_Type != Domain_Type) {
                M2_State_Temp.Reset_Num_Vars(Domain_Type);
            }
        }
        
        NLOPT_Optim_Algo(rec_N3_local_diff, id_diff, M2_State_Temp, Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
    }
    
    // First perform finite difference for second-order derivatives
    rec_dN3.d2_N3_111 = 0.0;
    rec_dN3.d2_N3_122 = 0.0;
    rec_dN3.d2_N3_123 = 0.0;
    
    rec_dN3.d3_N3_111_dgam_dN1_1 = 0.0;
    rec_dN3.d3_N3_122_dgam_dN1_1 = 0.0;
    rec_dN3.d3_N3_111_dnorm_f_dgam = 0.0;
    rec_dN3.d3_N3_122_dnorm_f_dgam = 0.0;
    rec_dN3.d4_N3_111_dnorm_f_dgam_dN1_1 = 0.0;
    rec_dN3.d4_N3_122_dnorm_f_dgam_dN1_1 = 0.0;
    
    int id_c, id_vals;
    for (int i = 0; i < n; i++) {
        id_c = n - i - 1;
        id_vals = nmax - i - 1;
        
        rec_dN3.d_N3_111 += c[id_c] * rec_N3_local_diff[id_vals].N3_111;
        rec_dN3.d_N3_122 += c[id_c] * rec_N3_local_diff[id_vals].N3_122;
        rec_dN3.d_N3_123 += c[id_c] * rec_N3_local_diff[id_vals].N3_123;
        
        rec_dN3.d2_N3_111_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].dN3_111_dN1_1;
        rec_dN3.d2_N3_122_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].dN3_122_dN1_1;
        
        if (flag_finite_diff_Sphere) {
            rec_dN3.d2_N3_111_dnorm_f_dgam += c[id_c] * rec_N3_local_diff[id_vals].dN3_111_dnorm_f;
            rec_dN3.d2_N3_122_dnorm_f_dgam += c[id_c] * rec_N3_local_diff[id_vals].dN3_122_dnorm_f;
            rec_dN3.d3_N3_111_dnorm_f_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].d2_N3_111_dnorm_f_dN1_1;
            rec_dN3.d3_N3_122_dnorm_f_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].d2_N3_122_dnorm_f_dN1_1;
        }
    }
    
    // Now perform finite difference for first-order derivatives
    order = 1;
    n = order + prec;
    switch(M2_State_Temp.finite_diff_type) {
        case FORWARD_FINITE_DIFFERENCE:
            // forward finite difference approximation to first derivative
            differ_forward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        case BACKWARD_FINITE_DIFFERENCE:
            // backward finite difference approximation to first derivative
            differ_backward ( M2_State_Temp.h_gam, order, prec, c, x );
            break;
        default:
            cout << "Finite difference type not specified" << endl;
            exit(0);
            break;
    }
    
    rec_dN3.d_N3_111 = 0.0;
    rec_dN3.d_N3_122 = 0.0;
    rec_dN3.d_N3_123 = 0.0;
    
    rec_dN3.d2_N3_111_dgam_dN1_1 = 0.0;
    rec_dN3.d2_N3_122_dgam_dN1_1 = 0.0;
    
    rec_dN3.d2_N3_111_dnorm_f_dgam = 0.0;
    rec_dN3.d2_N3_122_dnorm_f_dgam = 0.0;
    rec_dN3.d3_N3_111_dnorm_f_dgam_dN1_1 = 0.0;
    rec_dN3.d3_N3_122_dnorm_f_dgam_dN1_1 = 0.0;
    
    for (int i = 0; i < n; i++) {
        id_c = n - i - 1;
        switch(M2_State_Temp.finite_diff_type) {
            case FORWARD_FINITE_DIFFERENCE:
                // forward finite difference approximation to first derivative
                id_vals = n - i - 1;
                break;
            case BACKWARD_FINITE_DIFFERENCE:
                // backward finite difference approximation to first derivative
                id_vals = nmax - i - 1;
                break;
            default:
                cout << "Finite difference type not specified" << endl;
                exit(0);
                break;
        }
        
        rec_dN3.d_N3_111 += c[id_c] * rec_N3_local_diff[id_vals].N3_111;
        rec_dN3.d_N3_122 += c[id_c] * rec_N3_local_diff[id_vals].N3_122;
        rec_dN3.d_N3_123 += c[id_c] * rec_N3_local_diff[id_vals].N3_123;
        
        rec_dN3.d2_N3_111_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].dN3_111_dN1_1;
        rec_dN3.d2_N3_122_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].dN3_122_dN1_1;
        
        if (flag_finite_diff_Sphere) {
            rec_dN3.d2_N3_111_dnorm_f_dgam += c[id_c] * rec_N3_local_diff[id_vals].dN3_111_dnorm_f;
            rec_dN3.d2_N3_122_dnorm_f_dgam += c[id_c] * rec_N3_local_diff[id_vals].dN3_122_dnorm_f;
            rec_dN3.d3_N3_111_dnorm_f_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].d2_N3_111_dnorm_f_dN1_1;
            rec_dN3.d3_N3_122_dnorm_f_dgam_dN1_1 += c[id_c] * rec_N3_local_diff[id_vals].d2_N3_122_dnorm_f_dN1_1;
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] rec_N3_local_diff;
}

void Calculate_Higher_Order_Moments_Derivatives_Triangle_Vertices(record_N3 &rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const int &id_proc, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, bool flag_finite_diff_Sphere) {
    // At the vertices of the triangle, we are interested in derivatives up to second-order
    switch (M2_State.Domain_Type) {
        case BOUNDARY_GAM1_GAM2_EQ_0:
        case BOUNDARY_GAM1_GAM3_EQ_0:
            M2_State.face_triangle = FACE_A;
            break;
        case BOUNDARY_GAM2_GAM3_EQ_0:
            M2_State.face_triangle = FACE_B;
            break;
        default:
            cout << "Domain Type for triangle finite difference not specified" << endl;
            exit(0);
            break;
    };
    
    // Perform finite difference approximation to obtain derivative
    long double u1, u2, v1, v2, w1, w2;
    long double c1_1, c1_2, c1_3, c2_1, c2_2, c2_3, c3_1, c3_2, c3_3;
    long double h;
    record_Derivatives_Triangle rec_dN3_median, rec_dN3_edge_1, rec_dN3_edge_2;
    Finite_Diff_Data rec_Finite_Diff;
    
    long double determinant, det1, det2, det3;
    rec_Finite_Diff.order = 2;
    rec_Finite_Diff.prec = 2;
    rec_Finite_Diff.h = finite_diff_h_gam;
    
    // First perform finite differencing along the line passing through the vertex
    // of interest as well as the barycenter of the triangle
    Calculate_Higher_Order_Moments_Derivatives_Vertex_Median(rec_dN3_median, M2_State, Var_index, num_points, id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, flag_finite_diff_Sphere, rec_Finite_Diff);
    
    // Now perform finite differencing along each of the two edges passing through the vertex
    Calculate_Higher_Order_Moments_Derivatives_Vertex_Edges(rec_dN3_edge_1, M2_State, Var_index, num_points, id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, flag_finite_diff_Sphere, EDGE_ONE, rec_Finite_Diff);
    
    Calculate_Higher_Order_Moments_Derivatives_Vertex_Edges(rec_dN3_edge_2, M2_State, Var_index, num_points, id_proc, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, flag_finite_diff_Sphere, EDGE_TWO, rec_Finite_Diff);
    
    // *********************************************************
    // Compute first-derivatives at the vertices
    // *********************************************************
    
    // If we know the directional derivative along two different directions
    // We can then solve the resulting system of two equations for the two
    // unknowns that are the derivatives along the Cartesian axes
    switch (M2_State.Domain_Type) {
        // Direction vector along the line defined by a gam1 + b gam2 + c = 0 (edge )
        // is u1 = {b, -a}
        case BOUNDARY_GAM1_GAM2_EQ_0:
        case BOUNDARY_GAM1_GAM3_EQ_0:
            // Edge 1  corresponds to FACE_A in this case
            // For FACE_A (a = 1, b = 0)
            u1 = 0.0;
            u2 = 1.0;
            break;
        case BOUNDARY_GAM2_GAM3_EQ_0:
            // Edge 1  corresponds to FACE_B in this case
            // For FACE_B (a = 0, b = 1)
            u1 = 1.0;
            u2 = 0.0;
            break;
        default:
            cout << "Domain Type for triangle finite difference not specified" << endl;
            exit(0);
            break;
    };
    // Direction vector along the line defined by a gam1 + b gam2 + c = 0 (median)
    // is u2 = {b, -a}
    // For median we have equation gam2 = a_median * gam1 + c ==> (a = -a_median, b = 1)
    Median_Line_Unit_Vector(v1, v2, M2_State);
    
    // u is the unitt vector along the edge of interest and v is the unit
    // vector along the line passing through the barycenter of the triangle
    // and the point on the edge of the triangle
    
    // We then have the following system of equations
    // Directional derivative along u
    // u1 * dN3ijk_dgam1 + u2 * dN3ijk_dgam2 = d_N3_ijk_edge
    // Directional derivative along v
    // v1 * dN3ijk_dgam1 + v2 * dN3ijk_dgam2 = d_N3_ijk_median
    // The solution is
    // dN3ijk_dgam1 = (v2 * d_N3_ijk_edge - u2 * d_N3_ijk_median)/(u1 * v2 - u2 * v1);
    // dN3ijk_dgam2 = (u1 * d_N3_ijk_median - v1 * d_N3_ijk_edge)/(u1 * v2 - u2 * v1);
    determinant = u1*v2 - u2*v1;
    
    rec_N3.dN3_111_dgam1 = (v2*rec_dN3_edge_1.d_N3_111 - u2*rec_dN3_median.d_N3_111)/determinant;
    rec_N3.dN3_122_dgam1 = (v2*rec_dN3_edge_1.d_N3_122 - u2*rec_dN3_median.d_N3_122)/determinant;
    rec_N3.dN3_123_dgam1 = (v2*rec_dN3_edge_1.d_N3_123 - u2*rec_dN3_median.d_N3_123)/determinant;
    
    rec_N3.d2_N3_111_dgam1_dN1_1 = (v2*rec_dN3_edge_1.d2_N3_111_dgam_dN1_1 - u2*rec_dN3_median.d2_N3_111_dgam_dN1_1)/determinant;
    rec_N3.d2_N3_122_dgam1_dN1_1 = (v2*rec_dN3_edge_1.d2_N3_122_dgam_dN1_1 - u2*rec_dN3_median.d2_N3_122_dgam_dN1_1)/determinant;
    
    rec_N3.d2_N3_111_dnorm_f_dgam1 = (v2*rec_dN3_edge_1.d2_N3_111_dnorm_f_dgam - u2*rec_dN3_median.d2_N3_111_dnorm_f_dgam)/determinant;
    rec_N3.d2_N3_122_dnorm_f_dgam1 = (v2*rec_dN3_edge_1.d2_N3_122_dnorm_f_dgam - u2*rec_dN3_median.d2_N3_122_dnorm_f_dgam)/determinant;
    
    rec_N3.d3_N3_111_dnorm_f_dgam1_dN1_1 = (v2*rec_dN3_edge_1.d3_N3_111_dnorm_f_dgam_dN1_1 - u2*rec_dN3_median.d3_N3_111_dnorm_f_dgam_dN1_1)/determinant;
    rec_N3.d3_N3_122_dnorm_f_dgam1_dN1_1 = (v2*rec_dN3_edge_1.d3_N3_122_dnorm_f_dgam_dN1_1 - u2*rec_dN3_median.d3_N3_122_dnorm_f_dgam_dN1_1)/determinant;
    
    rec_N3.dN3_111_dgam2 = (u1*rec_dN3_median.d_N3_111 - v1*rec_dN3_edge_1.d_N3_111)/determinant;
    rec_N3.dN3_122_dgam2 = (u1*rec_dN3_median.d_N3_122 - v1*rec_dN3_edge_1.d_N3_122)/determinant;
    rec_N3.dN3_123_dgam2 = (u1*rec_dN3_median.d_N3_123 - v1*rec_dN3_edge_1.d_N3_123)/determinant;
    
//     rec_N3.d2_N3_111_dgam2_dN1_1 = (u1*rec_dN3_median.d2_N3_111_dgam_dN1_1 - v1*rec_dN3_edge_1.d2_N3_111_dgam_dN1_1)/determinant;
    rec_N3.d2_N3_122_dgam2_dN1_1 = (u1*rec_dN3_median.d2_N3_122_dgam_dN1_1 - v1*rec_dN3_edge_1.d2_N3_122_dgam_dN1_1)/determinant;
    
//     rec_N3.d2_N3_111_dnorm_f_dgam2 = (u1*rec_dN3_median.d2_N3_111_dnorm_f_dgam - v1*rec_dN3_edge_1.d2_N3_111_dnorm_f_dgam)/determinant;
    rec_N3.d2_N3_122_dnorm_f_dgam2 = (u1*rec_dN3_median.d2_N3_122_dnorm_f_dgam - v1*rec_dN3_edge_1.d2_N3_122_dnorm_f_dgam)/determinant;
    
//     rec_N3.d3_N3_111_dnorm_f_dgam2_dN1_1 = (u1*rec_dN3_median.d3_N3_111_dnorm_f_dgam_dN1_1 - v1*rec_dN3_edge_1.d3_N3_111_dnorm_f_dgam_dN1_1)/determinant;
    rec_N3.d3_N3_122_dnorm_f_dgam2_dN1_1 = (u1*rec_dN3_median.d3_N3_122_dnorm_f_dgam_dN1_1 - v1*rec_dN3_edge_1.d3_N3_122_dnorm_f_dgam_dN1_1)/determinant;
            
    if (M2_State.display) {
        if (M2_State.id_proc == M2_State.proc_display) {
            // cout << "median = " << rec_dN3_median.d_N3_111 << "  " << "edge = " << rec_dN3_edge_1.d_N3_111 << "  " << "median = " << rec_dN3_median.d_N3_122 << "  " << "edge = " << rec_dN3_edge_1.d_N3_122 << "  " << "dN3_111_dgam1 = " << rec_N3.dN3_111_dgam1 << "  " << "dN3_111_dgam2 = " << rec_N3.dN3_111_dgam2 << "  " << "dN3_122_dgam1 = " << rec_N3.dN3_122_dgam1 << "  " << "dN3_122_dgam2 = " << rec_N3.dN3_122_dgam2 << "  " << "dN3_123_dgam1 = " << rec_N3.dN3_123_dgam1 << "  " << "dN3_123_dgam2 = " << rec_N3.dN3_123_dgam2 << endl;
        }
    }
    
    
//     cout << "*********************************************************************" << endl;
//     cout << "u1 = " << u1 << "  " << "u2 = " << u2 << endl;
//     cout << "v1 = " << v1 << "  " << "v2 = " << v2 << endl;
//     cout << "AAAAAAAAAAAAAAAAAA" << "dN3_111_dgam1 = " << rec_N3.dN3_111_dgam1 << "  " << "edge_1 d_N3_111 = " << rec_dN3_edge_1.d_N3_111 << "  " << "edge_2 d_N3_111 = " << rec_dN3_edge_2.d_N3_111 << "  " << "median d_N3_111 = " << rec_dN3_median.d_N3_111 << endl;
//     cout << "*********************************************************************" << endl;
//     exit(0);
    
    // *********************************************************
    // Compute second-derivatives at the vertices
    // *********************************************************
    
    // Website for directional derivatives:
    // http://mathonline.wikidot.com/higher-order-directional-derivatives
    
    // If we know the directional second derivative along three different directions
    // We can then solve the resulting system of three equations for the three
    // unknowns that are the second derivatives along the Cartesian axes as well 
    // the mixed derivative
    switch (M2_State.Domain_Type) {
        // Direction vector along the line defined by a gam1 + b gam2 + c = 0 (edge )
        // is u1 = {b, -a}
        case BOUNDARY_GAM1_GAM2_EQ_0:
            // For FACE_A (a = 1, b = 0)
            u1 = 0.0;
            u2 = 1.0;
            // For FACE_B (a = 0, b = 1)
            v1 = 1.0;
            v2 = 0.0;
            break;
        case BOUNDARY_GAM1_GAM3_EQ_0:
            // For FACE_A (a = 1, b = 0)
            u1 = 0.0;
            u2 = 1.0;
            // For FACE_C (a = 1, b = 1)
            v1 = 1.0/sqrt(2.0);
            v2 = -1.0/sqrt(2.0);
            break;
        case BOUNDARY_GAM2_GAM3_EQ_0:
            // For FACE_B (a = 0, b = 1)
            u1 = 1.0;
            u2 = 0.0;
            // For FACE_C (a = 1, b = 1)
            v1 = 1.0/sqrt(2.0);
            v2 = -1.0/sqrt(2.0);
            break;
        default:
            cout << "Domain Type for triangle finite difference not specified" << endl;
            exit(0);
            break;
    };
    
    // Direction vector along the line defined by a gam1 + b gam2 + c = 0 (median)
    // is u2 = {b, -a}
    // For median we have equation gam2 = a_median * gam1 + c ==> (a = -a_median, b = 1)
    Median_Line_Unit_Vector(w1, w2, M2_State);
    
    c1_1 = w1*w1;
    c1_2 = 2.0*w1*w2;
    c1_3 = w2*w2;
    
    c2_1 = u1*u1;
    c2_2 = 2.0*u1*u1;
    c2_3 = u2*u2;
    
    c3_1 = v1*v1;
    c3_2 = 2.0*v1*v2;
    c3_3 = v2*v2;
    
    // We then have the following system of equations
    // Directional second derivative along u1
    // c1_1 * d2N3ijk_dgam1^2 + c1_2 * d2N3ijk_dgam1_dgam2 + c1_3 * d2N3ijk_dgam2^2 = d2_N3_ijk_median
    // Directional derivative along u2
    // c2_1 * d2N3ijk_dgam1^2 + c2_2 * d2N3ijk_dgam1_dgam2 + c2_3 * d2N3ijk_dgam2^2 = d2_N3_ijk_edge_1
    // Directional derivative along u2
    // c3_1 * d2N3ijk_dgam1^2 + c3_2 * d2N3ijk_dgam1_dgam2 + c3_3 * d2N3ijk_dgam2^2 = d2_N3_ijk_edge_2
    determinant = c1_1*c2_2*c3_3 - c1_1*c2_3*c3_2 - c1_2*c2_1*c3_3 + c1_2*c2_3*c3_1 + c1_3*c2_1*c3_2 - c1_3*c2_2*c3_1;
    
    det1 = c2_1*c3_3 - c2_3*c3_1;
    det2 = c1_1*c3_3 - c1_3*c3_1;
    det3 = c1_1*c2_3 - c1_3*c2_1;
    
//     rec_N3.d2_N3_111_dgam1_2 = (c2_2*c3_3 - c2_3*c3_2)*d2_N3_111_median - (c1_2*c3_3 - c1_3*c3_2)*d2_N3_111_edge_1 + (c1_2*c2_3 - c1_3*c2_2)*d2_N3_111_edge_2;
//     rec_N3.d2_N3_111_dgam1_2 /= determinant;
    
//     rec_N3.d2_N3_111_dgam2_2 = (c2_1*c3_2 - c2_2*c3_1)*d2_N3_111_median - (c1_1*c3_2 - c1_2*c3_1)*d2_N3_111_edge_1 + (c1_1*c2_2 - c1_2*c2_1)*d2_N3_111_edge_2;
//     rec_N3.d2_N3_111_dgam2_2 /= determinant;
    
    rec_N3.d2_N3_122_dgam1_dgam2 = -det1*rec_dN3_median.d2_N3_122 + det2*rec_dN3_edge_1.d2_N3_122 - det3*rec_dN3_edge_2.d2_N3_122;
    rec_N3.d2_N3_122_dgam1_dgam2 /= determinant;
    
    rec_N3.d3_N3_122_dgam1_dgam2_dN1_1 = -det1*rec_dN3_median.d3_N3_122_dgam_dN1_1 + det2*rec_dN3_edge_1.d3_N3_122_dgam_dN1_1 - det3*rec_dN3_edge_2.d3_N3_122_dgam_dN1_1;
    rec_N3.d3_N3_122_dgam1_dgam2_dN1_1 /= determinant;
    
    rec_N3.d3_N3_122_dnorm_f_dgam1_dgam2 = -det1*rec_dN3_median.d3_N3_122_dnorm_f_dgam + det2*rec_dN3_edge_1.d3_N3_122_dnorm_f_dgam - det3*rec_dN3_edge_2.d3_N3_122_dnorm_f_dgam;
    rec_N3.d3_N3_122_dnorm_f_dgam1_dgam2 /= determinant;
    
    rec_N3.d4_N3_122_dnorm_f_dgam1_dgam2_dN1_1 = -det1*rec_dN3_median.d4_N3_122_dnorm_f_dgam_dN1_1 + det2*rec_dN3_edge_1.d4_N3_122_dnorm_f_dgam_dN1_1- det3*rec_dN3_edge_2.d4_N3_122_dnorm_f_dgam_dN1_1;
    rec_N3.d4_N3_122_dnorm_f_dgam1_dgam2_dN1_1 /= determinant;
    
    
    rec_N3.d2_N3_123_dgam1_dgam2 = -det1*rec_dN3_median.d2_N3_123 + det2*rec_dN3_edge_1.d2_N3_123- det3*rec_dN3_edge_2.d2_N3_123;
    rec_N3.d2_N3_123_dgam1_dgam2 /= determinant;
    
    rec_N3.d3_N3_123_dgam1_dgam2_dN1_1 = -det1*rec_dN3_median.d3_N3_123_dgam_dN1_1 + det2*rec_dN3_edge_1.d3_N3_123_dgam_dN1_1 - det3*rec_dN3_edge_2.d3_N3_123_dgam_dN1_1;
    rec_N3.d3_N3_123_dgam1_dgam2_dN1_1 /= determinant;
    
    rec_N3.d3_N3_123_dN1_2_dgam1_dgam2 = -det1*rec_dN3_median.d3_N3_123_dN1_2_dgam + det2*rec_dN3_edge_1.d3_N3_123_dN1_2_dgam - det3*rec_dN3_edge_2.d3_N3_123_dN1_2_dgam;
    rec_N3.d3_N3_123_dN1_2_dgam1_dgam2 /= determinant;
    
    rec_N3.d3_N3_123_dN1_3_dgam1_dgam2 = -det1*rec_dN3_median.d3_N3_123_dN1_3_dgam + det2*rec_dN3_edge_1.d3_N3_123_dN1_3_dgam - det3*rec_dN3_edge_2.d3_N3_123_dN1_3_dgam;
    rec_N3.d3_N3_123_dN1_3_dgam1_dgam2 /= determinant;
    
    rec_N3.d4_N3_123_dN1_1_dN1_2_dgam1_dgam2 = -det1*rec_dN3_median.d4_N3_123_dN1_1_dN1_2_dgam + det2*rec_dN3_edge_1.d4_N3_123_dN1_1_dN1_2_dgam - det3*rec_dN3_edge_2.d4_N3_123_dN1_1_dN1_2_dgam;
    rec_N3.d4_N3_123_dN1_1_dN1_2_dgam1_dgam2 /= determinant;
    
    rec_N3.d4_N3_123_dN1_1_dN1_3_dgam1_dgam2 = -det1*rec_dN3_median.d4_N3_123_dN1_1_dN1_3_dgam + det2*rec_dN3_edge_1.d4_N3_123_dN1_1_dN1_3_dgam - det3*rec_dN3_edge_2.d4_N3_123_dN1_1_dN1_3_dgam;
    rec_N3.d4_N3_123_dN1_1_dN1_3_dgam1_dgam2 /= determinant;
    
    rec_N3.d4_N3_123_dN1_2_dN1_3_dgam1_dgam2 = -det1*rec_dN3_median.d4_N3_123_dN1_2_dN1_3_dgam + det2*rec_dN3_edge_1.d4_N3_123_dN1_2_dN1_3_dgam - det3*rec_dN3_edge_2.d4_N3_123_dN1_2_dN1_3_dgam;
    rec_N3.d4_N3_123_dN1_2_dN1_3_dgam1_dgam2 /= determinant;
    
    rec_N3.d5_N3_123_dN1_1_dN1_2_dN1_3_dgam1_dgam2 = -det1*rec_dN3_median.d5_N3_123_dN1_1_dN1_2_dN1_3_dgam + det2*rec_dN3_edge_1.d5_N3_123_dN1_1_dN1_2_dN1_3_dgam - det3*rec_dN3_edge_2.d5_N3_123_dN1_1_dN1_2_dN1_3_dgam;
    rec_N3.d5_N3_123_dN1_1_dN1_2_dN1_3_dgam1_dgam2 /= determinant;
    
    if (M2_State.display) {
        if (M2_State.id_proc == M2_State.proc_display) {
            //         cout << "median = " << d2_N3_111_median << "  " << "edge = " << d2_N3_111_edge_1 << "  " << "dN3_111_dgam1 = " << rec_N3.dN3_111_dgam1 << "  " << "dN3_111_dgam2 = " << rec_N3.dN3_111_dgam2 << "  " << "dN3_122_dgam1 = " << rec_N3.dN3_122_dgam1 << "  " << "dN3_122_dgam2 = " << rec_N3.dN3_122_dgam2 << "  " << "dN3_123_dgam1 = " << rec_N3.dN3_123_dgam1 << "  " << "dN3_123_dgam2 = " << rec_N3.dN3_123_dgam2 << endl;
        }
    }
}
