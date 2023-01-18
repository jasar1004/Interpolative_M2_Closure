#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#include "NG_MN_Model_3D_OPTIM.h"
#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED

#include <mpi.h>

long double L_vals[11] = {1.0e-6, 1.0e-1, 0.5, 1.0, 5.0, 1.0e1, 5.0e1, 1.0e2, 5.0e2, 1.0e3, 1.0e4};
long double E_vals[11] = {1.0e-6, 1.0e-1, 0.5, 1.0, 5.0, 1.0e1, 5.0e1, 1.0e2, 5.0e2, 1.0e3, 1.0e4};

long double r_l[4] = {0.0, 1.0e-8, 1.0e-6, 1.0e-4};
int Lebed_Rule_Set[6] = {10, 20, 32, 44, 50, 65};
int N_quad_points_Circle_Set[4] = {30, 64, 127, 129}; //{5, 10, 20, 64};//, 127, 255};
int M2_State_Param::Regime = GRAY;

const int M2_State_Param :: max_refinement_level;

int M2_State_Param::Max_Ent_Solution_Type = CLOSING_FLUX;

int DISPLAY_ID = 0;
int PRIMARY_ID = 0;

long double finite_diff_h_N1 = 5.0*1.0e-4;
long double finite_diff_h_gam = 5.0*1.0e-4;
long double tol_grad = 1.0*1.0e-4;

// Consider effect of constraints scaling and objective function tolerance on the
// convergence of the optimization algorithm
// Also consider the effect of maximum number of iterations of the optimization algorithm
// on the convergence, especially near the boundaries of the realizable space
                           
int Optimization_Algorithm(M2_State_Param &M2_State,
                           MN_Var_Num_Points *num_points, 
                           fstream &in_N3_out, 
                           const int &N_pts_Mob_Scale, 
                           const long double *Coefficients_Mobius_Scale_Fit) {
    record_N3 *rec_N3_local, *rec_N3_global;
    MPI_proc_params MPI_proc_parameters;
    MPI_Datatype rec_N3_type;
    
    num_points->N_Points_Triangle_gam1_gam2 = (num_points->gam1*(num_points->gam1+1))/2;
    
     // Create the rec_N3 datatype for MPI
    Create_MPI_Data_Type_rec_N3(rec_N3_type);
    
    // Setup parameters required for the purpose of parallel computing
    Setup_MPI_Processes(MPI_proc_parameters, M2_State, num_points);
    
    int Npts_total = Compute_Npts_Total(num_points);
    if (M2_State.id_proc == 0) {
        rec_N3_global = new record_N3[Npts_total];
    }
    
    rec_N3_local = new record_N3[MPI_proc_parameters.size_rec];
    
    // Call routine for entropy optimization at the points of interest
    Optimization_Algorithm_Array(rec_N3_local, M2_State, num_points, &MPI_proc_parameters, N_pts_Mob_Scale, Coefficients_Mobius_Scale_Fit);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Gather(rec_N3_local, MPI_proc_parameters.size_rec, rec_N3_type, rec_N3_global, MPI_proc_parameters.size_rec, rec_N3_type, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    int index;
    int num_proc_per_var_E = MPI_proc_parameters.num_proc_per_var_E;
    int num_proc_per_var_f = MPI_proc_parameters.num_proc_per_var_f;
    int num_proc_per_var_phi = MPI_proc_parameters.num_proc_per_var_phi;
    int num_proc_per_var_theta = MPI_proc_parameters.num_proc_per_var_theta;
    int num_proc_per_var_Triangle_gam1_gam2 = MPI_proc_parameters.num_proc_per_var_Triangle_gam1_gam2;
    
    if (M2_State.id_proc == 0) {
        if (M2_State.display) { 
            // cout << "id_count = " << id_count << endl;
        }
        for (int i_proc_E = 0; i_proc_E < num_proc_per_var_E; i_proc_E++) {
            for (int id_E = 0; id_E < num_points->E/num_proc_per_var_E; id_E++) {
                for (int i_proc_f = 0; i_proc_f < num_proc_per_var_f; i_proc_f++) {
                    for (int id_f = 0; id_f < num_points->f/num_proc_per_var_f; id_f++) {
                        for (int i_proc_phi = 0; i_proc_phi < num_proc_per_var_phi; i_proc_phi++) {
                            for (int id_phi = 0; id_phi < num_points->phi/num_proc_per_var_phi; id_phi++) {
                                for (int i_proc_theta = 0; i_proc_theta < num_proc_per_var_theta; i_proc_theta++) {
                                    for (int Id_mu = 0; Id_mu < num_points->theta/num_proc_per_var_theta; Id_mu++) {
                                        for (int i_proc_Triangle_gam1_gam2 = 0; i_proc_Triangle_gam1_gam2 < num_proc_per_var_Triangle_gam1_gam2; i_proc_Triangle_gam1_gam2++) {
                                            for (int Id_gam1_gam2 = 0; Id_gam1_gam2 < num_points->N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2; Id_gam1_gam2++) {
                                                    index = i_proc_E*num_proc_per_var_f + i_proc_f;
                                                    index = index*num_proc_per_var_phi + i_proc_phi;
                                                    index = index*num_proc_per_var_theta + i_proc_theta;
                                                    index = index*num_proc_per_var_Triangle_gam1_gam2 + i_proc_Triangle_gam1_gam2;
                                                    index = index*MPI_proc_parameters.size_rec;
                                                    
                                                    index = index + (((id_E*num_points->f/num_proc_per_var_f + id_f)*num_points->phi/num_proc_per_var_phi + id_phi)*num_points->theta/num_proc_per_var_theta + Id_mu)*(num_points->N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2) + Id_gam1_gam2;
                                                    write<record_N3>(in_N3_out, rec_N3_global[index]);
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
        if (M2_State.display) { 
            cout << "***************************** Writing to file Completed *****************************" << endl;
        }
    }
    
    Finalize:;
    
    if (M2_State.id_proc == 0) {
        delete[] rec_N3_global;
    }
    
    delete[] rec_N3_local;

    return 0;
}

int Optimization_Algorithm_Array(record_N3 *rec_N3, M2_State_Param &M2_State, const MN_Var_Num_Points *num_points, const MPI_proc_params *MPI_proc_parameters, const int &N_pts_Mob_Scale, const long double *Coefficients_Mobius_Scale_Fit) {
    int Domain_Type = 0;
    MN_Var_Num_Points Var_index;
    long double *x_Lebed = NULL;
    long double *y_Lebed = NULL;
    long double *z_Lebed = NULL;
    
    if (M2_State.Node_Dist_Phi_Theta == SPHERICAL_HARMONIC_DISTRIBUTION && M2_State.Rec_Moms_Test.flag_Moments_Test != true) {
//         if (num_points->phi != 2*num_points->theta) {
//             cout << "Spherical Harmonics Distribution not correct" << endl;
//             exit(0);
//         }
        x_Lebed = new long double[num_points->phi*num_points->theta];
        y_Lebed = new long double[num_points->phi*num_points->theta];
        z_Lebed = new long double[num_points->phi*num_points->theta];
        setup_spherical_harmonics_data ( num_points->phi, num_points->theta, x_Lebed, y_Lebed, z_Lebed, NULL );
    } else {
        x_Lebed = new long double[num_points->phi*num_points->theta];
        y_Lebed = new long double[num_points->phi*num_points->theta];
        z_Lebed = new long double[num_points->phi*num_points->theta];
    }
    
    M2_State_Param M2_State_Full_Triangle(M2_State.Dimension, GAM1_GAM2_GAM3, M2_State.Problem_Type);
    M2_State_Full_Triangle.Node_Dist_E = M2_State.Node_Dist_E;
    M2_State_Full_Triangle.Node_Dist_f = M2_State.Node_Dist_f;
    M2_State_Full_Triangle.Node_Dist_Phi_Theta = M2_State.Node_Dist_Phi_Theta;
    M2_State_Full_Triangle.Node_Dist_gam1 = M2_State.Node_Dist_gam1;
    M2_State_Full_Triangle.Triangle_Domain = M2_State.Triangle_Domain;
    M2_State_Full_Triangle.display = M2_State.display;
    M2_State_Full_Triangle.proc_display = M2_State.proc_display;
    M2_State_Full_Triangle.id_proc = M2_State.id_proc;
    M2_State_Full_Triangle.Rec_Moms_Test.Copy(M2_State.Rec_Moms_Test);
    
    M2_State_Param M2_State_Boundary_Gam1(M2_State.Dimension, BOUNDARY_GAM1, M2_State.Problem_Type);
    M2_State_Boundary_Gam1.Node_Dist_E = M2_State.Node_Dist_E;
    M2_State_Boundary_Gam1.Node_Dist_f = M2_State.Node_Dist_f;
    M2_State_Boundary_Gam1.Node_Dist_Phi_Theta = M2_State.Node_Dist_Phi_Theta;
    M2_State_Boundary_Gam1.Node_Dist_gam1 = M2_State.Node_Dist_gam1;
    M2_State_Boundary_Gam1.Triangle_Domain = M2_State.Triangle_Domain;
    M2_State_Boundary_Gam1.display = M2_State.display;
    M2_State_Boundary_Gam1.proc_display = M2_State.proc_display;
    M2_State_Boundary_Gam1.id_proc = M2_State.id_proc;
    M2_State_Boundary_Gam1.Rec_Moms_Test.Copy(M2_State.Rec_Moms_Test);
    
    M2_State_Param M2_State_Boundary_Gam2(M2_State.Dimension, BOUNDARY_GAM2, M2_State.Problem_Type);
    M2_State_Boundary_Gam2.Node_Dist_E = M2_State.Node_Dist_E;
    M2_State_Boundary_Gam2.Node_Dist_f = M2_State.Node_Dist_f;
    M2_State_Boundary_Gam2.Node_Dist_Phi_Theta = M2_State.Node_Dist_Phi_Theta;
    M2_State_Boundary_Gam2.Node_Dist_gam1 = M2_State.Node_Dist_gam1;
    M2_State_Boundary_Gam2.Triangle_Domain = M2_State.Triangle_Domain;
    M2_State_Boundary_Gam2.display = M2_State.display;
    M2_State_Boundary_Gam2.proc_display = M2_State.proc_display;
    M2_State_Boundary_Gam2.id_proc = M2_State.id_proc;
    M2_State_Boundary_Gam2.Rec_Moms_Test.Copy(M2_State.Rec_Moms_Test);
    
    M2_State_Param M2_State_Boundary_Gam3(M2_State.Dimension, BOUNDARY_GAM3, M2_State.Problem_Type);
    M2_State_Boundary_Gam3.Node_Dist_E = M2_State.Node_Dist_E;
    M2_State_Boundary_Gam3.Node_Dist_f = M2_State.Node_Dist_f;
    M2_State_Boundary_Gam3.Node_Dist_Phi_Theta = M2_State.Node_Dist_Phi_Theta;
    M2_State_Boundary_Gam3.Node_Dist_gam1 = M2_State.Node_Dist_gam1;
    M2_State_Boundary_Gam3.Triangle_Domain = M2_State.Triangle_Domain;
    M2_State_Boundary_Gam3.display = M2_State.display;
    M2_State_Boundary_Gam3.proc_display = M2_State.proc_display;
    M2_State_Boundary_Gam3.id_proc = M2_State.id_proc;
    M2_State_Boundary_Gam3.Rec_Moms_Test.Copy(M2_State.Rec_Moms_Test);
    
    M2_State_Param M2_State_Boundary_Point(M2_State.Dimension, GAM1_GAM2_GAM3, M2_State.Problem_Type);
    M2_State_Boundary_Point.Node_Dist_E = M2_State.Node_Dist_E;
    M2_State_Boundary_Point.Node_Dist_f = M2_State.Node_Dist_f;
    M2_State_Boundary_Point.Node_Dist_Phi_Theta = M2_State.Node_Dist_Phi_Theta;
    M2_State_Boundary_Point.Node_Dist_gam1 = M2_State.Node_Dist_gam1;
    M2_State_Boundary_Point.Triangle_Domain = M2_State.Triangle_Domain;
    M2_State_Boundary_Point.display = M2_State.display;
    M2_State_Boundary_Point.proc_display = M2_State.proc_display;
    M2_State_Boundary_Point.id_proc = M2_State.id_proc;
    M2_State_Boundary_Point.Rec_Moms_Test.Copy(M2_State.Rec_Moms_Test);
    
    int Num_pts_f_Mob_Scale, Order_SH_Mob_Scale, Num_pts_SH_Mob_Scale, Num_pts_gam1_Mob_Scale;
    if (num_points->Length_Scale_Dist_Type == LENGTH_SCALE_DIST_UNIF) {
        Num_pts_f_Mob_Scale = 1;
        Order_SH_Mob_Scale = 0;
        Num_pts_SH_Mob_Scale = 1;
        Num_pts_gam1_Mob_Scale = 1;
    } else if (num_points->Length_Scale_Dist_Type == LENGTH_SCALE_DIST_FIT) {
        Num_pts_f_Mob_Scale = num_points->f;
        Order_SH_Mob_Scale = num_points->Order_SH;
        Num_pts_SH_Mob_Scale = num_points->N_pts_SH;
        Num_pts_gam1_Mob_Scale = num_points->gam1;
    } else {
        cout << "Length_Scale_Dist_Type not spceified !!!!!!!!!!!!!!!!!!" << endl;
        exit(0);
    }
    
    Mobius_Scale_Parameters Mobius_Scale_Params(Num_pts_f_Mob_Scale, Order_SH_Mob_Scale, Num_pts_SH_Mob_Scale, Num_pts_gam1_Mob_Scale, Coefficients_Mobius_Scale_Fit);
    Mobius_Scale_Params.Length_Scale_Dist_Type = num_points->Length_Scale_Dist_Type;
    
    int num_proc_per_var_E, num_proc_per_var_f, num_proc_per_var_theta, num_proc_per_var_phi;
    int id_proc_E, id_proc_f, id_proc_theta, id_proc_phi;
    
    num_proc_per_var_E = MPI_proc_parameters->num_proc_per_var_E;
    num_proc_per_var_f = MPI_proc_parameters->num_proc_per_var_f;
    num_proc_per_var_phi = MPI_proc_parameters->num_proc_per_var_phi;
    num_proc_per_var_theta = MPI_proc_parameters->num_proc_per_var_theta;
    
    id_proc_E = MPI_proc_parameters->id_proc_E;
    id_proc_f = MPI_proc_parameters->id_proc_f;
    id_proc_phi = MPI_proc_parameters->id_proc_phi;
    id_proc_theta = MPI_proc_parameters->id_proc_theta;
    
    int id_count = 0;
    int id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, Id_mu_min, Id_mu_max, Id_Triangle_gam1_gam2_min, Id_Triangle_gam1_gam2_max;
    
    Compute_MPI_Processes_Max_Min_Indexes(*MPI_proc_parameters, num_points, id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, Id_mu_min, Id_mu_max, Id_Triangle_gam1_gam2_min, Id_Triangle_gam1_gam2_max);
    
    if (M2_State_Full_Triangle.id_proc == M2_State_Full_Triangle.proc_display) {
        cout << "id_E_min = " << id_E_min << "  " << "id_E_max = " << id_E_max << "  " << "id_f_min = " << id_f_min << "  " << "id_f_max = " << id_f_max << "  " << "id_phi_min = " << id_phi_min << "  " << "id_phi_max = " << id_phi_max << "  " << "Id_mu_min = " << Id_mu_min << "  " << "Id_mu_max = " << Id_mu_max << "  " << "Id_Triangle_gam1_gam2_min = " << Id_Triangle_gam1_gam2_min << "  " << "Id_Triangle_gam1_gam2_max = " << Id_Triangle_gam1_gam2_max << endl;
    }
    // cout << "M2_State.display = " << M2_State.display << "  " << "M2_State.proc_display = " << M2_State.proc_display << endl;
    
    int i_gam1_Square, i_gam2_Square;
    
    for (int index_e = id_E_min; index_e < id_E_max; index_e++){
        for (int i_f = id_f_min; i_f < id_f_max; i_f++) {
            for (int i_phi = id_phi_min; i_phi < id_phi_max; i_phi++) {
                for (int i_theta = Id_mu_min; i_theta < Id_mu_max; i_theta++) {
                    for (int i_Triangle_gam1_gam2 = Id_Triangle_gam1_gam2_min; i_Triangle_gam1_gam2 < Id_Triangle_gam1_gam2_max; i_Triangle_gam1_gam2++) {
                            Compute_Id_gam1_gam2_Triangle(i_gam1_Square, i_gam2_Square, i_Triangle_gam1_gam2, num_points); 
                            // index_e = 2;
                            // i_f = 2;
                            // i_phi = 1;
                            // i_theta = 1;
                            // i_gam1_Square = 0;
                            // i_gam2_Square = 0;
                            
//                             N1_1 = 0.333333   N1_2 = 0   N1_3 = 0   gam1 = 0.142857   gam2 = 3.70074e-17   N3_122 = 0  f_N3_122 = 1.70858  dN3_122 = 0.368617      
                        
                            Var_index.E = index_e;
                            Var_index.f = i_f;
                            Var_index.phi = i_phi;
                            Var_index.theta = i_theta;
                            Var_index.gam1 = i_gam1_Square;
                            Var_index.gam2 = i_gam2_Square;
                            
                            Check_Domain_Type(Domain_Type, &Var_index, num_points, M2_State_Full_Triangle);
                            
                            if (M2_State_Full_Triangle.flag_Realizability) {
                                if (M2_State_Full_Triangle.display) { 
                                    cout << ".......................Non Realizable Moments......................" << endl;
                                    cout << "index_e = " << Var_index.E << "     "  << "i_f = " << Var_index.f << "     "  << "i_phi = " << Var_index.phi << "     "  << "i_theta = " << Var_index.theta << "     "  << "i_gam1_Square = " << Var_index.gam1 << "     "  << "i_gam2_Square = " << Var_index.gam2 << endl;
                                    
                                    cout  << "I0 = " << M2_State_Full_Triangle.I0 << "   " << "N1_1 = " << M2_State_Full_Triangle.N1_1 << "    " << "N1_2 = " << M2_State_Full_Triangle.N1_2 << "   " << "N1_3 = " << M2_State_Full_Triangle.N1_3 << "     " << "gamma_1 = " << M2_State_Full_Triangle.gamma_1 << "    " << "gamma_2 = " << M2_State_Full_Triangle.gamma_2 << "    " << "id = " << M2_State_Full_Triangle.id_proc << "     "  << "Domain_Type = " << M2_State_Full_Triangle.Domain_Type << endl; 
                                }
                                exit(0);
                            }
    
                            switch (Domain_Type) {
                                case GAM1_GAM2_GAM3:
                                    NLOPT_Optim_Algo(rec_N3, id_count, M2_State_Full_Triangle, &Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
                                    break;
                                case BOUNDARY_GAM1:
                                    NLOPT_Optim_Algo(rec_N3, id_count, M2_State_Boundary_Gam1, &Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
                                    break;
                                case BOUNDARY_GAM2:
                                    NLOPT_Optim_Algo(rec_N3, id_count, M2_State_Boundary_Gam2, &Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
                                    break;
                                case BOUNDARY_GAM3:
                                    NLOPT_Optim_Algo(rec_N3, id_count, M2_State_Boundary_Gam3, &Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
                                    break;
                                default:
                                    M2_State_Boundary_Point.Domain_Type = Domain_Type;
                                    NLOPT_Optim_Algo(rec_N3, id_count, M2_State_Boundary_Point, &Var_index, num_points, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed);
                                    break;
                            }
                            id_count++;
                    }
                }
            }
        }
    }
    
    if (id_count != MPI_proc_parameters->size_rec) {
        cout << "id_count = " << id_count << "   " << "size_rec = " << MPI_proc_parameters->size_rec << endl;
        exit(0);
    }
     
    Finalize:;
    
    if (x_Lebed != NULL) {
        delete[] x_Lebed; x_Lebed = NULL;
    }
    if (y_Lebed != NULL) {
        delete[] y_Lebed; y_Lebed = NULL;
    }
    if (z_Lebed != NULL) {
        delete[] z_Lebed; z_Lebed = NULL;
    }

    return 0;
}

void NLOPT_Optim_Algo(record_N3 *rec_N3_local, const int &id_count, M2_State_Param &M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, Mobius_Scale_Parameters &Mobius_Scale_Params, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed) {
    long double *x, *tol_x, *Sk_final;
    long double minf = 0.0, min_grad = 0.0; /* the minimum objective value, upon return */
    long double norm_f_2;
    int max_regularization = 4;
    int max_iters_nlopt;
    
    max_iters_nlopt = 1500;
    
    nlopt_opt opt;
    opt = NULL;
    
    x = new long double[M2_State.NVARS];
    tol_x = new long double[M2_State.NVARS];
    Sk_final = new long double [M2_State.NVARS*M2_State.NVARS];
    
    for (int i_Sk = 0; i_Sk < M2_State.NVARS; i_Sk++) {
        for (int j_Sk = 0; j_Sk < M2_State.NVARS; j_Sk++) {
            Sk_final[i_Sk*M2_State.NVARS + j_Sk] = 0.0;
            if (i_Sk == j_Sk) {
                Sk_final[i_Sk*M2_State.NVARS + j_Sk] = 1.0;
            }
        }
    }
    
    // cout << "M2_State.display = " << M2_State.display << "  " << "M2_State.proc_display = " << M2_State.proc_display << endl;
    
    if (M2_State.display) { 
        if (M2_State.id_proc == M2_State.proc_display) {
            cout << "index_e = " << Var_index->E << "     "  << "i_f = " << Var_index->f << "     "  << "i_phi = " << Var_index->phi << "     "  << "i_theta = " << Var_index->theta << "     "  << "i_gam1_Square = " << Var_index->gam1 << "     "  << "i_gam2_Square = " << Var_index->gam2 << endl;
            if (M2_State.flag_finite_diff_Sphere || M2_State.flag_finite_diff_Triangle) {
                if (M2_State.flag_finite_diff_Sphere) {
                    cout << "Finite difference Sphere, i_diff = " << M2_State.id_finite_diff_norm_f << endl;
                    cout << endl;
                } else {
                    cout << "Finite difference triangle i_diff = " << M2_State.id_finite_diff_gams << endl;
                }
            }
        }
    }
            
    // Case where we are either in the free streaming limit or at one of the vertices of the triangle (P1 P2 P3)
    if ((M2_State.Domain_Type != BOUNDARY_GAM1) && 
        (M2_State.Domain_Type != BOUNDARY_GAM2) && 
        (M2_State.Domain_Type != BOUNDARY_GAM3) && 
        (M2_State.Domain_Type != GAM1_GAM2_GAM3)  ) {
        Set_Moments(&M2_State, Var_index, num_points, x_Lebed, y_Lebed, z_Lebed, &Mobius_Scale_Params);
        
        if (M2_State.Domain_Type == BOUNDARY_FREE_STREAMING) {
            if (M2_State.display) {
                if (M2_State.id_proc == M2_State.proc_display) {
                    cout << "....................... Free-streaming limit ......................" << endl;
                }
            }
            goto Continue;
        } else if ((M2_State.Domain_Type == BOUNDARY_GAM1_GAM2_EQ_0) || 
                   (M2_State.Domain_Type == BOUNDARY_GAM1_GAM3_EQ_0) || 
                   (M2_State.Domain_Type == BOUNDARY_GAM2_GAM3_EQ_0)   ) {
            
            if (M2_State.display) {
                if (M2_State.id_proc == M2_State.proc_display) {
                    cout << "....................... Triangle vertices ......................" << endl;
                }
            }
            goto Continue;
        } else {
            if (M2_State.display) { 
                if (M2_State.id_proc == M2_State.proc_display) {
                    cout << "No Domain Type: Fatal Error" << endl;   
                }
            }
            exit(0);
        }
    } else if (M2_State.flag_Taylor_Series_Expansion_Sphere) {
        if (M2_State.display) {
            if (M2_State.id_proc == M2_State.proc_display) {
                cout << "....................... Taylor Series near free-streaming limit ......................" << endl;
            }
        }
        goto Continue;
    } else if (M2_State.flag_Taylor_Series_Expansion_Triangle) {
        if (M2_State.display) {
            if (M2_State.id_proc == M2_State.proc_display) {
                cout << "....................... Taylor Series near Triangle Boundaries ......................" << endl;
            }
        }
        goto Continue;
    }
    
    for (int index_reg = 0; index_reg < max_regularization; index_reg++) {
        if (index_reg > 0) {
            if (M2_State.display) { 
                if (M2_State.id_proc == M2_State.proc_display) {
                    printf(".................Moment Regularization with r_l = %Le..................\n", r_l[index_reg]);
                }   
            }
        }
        
        // In this case we are either inside the realizable space for angular moments up to second-order, or along one of the edges of the triangle P1 P2 P3
        for ( int id_refinement = 0; id_refinement <= M2_State.max_refinement_level; id_refinement++ ) {
            // Set moments
            Set_Moments(&M2_State, Var_index, num_points, x_Lebed, y_Lebed, z_Lebed, &Mobius_Scale_Params);
            
            // Setup isotropic moments for the purpose of regularization
            Set_Iso_Moments(M2_State);
            
            // Regularize the given set of angular moments
            Regularize_Moment(r_l[index_reg], &M2_State);
            
            // Setup number of quadrature points for the purpose of numerical integration
            // over the sphere
            M2_State.Order_Quad_mu = 5;
            // Setup quadrature abscissas and weight for the purpose of numerical integration
            // over the sphere
            M2_State.Setup_Quadrature(id_refinement); // This function should be called after Set_Moments
            
            // Set initial guess and tolerance for the purpose of the numerical optimization
            Set_Initial_Guess_And_Tol(M2_State);
            for (int i = 0; i < M2_State.NVARS; i++) {
                x[i] = M2_State.x[i];
                tol_x[i] = M2_State.tol_x[i];
            }
            
            // create the non-linear optimization algorithm
            opt = nlopt_create(NLOPT_LD_SLSQP, M2_State.NVARS); /* algorithm and dimensionality */
            nlopt_set_min_objective_orthog(opt, myfunc, Q_func_Gram_Schmidt, &M2_State);
            
            nlopt_set_maxeval(opt, max_iters_nlopt);
            
            nlopt_set_xtol_abs(opt, tol_x);
            nlopt_set_ftol_abs(opt, tol_grad);
            
            my_constraint_data data[1];
            add_constraints(opt, data, M2_State);
            
            if (nlopt_optimize_orthog(opt, x, M2_State.I0, &minf, &min_grad, Sk_final) < 0) {
                if (M2_State.display) { 
//                     if (M2_State.id_proc == M2_State.proc_display) {
                        printf("....................nlopt failed!..................... min_grad = %0.10g \n", min_grad);
                        cout << "id_proc = " << M2_State.id_proc << endl;
                        cout << "Failure with set of moments "  << "ratio_I0 = " << M2_State.ratio_I0 << "   " << "I0 = " << M2_State.I0 << "   " << "N1_1 = " << M2_State.N1_1 << "    " << "N1_2 = " << M2_State.N1_2 << "   " << "N1_3 = " << M2_State.N1_3 << "     " << "gamma_1 = " << M2_State.gamma_1 << "    " << "gamma_2 = " << M2_State.gamma_2 << "    " << "id = " << M2_State.id_proc << "     " << "grad_norm = " << min_grad << endl;
//                     }
                }
                
                if (index_reg == max_regularization-1 && id_refinement == M2_State.max_refinement_level) {
                    if (M2_State.display) { 
                        if (M2_State.id_proc == M2_State.proc_display) {
                            printf("....................nlopt failed: Exiting!.....................\n");
                            cout << "Failure with set of moments "  << "I0 = " << M2_State.I0 << "   " << "N1_1 = " << M2_State.N1_1 << "    " << "N1_2 = " << M2_State.N1_2 << "   " << "N1_3 = " << M2_State.N1_3 << "     " << "gamma_1 = " << M2_State.gamma_1 << "    " << "gamma_2 = " << M2_State.gamma_2 << "    " << "id = " << M2_State.id_proc << "     " << "grad_norm = " << min_grad << endl;
                        }
                    }
                    exit(0);
                    goto Exiting;
                }
            } else {
                if (min_grad != min_grad ) {
                    printf("Gradient in nan, g = %0.10g................\n", min_grad);
                    
                    if (index_reg == max_regularization-1 && id_refinement == M2_State.max_refinement_level) {
                        if (M2_State.display) { 
                            if (M2_State.id_proc == M2_State.proc_display) {
                                printf("....................nlopt failed: Exiting!.....................\n");
                            }
                        }
                        exit(0);
                        goto Exiting;
                    }
                } else if (min_grad > tol_grad ) {
                    if (M2_State.display) { 
                        if (M2_State.id_proc == M2_State.proc_display) {
                            printf("***********Tolerance on gradient not satisfied g(");
                            for (int index_vars = 0; index_vars < M2_State.NVARS; index_vars++) {
                                if (index_vars < M2_State.NVARS - 1) {
                                    printf("%Lf,", x[index_vars]);
                                } else {
                                    printf("%Lf", x[index_vars]);
                                }
                            }
                            printf(") = %0.10Le................\n", min_grad);
                        }
                    }
                    
                    if (index_reg == max_regularization-1 && id_refinement == M2_State.max_refinement_level) {
                        if (M2_State.display) { 
                            if (M2_State.id_proc == M2_State.proc_display) {
                                printf("....................nlopt failed: Exiting!.....................\n");
                            }
                        }
                        exit(0);
                        goto Exiting;
                    }
                } else {
                    if (M2_State.display) {
                        // cout << "****** Checking gradient !!!!!!!!!!!!!!!!!!!! *******" << endl;
                        // Check_Gradient(x, Sk_final, M2_State, min_grad);
                        
                        if (M2_State.id_proc == M2_State.proc_display) {
                            printf("***********found minimum at f(");
                            for (int index_vars = 0; index_vars < M2_State.NVARS; index_vars++) {
                                if (index_vars < M2_State.NVARS - 1) {
                                    printf("%Lf,", x[index_vars]);
                                } else {
                                    printf("%Lf", x[index_vars]);
                                }
                            }
                            printf(") = %0.10Lf***********\n ", minf);
                            printf("***********Final gradient is min_grad = %0.10Le *********** \n\n", min_grad);
                        }
                    }
                    
                    goto Continue;
                }   
            }
            if (opt != NULL) {
                nlopt_destroy(opt);
                opt = NULL;
            }
        }
    }
    Continue:;
    
    // Now compute third-order closing fluxes or partial angular moments based on the 
    // regime of radiation transfer encountered
    switch (M2_State.Max_Ent_Solution_Type) {
        case CLOSING_FLUX:
            Compute_Third_Order_Closing_Fluxes(rec_N3_local, id_count, Var_index, num_points, M2_State, x, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, Sk_final);
            if (!M2_State.flag_Taylor_Series_Expansion_Sphere &&
                !M2_State.flag_Taylor_Series_Expansion_Triangle) {
                // We only need to compute derivatives separately in the case where we are not
                // performing a taylor series expansion
                // This is because flag_Taylor_Series_Expansion = true ==> flag_Taylor_Series_Expansion
                Compute_Third_Order_Closing_Fluxes_Derivatives(rec_N3_local, id_count, Var_index, num_points, M2_State, x, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, Sk_final);
            } else {
                cout << "exit for now!!!!!!!!" << endl;
                exit(0);
            }
            break;
        case PARTIAL_MOMENTS:
//             Compute_Partial_Moments(rec_Partial_Moments_local, id_count, Var_index, num_points, M2_State, x, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, Sk_final);
//             
//             if (!M2_State.flag_Taylor_Series_Expansion_Sphere &&
//                 !M2_State.flag_Taylor_Series_Expansion_Triangle) {
//                 // We only need to compute derivatives separately in the case where we are not
//                 // performing a taylor series expansion
//                 Compute_Partial_Moments_Derivatives(rec_Partial_Moments_local, id_count, Var_index, num_points, M2_State, x, Mobius_Scale_Params, x_Lebed, y_Lebed, z_Lebed, Sk_final);
//             }
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    }
    
    Exiting:;
    
    M2_State.Deallocate_Quad();
    
    if (x != NULL) {
        delete[] x; x = NULL;
    }
    if (tol_x != NULL) {
        delete[] tol_x; tol_x = NULL;
    }
    if (Sk_final != NULL) {
        delete[] Sk_final; Sk_final = NULL;
    }
    
    if (opt != NULL) {
        nlopt_destroy(opt);
        opt = NULL;
    }
}

void Triangle_to_Square_Mapping_Nodes(long double &zeta, long double &eta, const int &i, const int &j, const int &m, const int &Node_Dist_gam1) {
    long double v_x, v_y, v_z;
    int k;
    k = m - 1 - i - j;
    
    if (Node_Dist_gam1 == UNIFORM_DISTRIBUTION) {
        v_x = Uniform_Distribution_No_Endpoint(i, m - 1, 0.0, 1.0);
        v_y = Uniform_Distribution_No_Endpoint(j, m - 1, 0.0, 1.0);
        v_z = Uniform_Distribution_No_Endpoint(k, m - 1, 0.0, 1.0);
        
        // v_x = Uniform_Distribution(i, m - 1, 0.0, 1.0);
        // v_y = Uniform_Distribution(j, m - 1, 0.0, 1.0);
        // v_z = Uniform_Distribution(k, m - 1, 0.0, 1.0);
    } else {
        v_x = zeros_shifted(i, m, 0.0, 1.0, Node_Dist_gam1);
        v_y = zeros_shifted(j, m, 0.0, 1.0, Node_Dist_gam1);
        v_z = zeros_shifted(k, m, 0.0, 1.0, Node_Dist_gam1);
        
        // cout << "v_x = " << v_x << "  " << "v_y = " << v_y << "  " << "v_z = " << v_z << endl;
    }
    
    zeta = (1.0/3.0)*(1.0 + 2.0*v_x - v_y - v_z);
    eta = (1.0/3.0)*(1.0 + 2.0*v_y - v_x - v_z);
    
    // cout << "zeta = " << zeta << "  " << "eta = " << eta << endl;
}

int Check_Realizability(M2_State_Param *M2_State) {
    int flag_Realizability = 1;
    
    M2_State->flag_Taylor_Series_Expansion_Sphere = false;
    M2_State->flag_Taylor_Series_Expansion_Triangle = false;
    
    if ((1.0 - M2_State->N1) >= -TOLER                          && 
        M2_State->gamma_1 >= -TOLER                             &&
        M2_State->gamma_2 >= -TOLER                             && 
        (1.0 - M2_State->gamma_1 - M2_State->gamma_2) >= -TOLER) {
        flag_Realizability = 0;
    }
    
    if (fabs(M2_State->N1 - 1.0) < 1.0e-6) {
        // Then do nothing because we are in free-streaming in this case
    } else if (M2_State->N1 > 1.0 - finite_diff_h_N1) {
//         if (M2_State->flag_finite_diff_Sphere) {
//             cout << "Inconsistency in finite differencing for N1 !!!!!!!!!!!!!" << endl;
//             cout << "N1 = " << M2_State->N1 << "  " << "N1_1 = " << M2_State->N1_1 << "  " << "N1_2 = " << M2_State->N1_2 << "  " << "N1_3 = " << M2_State->N1_3 << endl;
//             exit(0);
//         } else {
//             M2_State->flag_Taylor_Series_Expansion_Sphere = true;
//         }
    }
    
//     if (M2_State->Domain_Type == GAM1_GAM2_GAM3) {
//         if (M2_State->gamma_1 < finite_diff_h_gam                            ||
//             M2_State->gamma_2 < finite_diff_h_gam                            ||
//             (1.0 - M2_State->gamma_1 - M2_State->gamma_2) < finite_diff_h_gam) {
//             
//             if (M2_State->flag_finite_diff_Triangle) {
//                 cout << "Inconsistency in finite differencing for triangle in GAM1_GAM2_GAM3 !!!!!!!!!!!!!" << endl;
//                 cout << "gamma_1 = " << M2_State->gamma_1 << "  " << "gamma_2 = " << M2_State->gamma_2 << endl;
//                 exit(0);
//             } else {
//                 M2_State->flag_Taylor_Series_Expansion_Triangle = true;
//             }
//         }
//     } else if (M2_State->Domain_Type == BOUNDARY_GAM1) {
//         if (M2_State->gamma_2 < finite_diff_h_gam       ||
//             (1.0 - M2_State->gamma_2) < finite_diff_h_gam) {
//             
//             if (M2_State->flag_finite_diff_Triangle) {
//                 cout << "Inconsistency in finite differencing for triangle in BOUNDARY_GAM1 !!!!!!!!!!!!!" << endl;
//                 cout << "gamma_1 = " << M2_State->gamma_1 << "  " << "gamma_2 = " << M2_State->gamma_2 << endl;
//                 exit(0);
//             } else {
//                 M2_State->flag_Taylor_Series_Expansion_Triangle = true;
//             }
//         }
//     } else if (M2_State->Domain_Type == BOUNDARY_GAM2) {
//         if (M2_State->gamma_1 < finite_diff_h_gam       ||
//             (1.0 - M2_State->gamma_1) < finite_diff_h_gam) {
//             
//             if (M2_State->flag_finite_diff_Triangle) {
//                 cout << "Inconsistency in finite differencing for triangle in BOUNDARY_GAM2 !!!!!!!!!!!!!" << endl;
//                 cout << "gamma_1 = " << M2_State->gamma_1 << "  " << "gamma_2 = " << M2_State->gamma_2 << endl;
//                 exit(0);
//             } else {
//                 M2_State->flag_Taylor_Series_Expansion_Triangle = true;
//             }
//         }
//     } else if (M2_State->Domain_Type == BOUNDARY_GAM3) {
//         if (M2_State->gamma_1 < finite_diff_h_gam       ||
//             (1.0 - M2_State->gamma_1) < finite_diff_h_gam) {
//             
//             if (M2_State->flag_finite_diff_Triangle) {
//                 cout << "Inconsistency in finite differencing for triangle in BOUNDARY_GAM3 !!!!!!!!!!!!!" << endl;
//                 cout << "gamma_1 = " << M2_State->gamma_1 << "  " << "gamma_2 = " << M2_State->gamma_2 << endl;
//                 exit(0);
//             } else {
//                 M2_State->flag_Taylor_Series_Expansion_Triangle = true;
//             }
//         }
//     }
    
    return flag_Realizability;
}

long double generate_Iso_Moments(const int &index) {
    long double Mom_Iso;
    
    switch (index) {
        case 0:
            Mom_Iso = 1.0;
            break;
        case 4:
        case 7:
            Mom_Iso = 1.0/3.0;
            break;
        default:
            Mom_Iso = 0.0;
            break;
    };
    return Mom_Iso;
}

void Set_Iso_Moments(M2_State_Param &M2_State) {
    for (int i = 0; i < M2_State.NVARS; i++) {
        M2_State.MOM_VEC_iso[i] = M2_State.MOM_VEC[0] * generate_Iso_Moments(M2_State.Index[i]);
    }
}

long double define_Moment_Set(const long double &E, const long double &fx_test, const long double &fy_test, const long double &fz_test, const long double &pxx_test, const long double &pyy_test, const int &index) {
    long double Moment;
    switch (index) {
        case 1:
            Moment = fx_test * E;
            break;
        case 2:
            Moment = fy_test * E;
            break;
        case 3:
            Moment = fz_test * E;
            break;
        case 4:
            Moment = pxx_test * E;
            break;
        case 5:
            Moment = fx_test * fy_test * E;
            break;
        case 6:
            Moment = fx_test * fz_test * E;
            break;
        case 7:
            Moment = pyy_test * E; 
            break;
        case 8:
            Moment = fy_test * fz_test * E;
            break;
        default:
            cout << "Invalid value for index in define_Moment_Set !!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
    return Moment;
}

long double define_Moment_Set_Not_RT(const long double &E, const long double &fx_test, const long double &fy_test, const long double &fz_test, const long double &pxx_test, const long double &pxy_test, const long double &pyy_test, const int &index) {
    long double Moment;
    switch (index) {
        case 1:
            Moment = fx_test * E;
            break;
        case 2:
            Moment = fy_test * E;
            break;
        case 3:
            Moment = fz_test * E;
            break;
        case 4:
            Moment = pxx_test * E;
            break;
        case 5:
            Moment = pxy_test * E;
            break;
        case 6:
            Moment = fx_test * fz_test * E;
            break;
        case 7:
            Moment = pyy_test * E; 
            break;
        case 8:
            Moment = fy_test * fz_test * E;
            break;
        default:
            cout << "Invalid value for index in define_Moment_Set !!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
    return Moment;
}

void Check_Domain_Type(int &Domain_Type, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, M2_State_Param &M2_State) {
    long double norm_f;
    long double gamma_1, gamma_2;
    
    if (M2_State.flag_finite_diff_Sphere) {
        switch (M2_State.finite_diff_domain_Sphere) {
            case RADIUS_SPHERE_FINITE_DIFFERENCE:
                norm_f = M2_State.norm_f_knot + M2_State.x_finite_diff_norm_f;
                // cout << "norm_f = " << norm_f << "   " << "norm_f_knot = " << M2_State.norm_f_knot << "   " << "x_finite_diff_norm_f = " << M2_State.x_finite_diff_norm_f << endl; 
                break;
            default:
                cout << "Type of finite difference on sphere not specified" << endl;
                exit(0);
                break;
        };
    } else {
        if (M2_State.Node_Dist_f == UNIFORM_DISTRIBUTION) {
            norm_f = Uniform_Distribution(Var_index->f, num_points->f - 1, 0.0, 1.0);
        } else {
            norm_f = zeros_shifted(num_points->f - 1 + Var_index->f, 2*(num_points->f - 1) + 1, -1.0, 1.0, M2_State.Node_Dist_f);
        }
    }
    M2_State.N1 = norm_f;
    
    if (M2_State.flag_finite_diff_Triangle) {
        gamma_1 = M2_State.gam1_knot + M2_State.x_finite_diff_gam1;
        gamma_2 = M2_State.gam2_knot + M2_State.x_finite_diff_gam2;
        switch (M2_State.finite_diff_domain_Triangle) {
            case EDGE_TRIANGLE_FINITE_DIFFERENCE:
                if (M2_State.face_triangle == FACE_A) {
                    if (gamma_1 > 1.0e-8) {
                        cout << "Incorrect value for gamma_1 along FACE_A, gamma_1 = " << gamma_1 << endl;
                        exit(0);
                    }
                } else if (M2_State.face_triangle == FACE_B) {
                    if (gamma_2 > 1.0e-8) {
                        cout << "Incorrect value for gamma_2 along FACE_B, gamma_2 = " << gamma_2 << endl;
                        exit(0);
                    }
                } else if (M2_State.face_triangle == FACE_C) {
                    if (fabs(1.0 - gamma_1 - gamma_2) > 1.0e-8) {
                        cout << "Incorrect value for gamma_1 and gamma_2 along FACE_C, gamma_1 + gamma_2 = " << gamma_1 + gamma_2 << endl;
                        exit(0);
                    }
                } else {
                    cout << "Triangle Face for finite difference not specified, face_triangle = " << M2_State.face_triangle << endl;
                    exit(0);
                }
                break;
            case MEDIAN_TRIANGLE_FINITE_DIFFERENCE:
                if (fabs(M2_State.gam1_knot - 1.0/3.0) < 1.0e-8 || fabs(M2_State.gam2_knot - 1.0/3.0) < 1.0e-8) {
                    
                } else if (fabs(M2_State.a_median * gamma_1 + M2_State.b_median - gamma_2) > 1.0e-8) {
                    cout << "Incorrect value for gamma_1 and gamma_2 along median line of interest, gamma_1 = " << gamma_1 << "  " << "gamma_2 = " << gamma_2 << "  " << "Predicted gamma_2 = " << M2_State.a_median * gamma_1 + M2_State.b_median << endl;
                    exit(0);
                }
                break;
            default:
                cout << "Finite difference domain over triangle not specified" << endl;
                exit(0);
                break;
        }
        
        if (fabs(gamma_1) < 1.0e-8) {
            gamma_1 = 0.0;
        }
        
        if (fabs(gamma_2) < 1.0e-8) {
            gamma_2 = 0.0;
        }
        
        if (gamma_1 < 0.0 || gamma_2 < 0.0) {
            cout << "One of the gammas is negative!!!! gamma_1 = " << gamma_1 << "  " << "gamma_2 = " << gamma_2 << endl;
            exit(0);
        }
    } else {
        if (M2_State.Triangle_Domain == FULL_TRIANGLE) {
            Triangle_to_Square_Mapping_Nodes(gamma_1, gamma_2, Var_index->gam1, Var_index->gam2, num_points->gam1, M2_State.Node_Dist_gam1);
        } else if (M2_State.Triangle_Domain == GAM1_EQ_GAM2_TRIANGLE) {
            if (M2_State.Node_Dist_gam1 == UNIFORM_DISTRIBUTION) {
                gamma_1 = Uniform_Distribution(Var_index->gam1, num_points->gam1 - 1, 0.0, 1.0/2.0);
            } else {
                gamma_1 = zeros_shifted(Var_index->gam1, num_points->gam1, 0.0, 1.0/2.0, M2_State.Node_Dist_gam1);
            }
            gamma_2 = gamma_1;
        } else if (M2_State.Triangle_Domain == GAM1_EQ_GAM2_EQ_0_5_TRIANGLE) {
            gamma_1 = 0.5;
            gamma_2 = gamma_1;
        } else if (M2_State.Triangle_Domain == GAM1_EQ_GAM2_EQ_GAM3_TRIANGLE) {
            gamma_1 = 1.0/3.0;
            gamma_2 = gamma_1;
        } else if (M2_State.Triangle_Domain == GAM3_EQ_0_TRIANGLE) {
            if (M2_State.Node_Dist_gam1 == UNIFORM_DISTRIBUTION) {
                gamma_1 = Uniform_Distribution(Var_index->gam1, num_points->gam1 - 1, 0.0, 1.0);
            } else {
                gamma_1 = zeros_shifted(Var_index->gam1, num_points->gam1, 0.0, 1.0, M2_State.Node_Dist_gam1);
            }
            gamma_2 = 1.0 - gamma_1;
        } else if (M2_State.Triangle_Domain == GAM2_EQ_0_TRIANGLE) {
            if (M2_State.Node_Dist_gam1 == UNIFORM_DISTRIBUTION) {
                gamma_1 = Uniform_Distribution(Var_index->gam1, num_points->gam1 - 1, 0.0, 1.0);
            } else {
                gamma_1 = zeros_shifted(Var_index->gam1, num_points->gam1, 0.0, 1.0, M2_State.Node_Dist_gam1);
            }
            gamma_2 = 0.0;
        } else {
            if (M2_State.display) {
                if (M2_State.id_proc == M2_State.proc_display) {
                    cout << " **************************** Triangle_Domain not Defined *************************" << endl;
                }
            }
            exit(0);
        }
    }
    
    M2_State.gamma_1 = gamma_1;
    M2_State.gamma_2 = gamma_2;
    M2_State.gamma_3 = 1.0 - gamma_1 - gamma_2;
    
    if (M2_State.Rec_Moms_Test.flag_Moments_Test) {
        M2_State.Rec_Moms_Test.Override_Moments_Tests(M2_State.MOM_VEC[0], M2_State.ratio_I0, norm_f, M2_State.x_SH, M2_State.y_SH, M2_State.z_SH, M2_State.gamma_1, M2_State.gamma_2, M2_State.Problem_Type);
    }
    
//     cout << "Domain Type Info =========>    " << "gamma_1 = " << gamma_1 << "    " << "gamma_2 = " << gamma_2 << endl;
    
    if (norm_f >= 1.0 - TOLER) {
        Domain_Type = BOUNDARY_FREE_STREAMING;
    } else {
        if (gamma_1 < TOLER) {
            if (gamma_2 < TOLER) {
                if (M2_State.display) {
                    if (M2_State.id_proc == M2_State.proc_display) {
                        cout << "Boundary Point" << endl;
                    }
                }
                Domain_Type = BOUNDARY_GAM1_GAM2_EQ_0;
            } else if (fabs(1.0 - gamma_1 - gamma_2) < TOLER) {
                if (M2_State.display) {
                    if (M2_State.id_proc == M2_State.proc_display) {
                        cout << "Boundary Point" << endl;
                    }
                }
                Domain_Type = BOUNDARY_GAM1_GAM3_EQ_0;
            } else {
                if (M2_State.display) {
                    if (M2_State.id_proc == M2_State.proc_display) {
                        cout << "BOUNDARY_GAM1" << endl;
                    }
                }
                Domain_Type = BOUNDARY_GAM1;
            }
        } else if (gamma_2 < TOLER) {
            if (gamma_1 > 1.0 - TOLER) {
                if (M2_State.display) {
                    if (M2_State.id_proc == M2_State.proc_display) {
                        cout << "Boundary Point" << endl;
                    }
                }
                Domain_Type = BOUNDARY_GAM2_GAM3_EQ_0;
            } else {
                if (M2_State.display) {
                    if (M2_State.id_proc == M2_State.proc_display) {
                        cout << "BOUNDARY_GAM2" << endl;
                    }
                }
                Domain_Type = BOUNDARY_GAM2;
            }
        } else if ( 1.0 - gamma_1 - gamma_2 < TOLER) {
            if (M2_State.display) {
                if (M2_State.id_proc == M2_State.proc_display) {
                    cout << "BOUNDARY_GAM3" << endl;
                }
            }
            Domain_Type = BOUNDARY_GAM3;
        } else {
            if (M2_State.display) {
                if (M2_State.id_proc == M2_State.proc_display) {
                    cout << "GAM1_GAM2_GAM3" << endl;
                }
            }
            Domain_Type = GAM1_GAM2_GAM3;
        }
    }
    
    M2_State.flag_Realizability = Check_Realizability(&M2_State);
}

void Set_Moments(M2_State_Param *M2_State, const MN_Var_Num_Points *Var_index, const MN_Var_Num_Points *num_points, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, Mobius_Scale_Parameters *Mobius_Scale_Params) {
    long double norm_f, norm_f_test_2, fx_test, fy_test, fz_test;
    long double Increment_f_test, Increment_Theta, Increment_Phi, Increment_gam1, Increment_gam2;
    long double pxx_test, pxy_test, pyy_test;
    long double gamma_1, gamma_2;
    long double Omega1, Omega2, Omega3;
    long double MOM_VEC_0;
    long double mu, Phi, Theta;
    
    if (M2_State->flag_finite_diff_Sphere) {
        switch (M2_State->finite_diff_domain_Sphere) {
            case RADIUS_SPHERE_FINITE_DIFFERENCE:
                norm_f = M2_State->norm_f_knot + M2_State->x_finite_diff_norm_f;
                break;
            default:
                cout << "Type of finite difference on sphere not specified" << endl;
                exit(0);
                break;
        };
    } else {
        if (M2_State->Node_Dist_f == UNIFORM_DISTRIBUTION) {
            norm_f = Uniform_Distribution(Var_index->f, num_points->f - 1, 0.0, 1.0);
        } else {
            norm_f = zeros_shifted(num_points->f - 1 + Var_index->f, 2*(num_points->f - 1) + 1, -1.0, 1.0, M2_State->Node_Dist_f);
        }
    }
    
//     if (fabs(norm_f) < 1.0e-8) {
//         norm_f = 0.1;
//     }
    
    if (M2_State->Node_Dist_Phi_Theta == UNIFORM_DISTRIBUTION) {
        Phi = Uniform_Distribution(Var_index->phi, num_points->phi - 1, 0.0, 2.0*PI);
        Theta = Uniform_Distribution(Var_index->theta, num_points->theta - 1, 0.0, PI);
        
        // Get rid of possible roundoff errors
        Phi = roundn(Phi, 12);
        Theta = roundn(Theta, 12);
        
        // Theta = PI/2.0;
        
        Omega1 = cos(Theta);
        Omega2 = sin(Theta)*cos(Phi);
        Omega3 = sin(Theta)*sin(Phi);
    } else if (M2_State->Node_Dist_Phi_Theta == SPHERICAL_HARMONIC_DISTRIBUTION) {
        Omega1 = x_Lebed[Var_index->phi*num_points->theta + Var_index->theta];
        Omega2 = y_Lebed[Var_index->phi*num_points->theta + Var_index->theta];
        Omega3 = z_Lebed[Var_index->phi*num_points->theta + Var_index->theta];
    } else if (M2_State->Node_Dist_Phi_Theta == FLUX_ONE_D) {
        Omega1 = 1.0;
        Omega2 = 0.0;
        Omega3 = 0.0;
    } else {
        cout << "Node Distribution for Phi and Theta not specified !!!!!!!!!!!!!!!!!!!" << endl;
        exit(0);
    }
    
    if (M2_State->flag_finite_diff_Triangle) {
        gamma_1 = M2_State->gam1_knot + M2_State->x_finite_diff_gam1;
        gamma_2 = M2_State->gam2_knot + M2_State->x_finite_diff_gam2;
        
        switch (M2_State->finite_diff_domain_Triangle) {
            case EDGE_TRIANGLE_FINITE_DIFFERENCE:
                if (M2_State->face_triangle == FACE_A) {
                    if (gamma_1 > 1.0e-8) {
                        cout << "Incorrect value for gamma_1 along FACE_A, gamma_1 = " << gamma_1 << endl;
                        exit(0);
                    }
                } else if (M2_State->face_triangle == FACE_B) {
                    if (gamma_2 > 1.0e-8) {
                        cout << "Incorrect value for gamma_2 along FACE_B, gamma_2 = " << gamma_2 << endl;
                        exit(0);
                    }
                } else if (M2_State->face_triangle == FACE_C) {
                    if (fabs(1.0 - gamma_1 - gamma_2) > 1.0e-8) {
                        cout << "Incorrect value for gamma_1 and gamma_2 along FACE_C, gamma_1 + gamma_2 = " << gamma_1 + gamma_2 << endl;
                        exit(0);
                    }
                } else {
                    cout << "Triangle Face for finite difference not specified, face_triangle = " << M2_State->face_triangle << endl;
                    exit(0);
                }
                break;
            case MEDIAN_TRIANGLE_FINITE_DIFFERENCE:
                if (fabs(M2_State->gam1_knot - 1.0/3.0) < 1.0e-8 || fabs(M2_State->gam2_knot - 1.0/3.0) < 1.0e-8) {
                    
                } else if (fabs(M2_State->a_median * gamma_1 + M2_State->b_median - gamma_2) > 1.0e-8) {
                    cout << "Incorrect value for gamma_1 and gamma_2 along median line of interest, gamma_1 = " << gamma_1 << "  " << "gamma_2 = " << gamma_2 << "  " << "Predicted gamma_2 = " << M2_State->a_median * gamma_1 + M2_State->b_median << endl;
                    exit(0);
                }
                break;
            default:
                cout << "Finite difference domain over triangle not specified" << endl;
                exit(0);
                break;
        }
    } else {
        if (M2_State->Triangle_Domain == FULL_TRIANGLE) {
            Triangle_to_Square_Mapping_Nodes(gamma_1, gamma_2, Var_index->gam1, Var_index->gam2, num_points->gam1, M2_State->Node_Dist_gam1);
        } else if (M2_State->Triangle_Domain == GAM1_EQ_GAM2_TRIANGLE) {
            if (M2_State->Node_Dist_gam1 == UNIFORM_DISTRIBUTION) {
                gamma_1 = Uniform_Distribution(Var_index->gam1, num_points->gam1 - 1, 0.0, 1.0/2.0);
            } else {
                gamma_1 = zeros_shifted(Var_index->gam1, num_points->gam1, 0.0, 1.0/2.0, M2_State->Node_Dist_gam1);
            }
            gamma_2 = gamma_1;
        } else if (M2_State->Triangle_Domain == GAM1_EQ_GAM2_EQ_0_5_TRIANGLE) {
            mu = zeros_shifted(num_points->theta - 1 + Var_index->theta, 2*(num_points->theta - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            fx_test = mu * norm_f;
            fy_test = sqrt(1.0 - pow(mu, 2))* norm_f;
            fz_test = 0.0;
            
            gamma_1 = 0.5;
            gamma_2 = gamma_1;
        } else if (M2_State->Triangle_Domain == GAM1_EQ_GAM2_EQ_GAM3_TRIANGLE) {
            gamma_1 = 1.0/3.0;
            gamma_2 = gamma_1;
        } else if (M2_State->Triangle_Domain == GAM3_EQ_0_TRIANGLE) {
            mu = zeros_shifted(num_points->theta - 1 + Var_index->theta, 2*(num_points->theta - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            fx_test = mu * norm_f;
            fy_test = sqrt(1.0 - pow(mu, 2))* norm_f;
            fz_test = 0.0;
            
            if (M2_State->Node_Dist_gam1 == UNIFORM_DISTRIBUTION) {
                gamma_1 = Uniform_Distribution(Var_index->gam1, num_points->gam1 - 1, 0.0, 1.0);
            } else {
                gamma_1 = zeros_shifted(Var_index->gam1, num_points->gam1, 0.0, 1.0, M2_State->Node_Dist_gam1);
            }
            gamma_2 = 1.0 - gamma_1;
        } else if (M2_State->Triangle_Domain == GAM2_EQ_0_TRIANGLE) {
            mu = zeros_shifted(num_points->theta - 1 + Var_index->theta, 2*(num_points->theta - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            fx_test = mu * norm_f;
            fy_test = 0.0;
            fz_test = sqrt(1.0 - pow(mu, 2))* norm_f;
            
            if (M2_State->Node_Dist_gam1 == UNIFORM_DISTRIBUTION) {
                gamma_1 = Uniform_Distribution(Var_index->gam1, num_points->gam1 - 1, 0.0, 1.0);
            } else {
                gamma_1 = zeros_shifted(Var_index->gam1, num_points->gam1, 0.0, 1.0, M2_State->Node_Dist_gam1);
            }
            gamma_2 = 0.0;
        } else {
            if (M2_State->display) {
                if (M2_State->id_proc == M2_State->proc_display) {
                    cout << " **************************** Triangle_Domain not Defined *************************" << endl;
                }
            }
            exit(0);
        }
    }
    
    bool flag_finite_diff_Energy_Spectrum = false;
    
    if (flag_finite_diff_Energy_Spectrum) {
        switch (M2_State->finite_diff_domain_Spectrum) {
            case SPECTRUM_ENERGY_FINITE_DIFFERENCE:
                M2_State->ratio_I0 = M2_State->ratio_I0_knot + M2_State->x_finite_diff_ratio_I0;
                break;
            default:
                cout << "Type of finite difference over the spectrum not specified" << endl;
                exit(0);
                break;
        };
    } else {
        if (M2_State->Node_Dist_E == DISCRETE_SET_ENERGY) {
            if (Var_index->E == 0) {
                M2_State->ratio_I0 = -1.0;
            } else if (Var_index->E == num_points->E - 1) {
                M2_State->ratio_I0 = 1.0;
            } else {
                M2_State->ratio_I0 = 0.0;
            }
            M2_State->MOM_VEC[0] = E_vals[Var_index->E];
        } else {
            if (M2_State->Node_Dist_E == UNIFORM_DISTRIBUTION) {
                M2_State->ratio_I0 = Uniform_Distribution(Var_index->E, num_points->E - 1, -1.0, 1.0);
            } else {
                M2_State->ratio_I0 = zeros_shifted(Var_index->E, num_points->E, -1.0, 1.0, M2_State->Node_Dist_E);
            }
            
            // Get rid of possible roundoff errors
            M2_State->ratio_I0 = roundn(M2_State->ratio_I0, 12);
            
            M2_State->MOM_VEC[0] = Inverse_Mobius_Transformation(M2_State->ratio_I0, Mobius_Scale_Params->Evaluate_Length_Scale(norm_f, Omega1, Omega2, Omega3, gamma_1, gamma_2));
            
            
//             if (M2_State->display) {
//                 if (M2_State->id_proc == M2_State->proc_display) {
//                     cout << "Length_Scale = " << Mobius_Scale_Params->Evaluate_Length_Scale(norm_f, Omega1, Omega2, Omega3, gamma_1, gamma_2) << endl;
//                     
//                     cout << "Coefficients_Mobius_Scale_Fit[0] = " << Mobius_Scale_Params->Coefficients_Mobius_Scale_Fit[0] << endl;
//                 }
//             }
        }
    }
    
    if (num_points->E == 1 && M2_State->Problem_Type == NON_GRAY && !M2_State->Rec_Moms_Test.flag_Moments_Test) {
        switch (num_points->Maximum_Entropy_Solution_Regime) {
            case HYPERBOLIC_LIMIT:
                M2_State->ratio_I0 = -1.0;
                break;
            case LOGARITHMIC_LIMIT:
                M2_State->ratio_I0 = 1.0;
                break;
            default:
                cout << "Invalid value for Maximum_Entropy_Solution_Regime" << endl;
                exit(0);
                break;
        }
    }
    
//     if (M2_State->display) {
//         if (M2_State->id_proc == M2_State->proc_display) {
//             cout << "ratio_I0 = " << M2_State->ratio_I0 << "   " << "Length_Scale = " << Mobius_Scale_Params->Evaluate_Length_Scale(norm_f, Omega1, Omega2, Omega3, gamma_1, gamma_2) << "   " << "MOM_VEC[0] = " << M2_State->MOM_VEC[0] << endl;
//         }
//     }
    
    if (M2_State->Rec_Moms_Test.flag_Moments_Test) {
        M2_State->Rec_Moms_Test.Override_Moments_Tests(M2_State->MOM_VEC[0], M2_State->ratio_I0, norm_f, Omega1, Omega2, Omega3, gamma_1, gamma_2, M2_State->Problem_Type);
    }
    
    // Get rid of possible roundoff errors
    norm_f = roundn(norm_f, 12);
    norm_f_test_2 = pow(norm_f,2);
    
    // Get rid of possible roundoff errors
    gamma_1 = roundn(gamma_1, 12);
    gamma_2 = roundn(gamma_2, 12);
    
    fx_test = Omega1*norm_f;
    fy_test = Omega2*norm_f;
    fz_test = Omega3*norm_f;
    
    pxx_test = pow(fx_test,2) + gamma_1*(1.0 - norm_f_test_2);
    
    pyy_test = pow(fy_test,2) + gamma_2*(1.0 - norm_f_test_2);
    
    Set_Regime(M2_State);
    
//     if (M2_State->Domain_Type == BOUNDARY_GAM1 || M2_State->Domain_Type == BOUNDARY_GAM2 || M2_State->Domain_Type == BOUNDARY_GAM3) {
//         M2_State->Regime = HYPERBOLIC_LIMIT; //LOGARITHMIC_LIMIT;
//     }
    
//     fx_test = 7.8631249355091637e-02;
//     fy_test = 1.6759814188776967e-01;
//     fz_test = 0.0;
//     pxx_test = 3.2843570390383575e-01;
//     pxy_test = 7.4016218727147030e-19;
//     pyy_test = 3.7104370900060629e-01;
    
    M2_State->I0 = M2_State->MOM_VEC[0];
    M2_State->N1 = norm_f;
    M2_State->N1_1 = fx_test;
    M2_State->N1_2 = fy_test;
    M2_State->N1_3 = fz_test;
    M2_State->x_SH = Omega1;
    M2_State->y_SH = Omega2;
    M2_State->z_SH = Omega3;
    M2_State->gamma_1 = gamma_1;
    M2_State->gamma_2 = gamma_2;
    M2_State->gamma_3 = 1.0 - gamma_1 - gamma_2;
    M2_State->N2_11 = pxx_test;
    M2_State->N2_12 = fx_test*fy_test;
    M2_State->N2_22 = pyy_test;
    
    // M2_State->N2_12 = pxy_test;
    
    M2_State->flag_Realizability = Check_Realizability(M2_State);
    
    for (int i = 1; i < M2_State->NVARS; i++) {
        // M2_State->MOM_VEC[i] = define_Moment_Set_Not_RT(M2_State->MOM_VEC[0], fx_test, fy_test, fz_test, pxx_test, pxy_test, pyy_test, M2_State->Index[i]);
        M2_State->MOM_VEC[i] = define_Moment_Set(M2_State->MOM_VEC[0], fx_test, fy_test, fz_test, pxx_test, pyy_test, M2_State->Index[i]);
    }
    
    if (M2_State->display) {
        if (M2_State->id_proc == M2_State->proc_display) {
            cout << "refinement_level = " << M2_State->refinement_level << "     " << "ratio_I0 = " << M2_State->ratio_I0 << "     " << "MOM_VEC[0] = " << M2_State->MOM_VEC[0] << "     " << "norm_f = " << norm_f << "     " << "fx_test = " << fx_test << "     " << "fy_test = " << fy_test << "     " << "fz_test = " << fz_test << "     "  << "gamma_1 = " << gamma_1 << "     "  << "gamma_2 = " << gamma_2 << endl;
            // cout << "flag_Moments_Test = " << M2_State->Rec_Moms_Test.flag_Moments_Test << endl;
            // cout << "Omega1 = " << Omega1 << "  " << "Omega2 = " << Omega2 << "  " << "Omega3 = " << Omega3 << endl;
        }
    }
}

 void add_constraints(nlopt_opt &opt, my_constraint_data *data, M2_State_Param &M2_State) {
     long double tol_constr = 0.0;
     // In tol_constr is negative then nlopt will return invalid arguments and not take
     // into account such constraints
     nlopt_add_inequality_constraint_orthog(opt, myconstraint, &M2_State, tol_constr);
 }

void generate_polynomials(long double* poly, const int &Id_angle, const M2_State_Param &M2_State) {
    for (int i = 0; i < M2_State.NVARS; i++) {
        poly[i] = generate_monomials_basis(M2_State.Omega1[Id_angle], M2_State.Omega2[Id_angle], M2_State.Omega3[Id_angle], M2_State.Index[i]);
    }
}

long double generate_monomials_basis(const long double &Omega1, const long double &Omega2, const long double &Omega3, const int &index) {
    long double poly;
    //     cout << "Omega1 = " << Omega1 << "     " << "Omega2 = " << Omega2 << "     " << "Omega3 = " << Omega3 << endl;
    
    switch (index) {
        case 0:
            poly = 1.0;
            break;
        case 1:
            poly = Omega1;
            break;
        case 2:
            poly = Omega2;
            break;
        case 3:
            poly = Omega3;
            break;
        case 4:
            poly = pow(Omega1,2);
            break;
        case 5:
            poly = Omega1*Omega2;
            break;
        case 6:
            poly = Omega1*Omega3;
            break;
        case 7:
            poly = pow(Omega2,2);
            break;
        case 8:
            poly = Omega2*Omega3; 
            break;
    };
    
    return poly;
}

long double generate_higher_order_monomials_basis(const long double &Omega1, const long double &Omega2, const long double &Omega3, const int &index) {
    long double poly;
    //     cout << "Omega1 = " << Omega1 << "     " << "Omega2 = " << Omega2 << "     " << "Omega3 = " << Omega3 << endl;
    
    switch (index) {
        case 0:
            poly = pow(Omega1,3);
            break;
        case 1:
            poly = Omega1*pow(Omega2,2);
            break;
        case 2:
            poly = Omega1*Omega2*Omega3;
            break;
        case 3:
            poly = pow(Omega1,2)*Omega2;
            break;
        case 4:
            poly = pow(Omega2,3);
            break;
        default:
            cout << "Invalid value for inded in generate_higher_order_monomials_basis !!!!!!!!!!!" << endl;
            exit(0);
            break;
    };
    
    return poly;
}

// long double generate_monomials_basis(const long double &Omega1, const long double &Omega2, const long double &Omega3, const int &index) {
//     long double poly;
// //     cout << "Omega1 = " << Omega1 << "     " << "Omega2 = " << Omega2 << "     " << "Omega3 = " << Omega3 << endl;
//     switch (index) {
//         case 0:
//             poly = 1.0;
//             break;
//         case 1:
//             poly = Omega1;
//             break;
//         case 2:
//             poly = Omega2;
//             break;
//         case 3:
//             poly = Omega3;
//             break;
//         case 4:
//             poly = pow(Omega1,2);
//             break;
//         case 5:
//             poly = Omega1*Omega2;
//             break;
//         case 6:
//             poly = Omega1*Omega3;
//             break;
//         case 7:
//             poly = pow(Omega2,2);
//             break;
//         case 8:
//             poly = Omega2*Omega3; 
//             break;
//     };
//     
//     return poly;
// }

int N3_M2_3D_DIRECTIONS(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     long double poly_xxx, poly_xxy, poly_xyy, poly_xyz, poly_yyy;
     long double Omega1, Omega2, Omega3;
     
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     
     coeff_1 = 1.0;
     
     long double *poly, *poly_Sk;
     poly = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     poly_Sk = new long double[M2_State->NVARS];
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
	 //cout << "x = " << M2_State->x[i] ;
     }
     //cout << endl;
     Omega1 = M2_State->Omega1[Id_angle];
     Omega2 = M2_State->Omega2[Id_angle];
     Omega3 = M2_State->Omega3[Id_angle];
     
//      Rotate_Frame_N3_Basis(Omega1, Omega2, Omega3, -M2_State->theta_N1, -M2_State->phi_N1);
     
     poly_xxx = pow(Omega1,3);
     poly_xxy = pow(Omega1,2)*Omega2;
     poly_xyy = Omega1*pow(Omega2,2);
     poly_xyz = Omega1*Omega2*Omega3;
     poly_yyy = pow(Omega2,3);
     
     //cout << "poly_xxx = " <<  << endl; 
     switch (M2_State->Regime) {
         case GRAY:
             MOMENTS[0] = coeff_1/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[1] = coeff_1*poly_xxx/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[2] = coeff_1*poly_xxy/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[3] = coeff_1*poly_xyy/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[4] = coeff_1*poly_xyz/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[5] = coeff_1*poly_yyy/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             break;
         case BOSE_EINSTEIN:
             MOMENTS[0] = coeff_1/(exp(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)) - 1.0);
             MOMENTS[1] = coeff_1*poly_xxx/(exp(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)) - 1.0);
             MOMENTS[2] = coeff_1*poly_xxy/(exp(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)) - 1.0);
             MOMENTS[3] = coeff_1*poly_xyy/(exp(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)) - 1.0);
             MOMENTS[4] = coeff_1*poly_xyz/(exp(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)) - 1.0);
             MOMENTS[5] = coeff_1*poly_yyy/(exp(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)) - 1.0);
             break;
         case LOGARITHMIC_LIMIT:
             MOMENTS[0] = coeff_1/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[1] = coeff_1*poly_xxx/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[2] = coeff_1*poly_xxy/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[3] = coeff_1*poly_xyy/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[4] = coeff_1*poly_xyz/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[5] = coeff_1*poly_yyy/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         case HYPERBOLIC_LIMIT:
             MOMENTS[0] = coeff_1*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[1] = coeff_1*poly_xxx*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[2] = coeff_1*poly_xxy*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[3] = coeff_1*poly_xyy*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[4] = coeff_1*poly_xyz*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[5] = coeff_1*poly_yyy*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         default:
             cout << "Invalid Regime Type for N3_M2_3D_DIRECTIONS !!!!!!!!!!!!!!!!" << endl;
             exit(0);
             break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
 void N3_M2_3D_RT(long double *N3, const long double *x, const long double *Sk, M2_State_Param &M2_State) {
    int NF_Chi;
    bool flag_realizable_closure;
    long double *Moments_Vals;
    NF_Chi = 5;
    Moments_Vals = new long double[NF_Chi+1];
    
    M2_State.set_x(x);
    M2_State.set_Sk(Sk);
    
    Lebedev_Quadrature_Func(Moments_Vals, NF_Chi+1, N3_M2_3D_DIRECTIONS, &M2_State);
    
    switch (M2_State.Domain_Type) {
        case GAM1_GAM2_GAM3:
            N3[0] = Moments_Vals[1]/Moments_Vals[0]; // ==> N3_111
            N3[1] = Moments_Vals[2]/Moments_Vals[0]; // ==> N3_112
            N3[2] = Moments_Vals[3]/Moments_Vals[0]; // ==> N3_122
            N3[3] = Moments_Vals[4]/Moments_Vals[0]; // ==> N3_123
            N3[4] = Moments_Vals[5]/Moments_Vals[0]; // ==> N3_222
            break;
        case BOUNDARY_GAM1:
            N3[0] = pow(M2_State.N1_1, 3);
            // N3[1] = ;
            N3[2] = M2_State.N1_1*(pow(M2_State.N1_2, 2) + M2_State.gamma_2*(1.0 - pow(M2_State.N1_1, 2) - pow(M2_State.N1_2, 2) - pow(M2_State.N1_3, 2)));
            N3[3] = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
            // N3[4] = ;
            break;
        case BOUNDARY_GAM2:
            N3[0] = Moments_Vals[1]/Moments_Vals[0];
            // N3[1] = ;
            N3[2] = M2_State.N1_1*pow(M2_State.N1_2, 2);
            N3[3] = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
            // N3[4] = ;
            break;
        case BOUNDARY_GAM3:
            N3[0] = Moments_Vals[1]/Moments_Vals[0];
            // N3[1] = ;
            N3[2] = Moments_Vals[2]/Moments_Vals[0];
            N3[3] = M2_State.N1_1*M2_State.N1_2*M2_State.N1_3;
            // N3[4] = ;
            break;
    };
    
//     flag_realizable_closure = Check_Realizability_Third_Order_Moment(M2_State.N1_1, M2_State.N2_11, M2_State.N2_22, N3[0], N3[1], N3[2]);
//     
//     if (!flag_realizable_closure) { 
//         cout << "N1_1 = " << M2_State.N1_1 << "  " << "N1_2 = " << M2_State.N1_2 << "  " << "N1_3 = " << M2_State.N1_3 << "  " << "gamma_1 = " << M2_State.gamma_1 << "  " << "gamma_2 = " << M2_State.gamma_2 << endl;
//         cout << "Non realizable third-order moment computed !!!!!!!!!!!!!" << endl;
//         exit(0);
//     }
    
//     if (fabs(Moments_Vals[0] - M2_State.I0) > 1.0e-6) {
    if (M2_State.display) {
        if (M2_State.id_proc == M2_State.proc_display) {
            cout << "E = " << Moments_Vals[0] << "    " << "N3_111 = " << N3[0] << "    " << "N3_112 = " << N3[1] << "    " << "N3_122 = " << N3[2] << "    " << "N3_123 = " << N3[3] << "    " << "N3_222 = " << N3[4] << "\n\n";
        }
    }
//         exit(0);
//     }
    
    delete[] Moments_Vals;
 }
 
 void Check_Gradient(const long double *x, const long double *Sk, M2_State_Param &M2_State, const long double &grad_SLSQP) {
    long double grad_norm;
    long double *gradients_vec, *MOM_VEC_new_basis;
    gradients_vec = new long double[M2_State.NVARS];
    MOM_VEC_new_basis = new long double[M2_State.NVARS];
    
    M2_State.set_x(x);
    M2_State.set_Sk(Sk);
    
    new_basis_Moments(M2_State.NVARS, MOM_VEC_new_basis, M2_State.MOM_VEC, Sk);
    
    Lebedev_Quadrature_Func(gradients_vec, M2_State.NVARS, gradient, &M2_State);
    
    
//     for (int i = 0; i < M2_State.NVARS; i++) {
//         for (int j = 0; j < M2_State.NVARS; j++) {
//             cout << "i = " << i << "    " << "j = " << j << "    "  << "Sk = " << Sk[i*M2_State.NVARS + j] << endl;
//         }
//     }
    
    for (int i = 0; i < M2_State.NVARS; i++) {
        gradients_vec[i] = gradients_vec[i] - MOM_VEC_new_basis[i];
        // cout << "i = " << i << "    "  << "grad_vals[i] = " << grad_vals[i] << "    " << "MOM_VEC_new_basis[i] = " << MOM_VEC_new_basis[i] << "    " << "grad[i] = " << grad[i] << endl;
    }
    
    
    grad_norm = 0.0;
    for (int i = 0; i < M2_State.NVARS; i++) {
        // printf("gradients_vec = %g \n", gradients_vec[i]);
        grad_norm += gradients_vec[i]*gradients_vec[i];
    }
    grad_norm = sqrt(grad_norm);
    
    grad_norm /= M2_State.I0;
    
    if (M2_State.display) {
        if (M2_State.id_proc == M2_State.proc_display) {
            cout << "grad_norm = " << grad_norm << "    " << "grad_norm_SLSQP = " << grad_SLSQP << "\n\n";
        }
    }
    
    delete[] gradients_vec;
 }
 
bool Check_Realizability_Third_Order_Moment(const long double &N1_1, const long double &N2_11, const long double &N2_22, const long double &N3_111, const long double &N3_122, const long double &N3_123) {
    bool flag_realizability = true;
    long double ub_N3_111, lb_N3_111;
    long double ub_N3_122, lb_N3_122;
    
    lb_N3_111 = -N2_11 + pow(N1_1 + N2_11, 2)/(1.0 + N1_1);
    ub_N3_111 = N2_11 - pow(N1_1 - N2_11, 2)/(1.0 - N1_1);
    
    lb_N3_122 = -N2_22;
    ub_N3_122 = N2_22;
    
    long double diff_lb_N3_111, diff_ub_N3_111;
    long double diff_lb_N3_122, diff_ub_N3_122;
    
    diff_lb_N3_111 = N3_111 - lb_N3_111;
    diff_ub_N3_111 = ub_N3_111 - N3_111;
    
    if (diff_lb_N3_111 < -1.0e-8 || diff_ub_N3_111 < -1.0e-8) {
        cout << "lb_N3_111 = " << lb_N3_111 << "  " << "N3_111 = " << N3_111 << "  " << "ub_N3_111 = " << ub_N3_111 << endl;
        cout << "diff_lb = " << N3_111 - lb_N3_111 << "  " << "diff_ub = " << ub_N3_111 - N3_111 << endl;
        flag_realizability = false;
    }
    
    diff_lb_N3_122 = N3_122 - lb_N3_122;
    diff_ub_N3_122 = ub_N3_122 - N3_122;
    
    if (diff_lb_N3_122 < -1.0e-8 || diff_ub_N3_122 < -1.0e-8) {
        cout << "lb_N3_122 = " << lb_N3_122 << "  " << "N3_122 = " << N3_122 << "  " << "ub_N3_122 = " << ub_N3_122 << endl;
        cout << "diff_lb = " << N3_122 - lb_N3_122 << "  " << "diff_ub = " << ub_N3_122 - N3_122 << endl;
        flag_realizability = false;
    }
    
    return flag_realizability;
}
 
void Regularize_Moment(const long double &r_l, M2_State_Param *M2_State) {
//     printf("r_l = %Le \n", r_l);
    // Applying regularization of moments by taking the convex combination
    // with the moments of the isotropic distribution
    for (int i = 1; i < M2_State->NVARS; i++) {
        M2_State->MOM_VEC[i] = (1.0 - r_l)*M2_State->MOM_VEC[i] + r_l * M2_State->MOM_VEC_iso[i];
    }
    
//     long double fx_test, fy_test, fz_test, pxx_test, pyy_test, gamma_1, gamma_2,norm_f_test_2;
//     fx_test = (1.0 - r_l)*M2_State->N1_1;
//     fy_test = (1.0 - r_l)*M2_State->N1_2;
//     fz_test = (1.0 - r_l)*M2_State->N1_3;
//     gamma_1 = (1.0 - r_l)*M2_State->gamma_1 + r_l*(1.0/3.0);
//     gamma_2 = (1.0 - r_l)*M2_State->gamma_2 + r_l*(1.0/3.0);
//     
//     norm_f_test_2 = pow(fx_test,2) + pow(fy_test,2) + pow(fz_test,2);
//     
//     pxx_test = pow(fx_test,2) + gamma_1*(1.0 - norm_f_test_2);
//     
//     pyy_test = pow(fy_test,2) + gamma_2*(1.0 - norm_f_test_2);
//     
//     for (int i = 1; i < M2_State->NVARS; i++) {
//         M2_State->MOM_VEC[i] = define_Moment_Set(M2_State->MOM_VEC[0], fx_test, fy_test, fz_test, pxx_test, pyy_test, M2_State->Index[i]);
//     }
}

long double sum_Lagrange(const long double *x, const long double *poly_Sk, const int &NVARS) {
    long double sum = 0.0;
    for (int i = 0; i < NVARS; i++) {
        sum += x[i]*poly_Sk[i];
    }
    
    return sum;
}

long double sum_Lag_Mom(const long double *x, const int &NVARS, const long double* MOMENTS) {
    long double sum;
    sum = 0.0;
    for (int i = 0 ; i < NVARS; i++){
        sum = sum + x[i]*MOMENTS[i];
    }
    return sum;
}
 
void new_basis_Moments(const int &n, long double *MOMENTS_NEW_BASIS, const long double *MOMENTS, const long double *Sk_cur) {
    for (int i = 0; i < n; i++) {
        MOMENTS_NEW_BASIS[i] = 0.0;
        for (int j = 0; j < n; j++) {
            MOMENTS_NEW_BASIS[i] += Sk_cur[i*n+j]*MOMENTS[j];
        }
    }
} // added by jojo

void Set_Initial_Guess_And_Tol(M2_State_Param &M2_State) {
    long double tolerance_x = 1.0e-32;
    long double coeff_int_Solid_Angle;
    
    switch (M2_State.Domain_Type){
        case BOUNDARY_GAM1:
        case BOUNDARY_GAM2:
        case BOUNDARY_GAM3:
            coeff_int_Solid_Angle = 2.0*PI;
            break;
        case GAM1_GAM2_GAM3:
            coeff_int_Solid_Angle = 4.0*PI;
            break;
        default:
            cout << "Domain type not speficied" << endl;
            exit(0);
            break;
    }
    
    switch (M2_State.Regime){
        case GRAY:
            M2_State.x[0] = -pow(coeff_int_Solid_Angle/M2_State.MOM_VEC[0], 1.0/4.0);
            // M2_State.x[0] = -sqrt(sqrt(coeff_int_Solid_Angle/M2_State.MOM_VEC[0]));
            break;
        case HYPERBOLIC_LIMIT:
            M2_State.x[0] = -log(coeff_int_Solid_Angle/M2_State.MOM_VEC[0]);
            break;
        case BOSE_EINSTEIN:
            M2_State.x[0] = -log(1.0 + coeff_int_Solid_Angle/M2_State.MOM_VEC[0]);
            break;
        case LOGARITHMIC_LIMIT:
            M2_State.x[0] = -coeff_int_Solid_Angle/M2_State.MOM_VEC[0];
            break;
        default:
            cout << "Regime type not specified" << endl;
            exit(0);
            break;
    };
    
    //cout << "x0 = " << M2_State.x[0] << "  " << "roundn(x0) = " << roundn(M2_State.x[0], 6) << endl;
    
    // why does the initial condition affect convergence so much
//     M2_State.x[0] = roundn(M2_State.x[0], 6);
    
    // cout << "M2_State.x[0] = " << M2_State.x[0] << endl;
    
//     if (M2_State.x[0] == 0.0) {
        M2_State.x[0] = -1.0;
//     }
    
    for (int i = 1; i < M2_State.NVARS; i++) {
        M2_State.x[i] = 0.0;  /* some initial guess */
    }
    
    for (int i = 0; i < M2_State.NVARS; i++) {
        M2_State.tol_x[i] = tolerance_x;
    }
}

void Set_Regime(M2_State_Param *M2_State) {
    if (M2_State->Problem_Type == NON_GRAY) {
        if (M2_State->ratio_I0 < -1.0 + 1.0e-6) {
            M2_State->Regime = HYPERBOLIC_LIMIT;
            if (M2_State->display) {
                if (M2_State->id_proc == M2_State->proc_display) {
                    cout << "Hyperbolic Limit" << endl;
                }
            }
            M2_State->MOM_VEC[0] = 1.0e0;
            M2_State->ratio_I0 = -1.0;
        } else if (M2_State->ratio_I0 < 1.0 - 1.0e-6) {
            M2_State->Regime = BOSE_EINSTEIN;
            if (M2_State->display) {
                if (M2_State->id_proc == M2_State->proc_display) {
                    cout << "Bose Einstein statistics" << endl;
                }
            }
        } else {
            M2_State->Regime = LOGARITHMIC_LIMIT;
            if (M2_State->display) { 
                if (M2_State->id_proc == M2_State->proc_display) {
                    cout << "Logarithmic Limit" << endl;
                }
            }
            M2_State->ratio_I0 = 1.0;
            M2_State->MOM_VEC[0] = 1.0e0;
        }
    } else if (M2_State->Problem_Type == GRAY) {
        M2_State->Regime = GRAY;
        M2_State->MOM_VEC[0] = 1.0;
    } else {
        if (M2_State->display) {
            if (M2_State->id_proc == M2_State->proc_display) {
                cout << "Problem Type not specified" << endl;
            }
        }
    }
}

 long double myfunc(unsigned n, const long double *x, const long double *Sk, long double *grad, void *my_func_data) {
     M2_State_Param *M2_State = (M2_State_Param *) my_func_data;
    long double *obj_func_val, *grad_vals, *MOM_VEC_new_basis;
    long double objective_function;
    int NF_obj;
    obj_func_val = new long double[1];
    grad_vals = new long double[M2_State->NVARS];
    MOM_VEC_new_basis = new long double[M2_State->NVARS];
    
    NF_obj = 1;
    
    new_basis_Moments(M2_State->NVARS, MOM_VEC_new_basis, M2_State->MOM_VEC, Sk);
    
    M2_State->set_x(x);
    M2_State->set_Sk(Sk);
    
    if (grad) {
        Lebedev_Quadrature_Func(grad_vals, M2_State->NVARS, gradient, M2_State);
        for (int i = 0; i < M2_State->NVARS; i++) {
            grad[i] = grad_vals[i] - MOM_VEC_new_basis[i];
            
            // cout << "i = " << i << "    "  << "grad_vals[i] = " << grad_vals[i] << "    " << "MOM_VEC_new_basis[i] = " << MOM_VEC_new_basis[i] << "    " << "grad[i] = " << grad[i] << endl;
        }
    }
    
    Lebedev_Quadrature_Func(obj_func_val, NF_obj, obj_function, M2_State);
    
    objective_function = obj_func_val[0] - sum_Lag_Mom(x,M2_State->NVARS,MOM_VEC_new_basis);
//         cout << "objective_function = " << objective_function << endl;
    
    delete[] obj_func_val;
    delete[] grad_vals;
    delete[] MOM_VEC_new_basis;
    
    return objective_function;   
}


// ******************************************************************************
// This routine computes the objective function of the dual maximum entropy 
// problem as well as the gradient of the latter for use in the NLOPT 
// optimization algorithm
// ******************************************************************************
int obj_function(long double *F_obj, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     
     coeff_1 = 1.0;
     
     long double *poly, *poly_Sk;
     poly = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     poly_Sk = new long double[M2_State->NVARS];
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
     }
     
     switch (M2_State->Regime) {
         case GRAY:
             F_obj[0] = (-1.0/3.0)/pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS),3);
             break;
         case BOSE_EINSTEIN:
             F_obj[0] = -coeff_1*log(1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)));
             break;
         case HYPERBOLIC_LIMIT:
             F_obj[0] = coeff_1*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             F_obj[0] = -coeff_1*log(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
// ******************************************************************************
// This routine computes the gradient of the dual maximum entropy 
// problem for use in the NLOPT optimization algorithm
// ******************************************************************************
int gradient(long double *grad, const int &NFUN, const int &Id_angle, void *fdata) {
    int coeff_1;
    M2_State_Param *M2_State = (M2_State_Param *) fdata;
    
    coeff_1 = 1.0;
    
    long double *poly, *poly_Sk;
    poly = new long double[M2_State->NVARS];
     
    generate_polynomials(poly, Id_angle, *M2_State);
     
    poly_Sk = new long double[M2_State->NVARS];
    for (int i = 0; i < M2_State->NVARS; i++) {
        poly_Sk[i] = 0.0;
        for (int j = 0; j < M2_State->NVARS; j++) {
            poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
    }
    
    switch (M2_State->Regime) {
        case GRAY:
            for (int i = 0; i < NFUN; i++) {
                grad[i] = 1.0*poly_Sk[i]/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
            }
            break;
         case BOSE_EINSTEIN:
             for (int i = 0; i < NFUN; i++) {
                 grad[i] = coeff_1*poly_Sk[i]*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))/(1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)));
             }
             break;
         case HYPERBOLIC_LIMIT:
             for (int i = 0; i < NFUN; i++) {
                 grad[i] = coeff_1*poly_Sk[i]*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             }
             break;
          case LOGARITHMIC_LIMIT:
              for (int i = 0; i < NFUN; i++) {
                  grad[i] = coeff_1*poly_Sk[i]/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
              }
              break;
    }
    
    delete[] poly;
    delete[] poly_Sk;
    return 0;
}

// ******************************************************************************
// This routine computes entries of Hessian matrix of second-derivatives of the 
// objective function associated with the dual maximum entropy problem for use 
// in the adaptive change of basis procedure employed in the NLOPT optimization 
// algorithm
// ******************************************************************************
 void Q_func_Gram_Schmidt(unsigned n, const long double *x, const long double *Sk, long double *Q_data, void *my_func_data, void *f_data_MN_State){
     M2_State_Param *M2_State_temp = (M2_State_Param *) f_data_MN_State; 
     NLOPT_State *NLP_State = (NLOPT_State *) my_func_data;  
     M2_State_temp->set_x(x);
     M2_State_temp->set_Sk(Sk);
     M2_State_temp->set_Sk_cur(NLP_State->Sk_cur);
     M2_State_temp->set_ak(NLP_State->ak);
     M2_State_temp->index_p = NLP_State->index_p;
     M2_State_temp->index_a = NLP_State->index_a;
     
     Lebedev_Quadrature_Func(Q_data, 1, H_ONE_Matrix_Entries, M2_State_temp);
 }

int Partial_Moments_DIRECTIONS(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata) {
    long double coeff_1;
    long double poly_N1_1, poly_N1_2, poly_N2_11, poly_N2_12;
    long double Omega1, Omega2, Omega3;
     
    M2_State_Param *M2_State = (M2_State_Param *) fdata;
     
    Omega1 = M2_State->Omega1[Id_angle];
    Omega2 = M2_State->Omega2[Id_angle];
    Omega3 = M2_State->Omega3[Id_angle];
     
     coeff_1 = 1.0;
     
     long double *poly, *poly_Sk;
     poly = new long double[M2_State->NVARS];
     
     generate_polynomials(poly, Id_angle, *M2_State);
     
     poly_Sk = new long double[M2_State->NVARS];
     for (int i = 0; i < M2_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M2_State->NVARS; j++) {
             poly_Sk[i] += M2_State->Sk[i*M2_State->NVARS+j]*poly[j];
        }
	 //cout << "x = " << M2_State->x[i] ;
     }
     poly_N1_1 = Omega1;
     poly_N1_2 = Omega2;
//      if (poly_N1_1 < 0.0)
//          cout << "poly_N1_1 = " << poly_N1_1 << endl;
     
     poly_N2_11 = pow(Omega1,2);
     poly_N2_12 = Omega1*Omega2;
     switch (M2_State->Regime) {
         case GRAY:
             MOMENTS[0] = coeff_1/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[1] = coeff_1*poly_N1_1/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[2] = coeff_1*poly_N1_2/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[3] = coeff_1*poly_N2_11/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             MOMENTS[4] = coeff_1*poly_N2_12/(pow(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS), 4));
             break;
         case BOSE_EINSTEIN:
             MOMENTS[0] = coeff_1*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))/(1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)));
             MOMENTS[1] = coeff_1*poly_N1_1*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))/(1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)));
             MOMENTS[2] = coeff_1*poly_N1_2*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))/(1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)));
             MOMENTS[3] = coeff_1*poly_N2_11*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))/(1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)));
             MOMENTS[4] = coeff_1*poly_N2_12*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS))/(1.0 - exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS)));
             break;
         case LOGARITHMIC_LIMIT:
             MOMENTS[0] = coeff_1/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[1] = coeff_1*poly_N1_1/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[2] = coeff_1*poly_N1_2/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[3] = coeff_1*poly_N2_11/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[4] = coeff_1*poly_N2_12/(-sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
         case HYPERBOLIC_LIMIT:
             MOMENTS[0] = coeff_1*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[1] = coeff_1*poly_N1_1*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[2] = coeff_1*poly_N1_2*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[3] = coeff_1*poly_N2_11*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             MOMENTS[4] = coeff_1*poly_N2_12*exp(sum_Lagrange(M2_State->x,poly_Sk,M2_State->NVARS));
             break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
 void Partial_Moments_M2_GRAY_3D_RT(long double *PARTIAL_MOMS, const long double *x, const long double *Sk, M2_State_Param &M2_State) {
    int NF_PM;
    NF_PM = 5;
    
    M2_State.set_x(x);
    M2_State.set_Sk(Sk);
     
    Partial_Lebedev_Quadrature_Func(PARTIAL_MOMS, NF_PM, Partial_Moments_DIRECTIONS, &M2_State);
    
    
//     if (id == 0) {
// 	    cout << "E_plus = " << PARTIAL_MOMS[0] << "    " << "fx_plus = " << PARTIAL_MOMS[1] << "    " << "pxx_plus = " << PARTIAL_MOMS[2] << "    " << "pxy_plus = " << PARTIAL_MOMS[3] << endl;
//     }
 }
 
 void Partial_Moments_M2_GRAY_3D(long double *PARTIAL_MOMS, const int &NVARS, const int &NFUN, const long double *x, const long double *Sk, const long double &x_Lebed, const long double &y_Lebed, const long double &z_Lebed) {
     M2_State_Param M2_State(NVARS);
     int available_Lebed_Table;
     long double *x_temp;
     int NF_PM;
     NF_PM = 5;
     x_temp = new long double[NVARS];
     M2_State.Lebedev_rule = 65;
     available_Lebed_Table = available_table ( M2_State.Lebedev_rule );
     
     if (available_Lebed_Table) {
         M2_State.Order_Quad_mu = order_table ( M2_State.Lebedev_rule );
         M2_State.Allocate_Quad();
    } else {
        cout << "No Table available for specified Lebedev Rule" << endl;
        exit(0);
    }
    
//     Rotate_Frame(x_temp, x, NVARS, x_Lebed, y_Lebed, z_Lebed);
    
    M2_State.set_x(x_temp);
    M2_State.set_Sk(Sk);
     
    Partial_Lebedev_Quadrature_Func(PARTIAL_MOMS, NF_PM, Partial_Moments_DIRECTIONS, &M2_State);
    
    M2_State.Deallocate_Quad();
    delete[] x_temp;
 }
 
long double Lebedev_Quadrature_Matrix(const long double *Matrix, const int &Quad_Rule) {
    int order;
    long double *x, *y, *z, *w;
    long double integral_approx;
    
    order = order_table ( Quad_Rule );
    
    x = new long double[order];
    y = new long double[order];
    z = new long double[order];
    w = new long double[order];
    
    cout << "Fix this first !!!!!!!!!!!!!!" << endl;
    
    // ld_by_order ( order, x, y, z, w );
    
    integral_approx = 0.0;
    
    for ( int m = 0; m < order; m++ ) {
        integral_approx = integral_approx + w[m] * Matrix[m];
    }
    
    integral_approx = 4.0*PI*integral_approx;
    
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] w;
    
    return integral_approx;
}

void Partial_Lebedev_Quadrature_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata) {
//     int order_lebed;
//     long double Temp_vals[NFUN];
//     M2_State_Param *M2_State = (M2_State_Param *) fdata;
//     order_lebed = M2_State->Order_Quad_mu;
//     for (int i = 0; i < NFUN; i++) {
//         Vals[i] = 0.0;
//     }
//     
//     for ( int m = 0; m < order_lebed; m++ ) {
//         if (M2_State->x_Lebed[m] >= 0.0) {
//             func(Temp_vals, NFUN, m, fdata);
//             for (int i = 0; i < NFUN; i++) {
//                 Vals[i] = Vals[i] + M2_State->w_Lebed[m] * Temp_vals[i];
//             }
//         }
//     }
//     
//     for (int i = 0; i < NFUN; i++) {
//         Vals[i] = 4.0*PI*Vals[i];
//     }
//     
}

// void Lebedev_Quadrature_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata) {
//     int order_lebed;
//     int Index_Angle;
//     long double Temp_vals[NFUN];
//     M2_State_Param *M2_State = (M2_State_Param *) fdata;
//     order_lebed = M2_State->Order_Quad_mu;
//     for (int i = 0; i < NFUN; i++) {
//         Vals[i] = 0.0;
//     }
//     
//     switch (M2_State->Domain_Type) {
//         case GAM1_GAM2_GAM3:
//             for ( int Id_Phi = 0; Id_Phi < 2*order_lebed; Id_Phi++ ) {
//                 for ( int Id_mu = 0; Id_mu < order_lebed; Id_mu++ ) {
//                     Index_Angle = M2_State->Id_Angle(Id_mu, Id_Phi);
//                     func(Temp_vals, NFUN, Index_Angle, fdata);
//                     for (int i = 0; i < NFUN; i++) {
//                         Vals[i] = Vals[i] + M2_State->weight_total[Index_Angle] * Temp_vals[i];
//                     }
//                 }
//             }
//             break;
//         case BOUNDARY_GAM1:
//         case BOUNDARY_GAM2:
//         case BOUNDARY_GAM3:
//             for ( int Id_Phi = 0; Id_Phi < 2*order_lebed; Id_Phi++ ) {
//                 func(Temp_vals, NFUN, Id_Phi, fdata);
//                 for (int i = 0; i < NFUN; i++) {
//                     Vals[i] = Vals[i] + M2_State->w_phi_quad[Id_Phi] * Temp_vals[i];
//                 }
//             }
//             break;
//         default:
//             break;
//     }
// }

void Lebedev_Quadrature_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata) {
    int order_lebed;
    long double Temp_vals[NFUN];
    M2_State_Param *M2_State = (M2_State_Param *) fdata;
    order_lebed = M2_State->Order_Quad_mu;
    for (int i = 0; i < NFUN; i++) {
        Vals[i] = 0.0;
    }
    
    long double mu_start, mu_end;
    long double phi_start, phi_end;
    int Index_Angle;
    
    switch (M2_State->Domain_Type) {
        case GAM1_GAM2_GAM3:
            for (int i_phi_quad = 0; i_phi_quad < M2_State->N_intervals_phi_quad; i_phi_quad++) {
                phi_start = M2_State->bounds_quad_phi_refin[i_phi_quad];
                phi_end = M2_State->bounds_quad_phi_refin[i_phi_quad+1];
                
                // if (M2_State->refinement_level > 1)
                //    cout << "i_phi_quad = " << i_phi_quad << "  " << "N_intervals_phi_quad = " << M2_State->N_intervals_phi_quad << "  " << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << endl;
                
                for (int i_mu_quad = 0; i_mu_quad < M2_State->N_intervals_mu_quad; i_mu_quad++) {
                    mu_start = M2_State->bounds_quad_mu_refin[i_mu_quad];
                    mu_end = M2_State->bounds_quad_mu_refin[i_mu_quad+1];
                    
                    M2_State->compute_Omegas(mu_start, mu_end, phi_start, phi_end);
                    
//                     if (M2_State->refinement_level > 2) {
//                        cout << "i_mu_quad = " << i_mu_quad << "  " << "N_intervals_mu_quad = " << M2_State->N_intervals_mu_quad << "  " << "mu_start = " << mu_start << "  " << "mu_end = " << mu_end << "  " << endl;
//                     }
                    
                    for ( int Id_Phi = 0; Id_Phi < 2*order_lebed; Id_Phi++ ) {
                        for ( int Id_mu = 0; Id_mu < order_lebed; Id_mu++ ) {
                            Index_Angle = M2_State->Id_Angle(Id_mu, Id_Phi);
                            
                            func(Temp_vals, NFUN, Index_Angle, fdata);
                            for (int i = 0; i < NFUN; i++) {
                                Vals[i] = Vals[i] + M2_State->weight_total[Index_Angle] * Temp_vals[i];
                            }
                        }
                    }
                }
            }
            break;
        case BOUNDARY_GAM1:
        case BOUNDARY_GAM2:
        case BOUNDARY_GAM3:
            for (int i_phi_quad = 0; i_phi_quad < M2_State->N_intervals_phi_quad; i_phi_quad++) {
                phi_start = M2_State->bounds_quad_phi_refin[i_phi_quad];
                phi_end = M2_State->bounds_quad_phi_refin[i_phi_quad+1];
                
                M2_State->compute_Omegas(mu_start, mu_end, phi_start, phi_end);
                
                // if (M2_State->refinement_level > 0)
                //    cout << "i_phi_quad = " << i_phi_quad << "  " << "N_intervals_phi_quad = " << M2_State->N_intervals_phi_quad << "  " << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << endl;
                    
                for ( int Id_Phi = 0; Id_Phi < 2*order_lebed; Id_Phi++ ) {
                    func(Temp_vals, NFUN, Id_Phi, fdata);
                    for (int i = 0; i < NFUN; i++) {
                        Vals[i] = Vals[i] + M2_State->weight_total[Id_Phi] * Temp_vals[i];
                    }
                }
            }
            break;
        default:
            break;
    }
}
 
long double Inverse_Mobius_Transformation(const long double &ratio, const long double &Length_Scale) {
     long double Val;
     Val = -Length_Scale*log((1.0 - ratio)/2.0);
//      Val = Length_Scale*(ratio + 1.0)/(1.0 - ratio);
     return Val;
 }
 
 long double Mobius_Transformation(long double &Moment, const long double &Length_Scale) {
     long double ratio;
     ratio = 1.0 - 2.0*exp(-Moment/Length_Scale);
//      ratio = (Moment - Length_Scale)/(Moment + Length_Scale);
     return ratio;
 }

//  long double myconstraint(unsigned n, const long double *x, const long double *Sk, long double *grad, void *data_constraints) {
//     long double temp_val = 0.0;
//     long double constr_val = 0.0, constr_val_max = -1.0e12;
//     my_constraint_data data_max, data;
//     M2_State_Param *M2_State_temp = (M2_State_Param *) data_constraints;
//     
//     int Id_Angle;
//     long double *poly;
//     poly = new long double[M2_State_temp->NVARS];
//     
//     if (grad) {
//         for (int i = 0; i < n; i++) {
//             grad[i] = 0.0;
//         }
//     }
//     
//     long double phi_start, phi_end;
//     long double mu_start, mu_end;
//     switch (M2_State_Param::Regime) {
//         case GRAY:
//         case BOSE_EINSTEIN:
//         case LOGARITHMIC_LIMIT:
//             for (int Id_Angle = 0; Id_Angle < M2_State_temp->N_dirs(); Id_Angle++) {
//                         generate_polynomials(poly, Id_Angle, *M2_State_temp);
//                         if (M2_State_temp->NVARS == 5) {
//                             data.a = poly[0];
//                             data.b = poly[1];
//                             data.c = poly[2];
//                             data.f = poly[3];
//                             data.g = poly[4];
//                             data.h = 0.0;
//                             data.i = 0.0;
//                             data.j = 0.0;
//                             data.k = 0.0;  
//                         } else if (M2_State_temp->NVARS == 9) {
//                             data.a = poly[0];
//                             data.b = poly[1];
//                             data.c = poly[2];
//                             data.f = poly[3];
//                             data.g = poly[4];
//                             data.h = poly[5];
//                             data.i = poly[6];
//                             data.j = poly[7];
//                             data.k = poly[8];
//                         }   
//                         
//                         constr_val = 0.0;
//                         for (int i = 0; i < n; i++) {
//                             for (int j = 0; j < n; j++) {
//                                 constr_val += Sk[i*n+j]*data[j]*x[i];
//                             }
//                         }
//                         
//                         constr_val_max = max(constr_val_max, constr_val);
//                         if (constr_val_max == constr_val) {
//                             data_max = data;
//                         }
//             }
//             
//             if (grad) {
//                 for (int i = 0; i < n; i++) {
//                     for (int j = 0; j < n; j++) {
//                         grad[i] += Sk[i*n+j]*data_max[j];
//                     }
//                     grad[i] *= 1.0e6;
//                 }
//             }
//             
//             temp_val = 1.0e6*constr_val_max;
//             break;
//         case HYPERBOLIC_LIMIT:
//             // Positivity of the distribution does not impose any constraint
//             // on the Lagrange multipliers for the Boltzmann distribution
//             if (grad) {
//                 for (int i = 0; i < n; i++) {
//                     grad[i] = 0.0;
//                 }
//             }
//             temp_val = 0.0;
//             break;
//     }
//     
//     delete[] poly;
//             
//     return temp_val;
//  }
 
long double myconstraint(unsigned n, const long double *x, const long double *Sk, long double *grad, void *data_constraints) {
    long double temp_val = 0.0;
    long double constr_val = 0.0, constr_val_max = -1.0e12;
    my_constraint_data data_max, data;
    M2_State_Param *M2_State_temp = (M2_State_Param *) data_constraints;
    
    int Id_Angle;
    long double *poly;
    poly = new long double[M2_State_temp->NVARS];
    
    long double phi_start, phi_end;
    long double mu_start, mu_end;
    switch (M2_State_Param::Regime) {
        case GRAY:
        case BOSE_EINSTEIN:
        case LOGARITHMIC_LIMIT:
            for (int i_phi_quad = 0; i_phi_quad < M2_State_temp->N_intervals_phi_quad; i_phi_quad++) {
                phi_start = M2_State_temp->bounds_quad_phi_refin[i_phi_quad];
                phi_end = M2_State_temp->bounds_quad_phi_refin[i_phi_quad+1];
                
                // cout << "i_phi_quad = " << i_phi_quad << "  " << "N_intervals_phi_quad = " << M2_State_temp->N_intervals_phi_quad << "  " << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << endl;
                
                for (int i_mu_quad = 0; i_mu_quad < M2_State_temp->N_intervals_mu_quad; i_mu_quad++) {
                    mu_start = M2_State_temp->bounds_quad_mu_refin[i_mu_quad];
                    mu_end = M2_State_temp->bounds_quad_mu_refin[i_mu_quad+1];
                    
                    // cout << "i_mu_quad = " << i_mu_quad << "  " << "N_intervals_mu_quad = " << M2_State_temp->N_intervals_mu_quad << "  " << "mu_start = " << mu_start << "  " << "mu_end = " << mu_end << "  " << endl;
                    
                    
                    M2_State_temp->compute_Omegas(mu_start, mu_end, phi_start, phi_end);
                    
                    for (int Id_Angle = 0; Id_Angle < M2_State_temp->N_dirs(); Id_Angle++) {
                        generate_polynomials(poly, Id_Angle, *M2_State_temp);
                        if (M2_State_temp->NVARS == 5) {
                            data.a = poly[0];
                            data.b = poly[1];
                            data.c = poly[2];
                            data.f = poly[3];
                            data.g = poly[4];
                            data.h = 0.0;
                            data.i = 0.0;
                            data.j = 0.0;
                            data.k = 0.0;  
                        } else if (M2_State_temp->NVARS == 9) {
                            data.a = poly[0];
                            data.b = poly[1];
                            data.c = poly[2];
                            data.f = poly[3];
                            data.g = poly[4];
                            data.h = poly[5];
                            data.i = poly[6];
                            data.j = poly[7];
                            data.k = poly[8];
                        }   
                        
                        constr_val = 0.0;
                        for (int i = 0; i < n; i++) {
                            for (int j = 0; j < n; j++) {
                                constr_val += Sk[i*n+j]*data[j]*x[i];
                            }
                        }
                        
                        constr_val_max = max(constr_val_max, constr_val);
                        if (constr_val_max == constr_val) {
                            data_max = data;
                            
                        }
                    }
                }
            }
            
            if (grad) {
                for (int i = 0; i < n; i++) {
                    grad[i] = 0.0;
                    // Since the location of constr_val_max changes, the gradient should be set
                    // to zero to avoid convergence issues??
//                     for (int j = 0; j < n; j++) {
//                         grad[i] += Sk[i*n+j]*data_max[j];
//                     }
//                     // Check what is wrong with SLSQP sign convention for gradients of inequality constraints ??
//                     grad[i] = -grad[i];
//                     // grad[i] *= 1.0e6;
                }
            }
            
            temp_val = constr_val_max + 1.0e-22;
            // temp_val *= 1.0e6;
            break;
        case HYPERBOLIC_LIMIT:
            // Positivity of the distribution does not impose any constraint
            // on the Lagrange multipliers for the Boltzmann distribution
            if (grad) {
                for (int i = 0; i < n; i++) {
                    grad[i] = 0.0;
                }
            }
            temp_val = 0.0;
            break;
        default:
            cout << "Invalid value specified for M2_State_Param::Regime !!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    delete[] poly;
    
    // cout << "temp_val = " << temp_val << endl;
            
    return temp_val;
 }
 
 long double Rotation_Matrix_X_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle) {
     long double mat = 0.0;
     switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = 1.0;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 1:
                    mat = cos_angle;
                    break;
                case 2:
                    mat = -sin_angle;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 1:
                    mat = sin_angle;
                    break;
                case 2:
                    mat = cos_angle;
                    break;
            };
            break;
    };
    
    return mat;
}

long double Rotation_Matrix_Y_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle) {
    long double mat = 0.0;
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = cos_angle;
                    break;
                case 2:
                    mat = -sin_angle;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 1:
                    mat = 1.0;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 0:
                    mat = sin_angle;
                    break;
                case 2:
                    mat = cos_angle;
                    break;
            };
            break;
    };
    
    return mat;
} 

long double Rotation_Matrix_Z_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle) {
    long double mat = 0.0;
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = cos_angle;
                    break;
                case 1:
                    mat = sin_angle;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 0:
                    mat = -sin_angle;
                    break;
                case 1:
                    mat = cos_angle;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 2:
                    mat = 1.0;
                    break;
            };
            break;
    };
    
    return mat;
}

void Rotate_Frame_N3_Basis(long double &Omega1, long double &Omega2, long double &Omega3, const long double &theta, const long double &phi) {
    long double *x_rot_interm, *x_rot_temp;;
    long double cos_angle, sin_angle;
    
    x_rot_interm = new long double[3];
    x_rot_temp = new long double[3];
    
    x_rot_interm[0] = Omega1;
    x_rot_interm[1] = Omega2;
    x_rot_interm[2] = Omega3;
    
    cos_angle = cos(-phi);
    sin_angle = sin(-phi);
    
    // Rotating Vector around X axis
    for (int i = 0; i < 3; i++) {
        x_rot_temp[i] = 0.0;
        for (int j = 0; j < 3; j++) {
//             cout << "i = " << i << "   " << "j = " << j << "   " << "Rotation_Matrix_X_axis = " << Rotation_Matrix_X_axis(i, j, cos_angle, sin_angle) << endl;
            x_rot_temp[i] += x_rot_interm[j]*Rotation_Matrix_X_axis(i, j, cos_angle, sin_angle);
        }
//         cout << "i = " << i << "   " << "x_rot_interm[1+i] = " << x_rot_interm[1+i] << "   " << "M2_State.MOM_VEC[1+i] = " << M2_State.MOM_VEC[1+i] << endl;
    }
    
    
    cos_angle = cos(-theta);
    sin_angle = sin(-theta);
    
    // Rotating Vector around Z axis
    for (int i = 0; i < 3; i++) {
        x_rot_interm[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            x_rot_interm[i] += x_rot_temp[j]*Rotation_Matrix_Z_axis(i, j, cos_angle, sin_angle);
        }
    }
    
    Omega1 = x_rot_interm[0];
    Omega2 = x_rot_interm[1];
    Omega3 = x_rot_interm[2];
    
    delete[] x_rot_interm;
    delete[] x_rot_temp;
}

void Rotate_Frame(M2_State_Param &M2_State, const long double &theta, const long double &phi) {
    long double *x_rot_interm;
    long double cos_angle, sin_angle;
    int index_mat_to_vec, index_temp;
    long double trace;
    
    cout << "theta = " << theta << "   " << "phi = " << phi << endl;
    
    x_rot_interm = new long double[M2_State.NVARS];
    
    cos_angle = cos(-phi);
    sin_angle = sin(-phi);
    
    // Rotating Vector around X axis
    for (int i = 0; i < 3; i++) {
        x_rot_interm[1+i] = 0.0;
        for (int j = 0; j < 3; j++) {
//             cout << "i = " << i << "   " << "j = " << j << "   " << "Rotation_Matrix_X_axis = " << Rotation_Matrix_X_axis(i, j, cos_angle, sin_angle) << endl;
            x_rot_interm[1+i] += M2_State.MOM_VEC[1+j]*Rotation_Matrix_X_axis(i, j, cos_angle, sin_angle);
        }
//         cout << "i = " << i << "   " << "x_rot_interm[1+i] = " << x_rot_interm[1+i] << "   " << "M2_State.MOM_VEC[1+i] = " << M2_State.MOM_VEC[1+i] << endl;
    }
    
    // Rotating Second-order tensor around X axis
    trace = 1.0 - M2_State.MOM_VEC[4] - M2_State.MOM_VEC[7];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i < 2 || j < 2) {
                index_mat_to_vec = Matrix_to_Vector_Indexing(i, j, 3);
                
//                 cout << "index_mat_to_vec = " << index_mat_to_vec << endl;
                
                x_rot_interm[4+index_mat_to_vec] = 0.0;
                
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        index_temp = Matrix_to_Vector_Indexing(k, l, 3);
                        if (k == l && k == 2) {
                            x_rot_interm[4+index_mat_to_vec] += trace*Rotation_Matrix_X_axis(i, k, cos_angle, sin_angle)*Rotation_Matrix_X_axis(j, l, cos_angle, sin_angle);
                        } else {
                            x_rot_interm[4+index_mat_to_vec] += M2_State.MOM_VEC[4+index_temp]*Rotation_Matrix_X_axis(i, k, cos_angle, sin_angle)*Rotation_Matrix_X_axis(j, l, cos_angle, sin_angle);
                        }
                    }
                }
            }
        }
    }
    
    cos_angle = cos(-theta);
    sin_angle = sin(-theta);
    
    // Rotating Vector around Z axis
    for (int i = 0; i < 3; i++) {
        M2_State.MOM_VEC[1+i] = 0.0;
        for (int j = 0; j < 3; j++) {
            M2_State.MOM_VEC[1+i] += x_rot_interm[1+j]*Rotation_Matrix_Z_axis(i, j, cos_angle, sin_angle);
        }
    }
    
    // Rotating Second-order tensor around Z axis
    trace = 1.0 - x_rot_interm[4] - x_rot_interm[7];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i < 2 || j < 2) {
                index_mat_to_vec = Matrix_to_Vector_Indexing(i, j, 3);
                M2_State.MOM_VEC[4+index_mat_to_vec] = 0.0;
                
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        index_temp = Matrix_to_Vector_Indexing(k, l, 3);
                        if (k == l && k == 2) {
                            M2_State.MOM_VEC[4+index_mat_to_vec] += trace*Rotation_Matrix_Z_axis(i, k, cos_angle, sin_angle)*Rotation_Matrix_Z_axis(j, l, cos_angle, sin_angle);
                        } else {
                            M2_State.MOM_VEC[4+index_mat_to_vec] += x_rot_interm[4+index_temp]*Rotation_Matrix_Z_axis(i, k, cos_angle, sin_angle)*Rotation_Matrix_Z_axis(j, l, cos_angle, sin_angle);
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < M2_State.NVARS; i++) {
        cout << "i = " << i << "   " << "x_rot_interm[i] = " << x_rot_interm[i] << "   " << "M2_State.MOM_VEC[i] = " << M2_State.MOM_VEC[i] << endl;
    }
    
    delete[] x_rot_interm;
}

int Matrix_to_Vector_Indexing(const int &i, const int &j, const int &NVARS) {
    int index_vec;
    
    if (i > j) {
        index_vec = NVARS - (NVARS - j)*(NVARS - 1 - j)/2 + i;
    } else { // In this case i switches to j
        index_vec = NVARS - (NVARS - i)*(NVARS - 1 - i)/2 + j;
    }
    
    return index_vec;
}

void Modified_Gram_Schmidt_Factorization(const int &n, const long double *A, long double *Q, long double *R) {
    // Modified Gram Schmidt algorithm wih reorthogonalization
    long double Q_val_ak_pm, Q_val_pm2, *Sk_prev;
    long double norm_a;
    long double proj_Q_A, proj_Q_Q;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Q[i*n+j] = 0.0;
        }
    }

    for (int k = 0; k < n; k++){ // k represents the number of polynomials in the non-orthogonal basis
        for (int i = 0; i < n; i++) {
            Q[i*n + k] = A[i*n + k];
        }  // end for loop i
        for (int l = 0; l < 2; l++) { // loop for reorthogonalization
            for(int m = 0; m < k; m++) {
                // evaluate norm`
                proj_Q_A = 0.0;
                proj_Q_Q = 0.0;
                for (int i = 0; i < n; i++) {
                    proj_Q_A += Q[i*n + k]*Q[i*n + m];
                    proj_Q_Q += Q[i*n + m]*Q[i*n + m];
                }  // end for loop i
                for (int i = 0; i < n; i++) {
                    // orthogonalize k^{th} column of A with respect to columns up to (k-1) of Q
                    Q[i*n + k] = Q[i*n + k] - (proj_Q_A/proj_Q_Q)*Q[i*n + m];
                }  // end for loop i
            } // end for loop m
        } // end for loop l
        
        proj_Q_Q = 0.0;
        for (int i = 0; i < n; i++) {
            proj_Q_Q += Q[i*n + k]*Q[i*n + k];
        }  // end for loop i
        for (int i = 0; i < n; i++) {
            Q[i*n+k] = Q[i*n+k]/sqrt(proj_Q_Q); // Orthonormalize
//             cout << "i = " << i << "   " << "k = " << k << "   " << "Q = " << Q[i*n+k] << "   " << "proj_Q_Q = " << proj_Q_Q << endl;
        }
    } // end for loop k
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            R[i*n+j] = 0.0;
            for (int k = 0; k < n; k++) {
                R[i*n+j] += Q[k*n+i]*A[k*n+j];
            }
//             cout << "i = " << i << "   " << "j = " << j << "   " << "R = " << R[i*n+j] << endl;
        }
    }
}

void setup_spherical_harmonics_data (const int &N_Points_Phi, const int &N_Points_Theta, long double *x_Lebed, long double *y_Lebed, long double *z_Lebed, long double *w_quad) {
    long double *phi_quad, *w_phi_quad, *mu_quad, *w_Leg;
    int index;
    phi_quad = new long double[N_Points_Phi];
    w_phi_quad = new long double[N_Points_Phi];
    mu_quad = new long double[N_Points_Theta];
    w_Leg = new long double[N_Points_Theta];
        
    circle_rule ( N_Points_Phi, w_phi_quad, phi_quad );
    legendre_set ( N_Points_Theta, mu_quad, w_Leg );
    
    // cout << "Consider Exploiting evenness of closing fluxes with respect to spherical harmonics" << endl;
    
    for (int Id_Phi = 0; Id_Phi < N_Points_Phi; Id_Phi++) {
        for (int Id_mu = 0; Id_mu < N_Points_Theta; Id_mu++) {
            index = Id_Phi*N_Points_Theta + Id_mu;
            z_Lebed[index] = mu_quad[Id_mu];
            x_Lebed[index] = sqrt(1.0 - pow(mu_quad[Id_mu],2))*cos(phi_quad[Id_Phi]);
            y_Lebed[index] = sqrt(1.0 - pow(mu_quad[Id_mu],2))*sin(phi_quad[Id_Phi]);
            if (w_quad != NULL) {
                w_quad[index] = 2.0*PI*w_phi_quad[Id_Phi]*w_Leg[Id_mu];
            }
        }
    }
    
    delete[] phi_quad;
    delete[] mu_quad;
    delete[] w_phi_quad;
    delete[] w_Leg;
}

// void setup_spherical_harmonics_data (const int &N_Points_Theta, long double *x_Lebed, long double *y_Lebed, long double *z_Lebed, long double *w_quad) {
//     long double *phi_quad, *w_phi_quad, *mu_quad, *w_Leg;
//     int index;
//     phi_quad = new long double[2*N_Points_Theta];
//     w_phi_quad = new long double[2*N_Points_Theta];
//     mu_quad = new long double[2*(N_Points_Theta - 1) + 1];
//     w_Leg = new long double[2*(N_Points_Theta - 1) + 1];
//         
//     circle_rule ( 2*N_Points_Theta, w_phi_quad, phi_quad );
//     legendre_set ( 2*(N_Points_Theta - 1) + 1, mu_quad, w_Leg );
//     
//     M2_State_Param M2_State(9);
//     
//     // cout << "Consider Exploiting evenness of closing fluxes with respect to spherical harmonics" << endl;
//     int Id_mu_pos;
//     long double phi_val, dphi_val;
//     for (int Id_Phi = 0; Id_Phi < 2*N_Points_Theta; Id_Phi++) {
//         for (int Id_mu = 0; Id_mu < N_Points_Theta; Id_mu++) {
//             Id_mu_pos = (N_Points_Theta - 1) + Id_mu;
//             index = Id_Phi*N_Points_Theta + Id_mu;
//             
//             phi_val = M2_State.linear_interpolation(0.0, PI, phi_quad[Id_Phi], 0.0, 2.0*PI);
//             dphi_val = M2_State.diff_linear_interpolation(0.0, PI, 0.0, 2.0*PI);
//             
//             // cout << "Id_Phi = " << Id_Phi << "  " << "Id_mu = " << Id_mu << "  " << "mu = " << mu_quad[Id_mu_pos] << "  " << "phi_val = " << phi_val << endl;
//             
//             x_Lebed[index] = sqrt(1.0 - pow(mu_quad[Id_mu_pos],2))*cos(phi_val);
//             y_Lebed[index] = sqrt(1.0 - pow(mu_quad[Id_mu_pos],2))*sin(phi_val);
//             z_Lebed[index] = mu_quad[Id_mu_pos];
//             if (w_quad != NULL) {
//                 w_quad[index] = 2.0*PI*dphi_val*w_phi_quad[Id_Phi]*w_Leg[Id_mu_pos];
//             }
//         }
//     }
//     
//     delete[] phi_quad;
//     delete[] mu_quad;
//     delete[] w_phi_quad;
//     delete[] w_Leg;
// }

long double roundval( const long double &val ) {
    if ( val < 0.0 ) {
        return -floor(fabs(val) + 0.5);
    }
    return floor(val + 0.5);
}

long double roundn( const long double &val, const int &n ) {
    long double value, temp_val;
    long double prec;
    prec = pow(10.0, n);
    temp_val = val * prec;
    
    value = roundval( temp_val ) / prec;
    
    return value;
}

void Rotate_Moments(M2_State_Param &M2_State, const long double &phi) {
//     double fx_val, fy_val, pxx_val, pxy_val, pyy_val;
//     double cos_tmp, sin_tmp;
//     cos_tmp = cos(phi);
//     sin_tmp = sin(phi);
//     
//     fx_val = M2_State.N1_1;
//     fy_val = M2_State.N1_2;
//     pxx_val = M2_State.N2_11;
//     pxy_val = M2_State.N1_1*M2_State.N1_2;
//     pyy_val = M2_State.N2_22;
//     
//     // Radiative energy density remains the same
//     
//     // Rotate flux vector
//     M2_State.N1_1 = fx_val*cos_tmp + fy_val*sin_tmp;
//     M2_State.N1_2 = -fx_val*sin_tmp + fy_val*cos_tmp;
//     
//     // Rotate Pressure tensor
//     M2_State.N2_11 = pxx_val*sqr(cos_tmp) + 2.0*cos_tmp*sin_tmp*pxy_val + pyy_val*sqr(sin_tmp);
//     M2_State.N2_12 = -cos_tmp*sin_tmp*pxx_val + (sqr(cos_tmp)-sqr(sin_tmp))*pxy_val + cos_tmp*sin_tmp*pyy_val;
//     M2_State.N2_22 = pxx_val*sqr(sin_tmp) - 2.0*cos_tmp*sin_tmp*pxy_val + pyy_val*sqr(cos_tmp);
}
