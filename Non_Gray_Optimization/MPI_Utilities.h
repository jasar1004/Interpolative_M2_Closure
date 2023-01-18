#ifndef _MPI_UTILITIES_H_INCLUDED
#define _MPI_UTILITIES_H_INCLUDED

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

#ifndef _M2_STATE_PARAMETERS_H_INCLUDED
#include "./M2_State_Parameters.h"
#endif // _M2_STATE_PARAMETERS_H_INCLUDED

using namespace std;
using namespace nlopt;

struct MPI_proc_params {
    int num_proc_per_var_E = 1, num_proc_per_var_f = 1, num_proc_per_var_phi = 1, num_proc_per_var_theta = 1, num_proc_per_var_Triangle_gam1_gam2 = 1;
    int id_proc_E = 0, id_proc_f = 0, id_proc_phi = 0, id_proc_theta = 0, id_proc_Triangle_gam1_gam2 = 0;
    int size_rec = 1;
    int N_Points_E_per_proc, N_Points_f_per_proc, N_Points_Phi_per_proc, N_Points_Theta_per_proc, N_Points_Triangle_gam1_gam2_Per_Proc;
    
    void Copy(const MPI_proc_params &MPI_proc_parameters) {
        num_proc_per_var_E = MPI_proc_parameters.num_proc_per_var_E;
        num_proc_per_var_f = MPI_proc_parameters.num_proc_per_var_f;
        num_proc_per_var_phi = MPI_proc_parameters.num_proc_per_var_phi;
        num_proc_per_var_theta = MPI_proc_parameters.num_proc_per_var_theta;
        num_proc_per_var_Triangle_gam1_gam2 = MPI_proc_parameters.num_proc_per_var_Triangle_gam1_gam2;
        
        id_proc_E = MPI_proc_parameters.id_proc_E;
        id_proc_f = MPI_proc_parameters.id_proc_f;
        id_proc_phi = MPI_proc_parameters.id_proc_phi;
        id_proc_theta = MPI_proc_parameters.id_proc_theta;
        id_proc_Triangle_gam1_gam2 = MPI_proc_parameters.id_proc_Triangle_gam1_gam2;
        
        size_rec = MPI_proc_parameters.size_rec;
        
        N_Points_E_per_proc = MPI_proc_parameters.N_Points_E_per_proc;
        N_Points_f_per_proc = MPI_proc_parameters.N_Points_f_per_proc;
        N_Points_Phi_per_proc = MPI_proc_parameters.N_Points_Phi_per_proc;
        N_Points_Theta_per_proc = MPI_proc_parameters.N_Points_Theta_per_proc;
        N_Points_Triangle_gam1_gam2_Per_Proc = MPI_proc_parameters.N_Points_Triangle_gam1_gam2_Per_Proc;
    }
};

///////////////////////////////////////////////////////
// Routines for MPI Utilities
///////////////////////////////////////////////////////
int Compute_Npts_Total(const MN_Var_Num_Points *num_points);

void Setup_MPI_Processes(MPI_proc_params &MPI_proc_parameters, M2_State_Param &M2_State, const MN_Var_Num_Points *num_points);
void Compute_Id_gam1_gam2_Triangle(int &i_gam1, int &i_gam2, const int &i_Triangle_gam1_gam2, const MN_Var_Num_Points *num_points);  
void Compute_MPI_Processes_Max_Min_Indexes(const MPI_proc_params &MPI_proc_parameters, const MN_Var_Num_Points *num_points, int &id_E_min, int &id_E_max, int &id_f_min, int &id_f_max, int &id_phi_min, int &id_phi_max, int &id_theta_min, int &id_theta_max, int &id_gam1_min, int &id_gam1_max);

void Create_MPI_Data_Type_rec_N3(MPI_Datatype &rec_N3_type);

#endif // _MPI_UTILITIES_H_INCLUDED
