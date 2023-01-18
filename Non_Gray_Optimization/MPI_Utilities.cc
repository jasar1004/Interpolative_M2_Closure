#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#include "NG_MN_Model_3D_OPTIM.h"
#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED

//*****************************************************************
// Compute total number of points
//*****************************************************************
int Compute_Npts_Total(const MN_Var_Num_Points *num_points) {
    int Npts_total;
    
    Npts_total = num_points->E*num_points->f*num_points->theta*num_points->phi*num_points->gam1*(num_points->gam1+1)/2;
    
    return Npts_total;
}

//******************************************************************************************
// Setup parameters required for the purpose of parallel computing
//******************************************************************************************
void Setup_MPI_Processes(MPI_proc_params &MPI_proc_parameters, M2_State_Param &M2_State, const MN_Var_Num_Points *num_points) {
    int num_proc_per_var_E, num_proc_per_var_f, num_proc_per_var_phi, num_proc_per_var_theta, num_proc_per_var_Triangle_gam1_gam2;
    int id_proc_E, id_proc_f, id_proc_phi, id_proc_theta, id_proc_Triangle_gam1_gam2;
    int size_rec;
    
    if ((M2_State.Problem_Type == NON_GRAY) && 
        (num_points->Maximum_Entropy_Solution_Regime != HYPERBOLIC_LIMIT) && 
        (num_points->Maximum_Entropy_Solution_Regime != LOGARITHMIC_LIMIT)) {
        
        if (M2_State.num_proc > 1) {
            num_proc_per_var_E = M2_State.num_proc_E;
            num_proc_per_var_f = M2_State.num_proc_f;
            num_proc_per_var_phi = M2_State.num_proc_Phi;
            num_proc_per_var_theta = M2_State.num_proc_Theta;
            num_proc_per_var_Triangle_gam1_gam2 = M2_State.num_proc_Triangle_gam1_gam2;
        } else {
            num_proc_per_var_E = 1;
            num_proc_per_var_f = 1;
            num_proc_per_var_phi = 1;
            num_proc_per_var_theta = 1;
            num_proc_per_var_Triangle_gam1_gam2 = 1;
        }
        
        // Determine id for each variable based on the current id_proc
        id_proc_E = floor(M2_State.id_proc/(num_proc_per_var_f*num_proc_per_var_theta*num_proc_per_var_phi*num_proc_per_var_Triangle_gam1_gam2));
        id_proc_f = floor(M2_State.id_proc/(num_proc_per_var_theta*num_proc_per_var_phi*num_proc_per_var_Triangle_gam1_gam2)) - id_proc_E*num_proc_per_var_f;
        id_proc_phi = floor(M2_State.id_proc/(num_proc_per_var_theta*num_proc_per_var_Triangle_gam1_gam2)) - ((id_proc_E*num_proc_per_var_f+id_proc_f)*num_proc_per_var_phi);
        id_proc_theta = floor(M2_State.id_proc/num_proc_per_var_Triangle_gam1_gam2) - ((id_proc_E*num_proc_per_var_f+id_proc_f)*num_proc_per_var_phi+id_proc_phi)*num_proc_per_var_theta;
        id_proc_Triangle_gam1_gam2 = M2_State.id_proc - (((id_proc_E*num_proc_per_var_f + id_proc_f)*num_proc_per_var_phi + id_proc_phi)*num_proc_per_var_theta+id_proc_theta)*num_proc_per_var_Triangle_gam1_gam2;
    } else if ((M2_State.Problem_Type == GRAY) ||
               (num_points->Maximum_Entropy_Solution_Regime != HYPERBOLIC_LIMIT) ||
               (num_points->Maximum_Entropy_Solution_Regime != LOGARITHMIC_LIMIT)) {
        
        if (M2_State.num_proc > 1) {
            num_proc_per_var_E = M2_State.num_proc_E;
            num_proc_per_var_f = M2_State.num_proc_f;
            num_proc_per_var_phi = M2_State.num_proc_Phi;
            num_proc_per_var_theta = M2_State.num_proc_Theta;
            num_proc_per_var_Triangle_gam1_gam2 = M2_State.num_proc_Triangle_gam1_gam2;
        } else {
            num_proc_per_var_E = 1;
            num_proc_per_var_f = 1;
            num_proc_per_var_phi = 1;
            num_proc_per_var_theta = 1;
            num_proc_per_var_Triangle_gam1_gam2 = 1;
        }
    
        // Determine id for each variable based on the current id_proc
        id_proc_E = 0;
        id_proc_f = floor(M2_State.id_proc/(num_proc_per_var_theta*num_proc_per_var_phi*num_proc_per_var_Triangle_gam1_gam2));   
        id_proc_phi = floor(M2_State.id_proc/(num_proc_per_var_theta*num_proc_per_var_Triangle_gam1_gam2)) - id_proc_f*num_proc_per_var_phi;
        id_proc_theta = floor(M2_State.id_proc/num_proc_per_var_Triangle_gam1_gam2) - (id_proc_f*num_proc_per_var_phi+id_proc_phi)*num_proc_per_var_theta;
        id_proc_Triangle_gam1_gam2 = M2_State.id_proc - ((id_proc_f*num_proc_per_var_phi+id_proc_phi)*num_proc_per_var_theta+id_proc_theta)*num_proc_per_var_Triangle_gam1_gam2;
    } else {
        cout << "Invalid value for Problem_Type" << endl;
        exit(0);
    }
    
    MPI_proc_parameters.num_proc_per_var_E = num_proc_per_var_E;
    MPI_proc_parameters.num_proc_per_var_f = num_proc_per_var_f;
    MPI_proc_parameters.num_proc_per_var_phi = num_proc_per_var_phi;
    MPI_proc_parameters.num_proc_per_var_theta = num_proc_per_var_theta;
    MPI_proc_parameters.num_proc_per_var_Triangle_gam1_gam2 = num_proc_per_var_Triangle_gam1_gam2;
    
    MPI_proc_parameters.id_proc_E = id_proc_E;
    MPI_proc_parameters.id_proc_f = id_proc_f;
    MPI_proc_parameters.id_proc_phi = id_proc_phi;
    MPI_proc_parameters.id_proc_theta = id_proc_theta;
    MPI_proc_parameters.id_proc_Triangle_gam1_gam2 = id_proc_Triangle_gam1_gam2;
    
    cout << "num_proc_per_var_E = " << num_proc_per_var_E << "  " << "num_proc_per_var_f = " << num_proc_per_var_f << "  " << "num_proc_per_var_phi = " << num_proc_per_var_phi << "  " << "num_proc_per_var_theta = " << num_proc_per_var_theta << "  " << "num_proc_per_var_Triangle_gam1_gam2 = " << num_proc_per_var_Triangle_gam1_gam2 << endl;
    
//         if (id_proc == 1) {
        cout << "id_proc = " << M2_State.id_proc << "  " << "id_proc_E = " << id_proc_E << "  " << "id_proc_f = " << id_proc_f << "  " << "id_proc_phi = " << id_proc_phi << "  " << "id_proc_theta = " << id_proc_theta << "  " << "id_proc_Triangle_gam1_gam2 = " << id_proc_Triangle_gam1_gam2 << endl;
        cout << "num_proc_per_var_E = " << num_proc_per_var_E << "  " << "num_proc_per_var_f = " << num_proc_per_var_f << "  " << "num_proc_per_var_phi = " << num_proc_per_var_phi << "  " << "num_proc_per_var_theta = " << num_proc_per_var_theta << "  " << "num_proc_per_var_Triangle_gam1_gam2 = " << num_proc_per_var_Triangle_gam1_gam2 << endl;
//     }
    
    if (num_points->E % num_proc_per_var_E != 0) {
        if (M2_State.id_proc == 0) {
            cout << "The number num_points->E = " << num_points->E << "is not divisible by num_proc_per_var_E = " << num_proc_per_var_E << "the quotient is " << num_points->E/num_proc_per_var_E << endl;
        }
        exit(0);
    }
    
    if (num_points->f % num_proc_per_var_f != 0) {
        if (M2_State.id_proc == 0) {
            cout << "The number num_points->f = " << num_points->f << "is not divisible by num_proc_per_var_f = " << num_proc_per_var_f << "the quotient is " << num_points->f/num_proc_per_var_f << endl;
        }
        exit(0);
    }
    
    if (num_points->phi % num_proc_per_var_phi != 0) {
        if (M2_State.id_proc == 0) {
            cout << "The number num_points->phi = " << num_points->phi << "is not divisible by num_proc_per_var_phi = " << num_proc_per_var_phi << "the quotient is " << num_points->phi/num_proc_per_var_phi << endl;
        }
        exit(0);
    }
    
    if (num_points->theta % num_proc_per_var_theta != 0) {
        if (M2_State.id_proc == 0) {
            cout << "The number num_points->theta = " << num_points->theta << "is not divisible by num_proc_per_var_theta = " << num_proc_per_var_theta << "the quotient is " << num_points->theta/num_proc_per_var_theta << endl;
        }
        exit(0);
    }
    
    if (num_points->N_Points_Triangle_gam1_gam2 % num_proc_per_var_Triangle_gam1_gam2 != 0) {
        if (M2_State.id_proc == 0) {
            cout << "The number N_Points_Triangle_gam1_gam2 = " << num_points->N_Points_Triangle_gam1_gam2 << "is not divisible by num_proc_per_var_Triangle_gam1_gam2 = " << num_proc_per_var_Triangle_gam1_gam2 << "the quotient is " << num_points->N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2 << endl;
        }
        exit(0);
    }
    
    size_rec = Compute_Npts_Total(num_points);
    size_rec /= num_proc_per_var_E*num_proc_per_var_f*num_proc_per_var_theta*num_proc_per_var_phi*num_proc_per_var_Triangle_gam1_gam2;
    
    MPI_proc_parameters.size_rec = size_rec;
}

void Compute_MPI_Processes_Max_Min_Indexes(const MPI_proc_params &MPI_proc_parameters, const MN_Var_Num_Points *num_points, int &id_E_min, int &id_E_max, int &id_f_min, int &id_f_max, int &id_phi_min, int &id_phi_max, int &id_theta_min, int &id_theta_max, int &id_Triangle_gam1_gam2_min, int &id_Triangle_gam1_gam2_max) {
    int num_proc_per_var_E, num_proc_per_var_f, num_proc_per_var_phi, num_proc_per_var_theta, num_proc_per_var_Triangle_gam1_gam2;
    int id_proc_E, id_proc_f, id_proc_phi, id_proc_theta, id_proc_Triangle_gam1_gam2;
    
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
    
    id_E_min = id_proc_E*num_points->E/num_proc_per_var_E;
    id_E_max = (id_proc_E+1)*num_points->E/num_proc_per_var_E;
    
    id_f_min = id_proc_f*num_points->f/num_proc_per_var_f;
    id_f_max = (id_proc_f+1)*num_points->f/num_proc_per_var_f;
    
    id_phi_min = id_proc_phi*num_points->phi/num_proc_per_var_phi;
    id_phi_max = (id_proc_phi+1)*num_points->phi/num_proc_per_var_phi;
    
    id_theta_min = id_proc_theta*num_points->theta/num_proc_per_var_theta;
    id_theta_max = (id_proc_theta+1)*num_points->theta/num_proc_per_var_theta;
    
    id_Triangle_gam1_gam2_min = id_proc_Triangle_gam1_gam2*num_points->N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2;
    id_Triangle_gam1_gam2_max = (id_proc_Triangle_gam1_gam2+1)*num_points->N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2;
}

void Create_MPI_Data_Type_rec_N3(MPI_Datatype &rec_N3_type) {
    // Create the rec_N3 datatype for MPI
    int lengths[1] = { 98 };
    const MPI_Aint displacements[1] = { 0 };
    MPI_Datatype types[1] = { MPI_LONG_DOUBLE };
    MPI_Type_create_struct(1, lengths, displacements, types, &rec_N3_type);
    MPI_Type_commit(&rec_N3_type);
}

void Compute_Id_gam1_gam2_Triangle(int &i_gam1, int &i_gam2, const int &i_Triangle_gam1_gam2, const MN_Var_Num_Points *num_points) {
    int id_iter_Triangle;
    int flag_break = 0;
    
    id_iter_Triangle = 0;
    for (int i_gam1_iter = 0; i_gam1_iter < num_points->gam1; i_gam1_iter++) {
        for (int i_gam2_iter = 0; i_gam2_iter < num_points->gam1 - i_gam1_iter; i_gam2_iter++) {
            // cout << "id_iter_Triangle = " << id_iter_Triangle << endl;
            id_iter_Triangle++;
            if (i_Triangle_gam1_gam2 == id_iter_Triangle - 1) {
                // cout << "AAAAAAAAAAAAAAAAAAAA" << endl;
                i_gam1 = i_gam1_iter;
                i_gam2 = i_gam2_iter;
                flag_break = 1;
                break;
            }
        }
        if (flag_break) {
            break;
        }
    }
    
    // cout << "i_Triangle_gam1_gam2 = " << i_Triangle_gam1_gam2 << "  " << "id_iter_Triangle = " << id_iter_Triangle-1 << endl;
}
