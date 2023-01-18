#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#include "../Non_Gray_Optimization/NG_MN_Model_3D_OPTIM.h"
#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED

#include <mpi.h>

int main (int argc, char *argv[]){
    int Optim_Flag, Dimension, Problem_Type;
    char path1_out[256], prefix[256], extension[256];
    MN_Var_Num_Points Num_points;
    int Node_Distribution_E, Node_Distribution_f, Node_Distribution_Phi_Theta, Node_Distribution_gam1;
    long double Length_Scale_Mobius = 0.0;
    int id_max = 1;
    M2_State_Param M2_State(9);
    
    strcpy(path1_out, getenv(PATHVAR));
    strcat(path1_out, "/M2_Model/Gray_M2_Model/Uniform_N3_Gray_M2_3D");
    
    sprintf(extension, "_%.6d", 0);
   
    strcat(extension, ".dat");
    strcat(path1_out, extension);
    
    fstream in1_out;
    in1_out.open(path1_out, ios::out|ios::binary);
    
    if (!in1_out) {
        cout << "Uniform_N3_Non_Gray_M2_3D_BE.dat could not be accessed!" << endl;
    }
    
    int ierr = 0, p = 0, id_proc = 0;
    //     Initialize MPI.
    ierr = MPI_Init ( &argc, &argv );
    //     Get the number of processes.
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );
    //     Get the ID of this process.
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id_proc );
    
    if (id_proc == 0) {
        cout << "**************************** Starting Calculations on all Processors ******************************" << endl;
    }
    
    Node_Distribution_E = UNIFORM_DISTRIBUTION; //CHEBYSHEV_SECOND_KIND_DISTRIBUTION;
    Node_Distribution_f = UNIFORM_DISTRIBUTION; //CHEBYSHEV_SECOND_KIND_DISTRIBUTION;
    Node_Distribution_Phi_Theta = UNIFORM_DISTRIBUTION;
    Node_Distribution_gam1 = UNIFORM_DISTRIBUTION; //CHEBYSHEV_SECOND_KIND_DISTRIBUTION;
    
    Problem_Type = GRAY;
    Dimension = 3;
    Num_points.E = 1; //50;
    Num_points.f = 10; //50;
    Num_points.phi = 10;
    Num_points.theta = 5; //50;
    Num_points.gam1 = 10; //50;
    
    M2_State.Dimension = Dimension;
    M2_State.Problem_Type = Problem_Type;
    M2_State.id_proc = id_proc;
    M2_State.num_proc = p;
    M2_State.Node_Dist_E = Node_Distribution_E;
    M2_State.Node_Dist_f = Node_Distribution_f;
    M2_State.Node_Dist_Phi_Theta = Node_Distribution_Phi_Theta;
    M2_State.Node_Dist_gam1 = Node_Distribution_gam1;
    M2_State.Triangle_Domain = FULL_TRIANGLE;
    M2_State.display = true; //false; //true;
    
    if (in1_out.good()) {
        cout << "Number of Points ..........................." << "N_Points_E = " << Num_points.E << "     "  << "N_Points_f = " << Num_points.f << "     "  << "N_Points_phi = " << Num_points.phi << "     "  << "N_Points_theta = " << Num_points.theta << "     "  << "N_Points_gam1 = " << Num_points.gam1 << endl;
        
        Optim_Flag = Optimization_Algorithm(M2_State, &Num_points, in1_out, 1, &Length_Scale_Mobius);
        
        in1_out.close();
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    
    // Termminate MPI
    MPI_Finalize(); // The code was not working when I did not include this
    //  Terminate

    if (id_proc == 0) {
        cout << "**************************** Calculations Terminated on all Processors ******************************" << endl;
    }
    
    return Optim_Flag;
}
