#include <mpi.h>

#ifndef _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED
#include "../M2_Model_Interp/N3_Non_Gray_M2_3D_RT_Cheby.h"
#endif // _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED

int main(int argc, char **argv) {
    int N_pts_Mobius_Scale;
    int N_Points_E, N_Points_f, N_Points_Theta, N_Points_Phi, N_Points_gam1;
    int N_Coeffs_E, N_Coeffs_f, N_Coeffs_Phi, N_Coeffs_Theta, Order_SH, N_Coeffs_gam1;
    char path_N3_ijk_out[256];
    ofstream out_N3_ijk;
    
    // Setup static variables for the N3_Non_Gray_M2_3D_RT_Cheby data structure
    N3_Non_Gray_M2_3D_RT_Cheby :: Problem_Type = GRAY;
    
    // Setup MPI
    int num_proc, id_proc, ierr, done_already;
    MPI_Initialized(&done_already);
    if (!done_already) {
        //     Initialize MPI.
        ierr = MPI_Init ( &argc, &argv );
        //     Get the number of processes.
        ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_proc );
        //     Get the ID of this process.
        ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id_proc );
    }
    
    // Setup paths
    strcpy(path_N3_ijk_out, getenv(PATHVAR));
    strcat(path_N3_ijk_out, "/M2_Model/Gray_M2_Model/Coefficients_Fit_N3_ijk_Gray_BE.dat");
                                
    out_N3_ijk.open(path_N3_ijk_out);    
    if (!out_N3_ijk) {
        cout << "Coefficients_Fit_N3_ijk_Gray_BE.dat output file could not be accessed" << endl;
    }
    
    N_pts_Mobius_Scale = 1;
    N_Points_E = 1;
    N_Points_f = 10;
    N_Points_Phi = 10;
    N_Points_Theta = 5;
    N_Points_gam1 = 10;
    
    // Set up data for M2 closure in Hyperbolic limit
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_Unif_HL(1, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, 1);
    // Set up data for M2 closure in Logarithmic limit
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_Unif_LL(1, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, 1);
    
    // Set up data for M2 closure in Bose-Einstein regime
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_Unif(N_Points_E, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mobius_Scale);
    
    // Setup id_proc and num_proc for the data structure used for the 
    // polynomial interpolation procedure
    N3_M2_3D_RT_Unif.id_proc = id_proc;
    N3_M2_3D_RT_Unif.num_proc = num_proc;
    
    N3_M2_3D_RT_Unif.num_proc_E = 1;
    N3_M2_3D_RT_Unif.num_proc_f = 1;
    N3_M2_3D_RT_Unif.num_proc_Phi = 10;
    N3_M2_3D_RT_Unif.num_proc_Theta = 5;
    N3_M2_3D_RT_Unif.num_proc_Triangle_gam1_gam2 = 2;
    
    // Setup parameters for MPI-based parallel calculations for N3_M2_3D_RT_Unif
    N3_M2_3D_RT_Unif.Create_MPI_Data_Type_rec_N3();
    N3_M2_3D_RT_Unif.Setup_MPI_Processes();
    
    N3_M2_3D_RT_Unif.OpenInputFile("Gray_M2_Model/Uniform_N3_Gray_M2_3D");
    N3_M2_3D_RT_Unif.ReadInputData();
    N3_M2_3D_RT_Unif.CloseInputFile();
    // N3_M2_3D_RT_Unif.SetupInterpolant_Values_BE();
    // N3_M2_3D_RT_Unif.SetupInterpolant_Values_BE(N3_M2_3D_RT_Unif_HL, N3_M2_3D_RT_Unif_LL);
    // N3_M2_3D_RT_Unif.Write_Data_Matlab(BOSE_EINSTEIN);
    
    // Setup number of coefficients for the purpose of the non-gray M2 closure
    // polynomial interpolation
    N_Coeffs_E = 1;
    N_Coeffs_f = 4;
    Order_SH = 8;
    N_Coeffs_Phi = 20; // 2*(Order_SH + 1);
    N_Coeffs_Theta = 10; //(Order_SH + 1);
    N_Coeffs_gam1 = 4;
    
    // Setup Data structure for the purpose of the polynomial interpolation
    // for the non-gray M2 closure
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_Cheby(N_Coeffs_E, N_Coeffs_f, N_Coeffs_Phi, N_Coeffs_Theta, N_Coeffs_gam1, 1, Order_SH);
    
    // Setup id_proc and num_proc for the data structure used for the 
    // polynomial interpolation procedure
    N3_M2_3D_RT_Cheby.id_proc = id_proc;
    N3_M2_3D_RT_Cheby.num_proc = num_proc;
    
    N3_M2_3D_RT_Cheby.num_proc_E = 1; //N_Coeffs_E;
    N3_M2_3D_RT_Cheby.num_proc_f = 1; //N_Coeffs_f;
    N3_M2_3D_RT_Cheby.num_proc_Phi = 10; //N_Coeffs_Theta;
    N3_M2_3D_RT_Cheby.num_proc_Theta = 5; //N_Coeffs_Theta;
    N3_M2_3D_RT_Cheby.num_proc_Triangle_gam1_gam2 = 2;
    
    N3_M2_3D_RT_Cheby.Polynomial_Interpolation_Gray_M2_Closure(N3_M2_3D_RT_Unif, out_N3_ijk);
    
    // Close output file for the coefficients
    if (id_proc == 0) {
        out_N3_ijk.close();
    }
    
    // Wait for all the processors running to get to this point before 
    // finalizing the MPI call
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (!done_already) {
        // Termminate MPI
        MPI_Finalize(); // The code was not working when I did not include this
        //  Terminate
    }
    
    return 0;                         
}
