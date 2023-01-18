#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#include "./NG_MN_Model_3D_OPTIM.h"
#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED

#include <mpi.h>

int DIRECTIONS_COSINES(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata);
void DIRECTIONS_COSINES_TEST(M2_State_Param &M2_State);

int main (int argc, char *argv[]){
    M2_State_Param M2_State(9);
    int id_refinement = 2;
    long double phi_start, phi_end;
    long double mu_start, mu_end;
    
    M2_State.N1_1 = 0.4847;
    M2_State.N1_2 = 0.839525;
    M2_State.N1_3 = 0.0;
    
    M2_State.gamma_1 = 0.2;
    M2_State.gamma_2 = 0.2;
    
    M2_State.Order_Quad_mu = 5;
    M2_State.Allocate_Quad();
    M2_State.compute_Omegas(0.0, PI, 0.0, 2.0*PI);
    
    M2_State.set_bounds_quad_refin(id_refinement);
    
    DIRECTIONS_COSINES_TEST(M2_State);
    
    cout << endl;
    cout << endl;
    cout << "********************************************************************" << endl;
    cout << "N_peaks_mu = " << M2_State.N_peaks_mu << "  "  << "N_peaks_phi = " << M2_State.N_peaks_phi << endl;
    
    for (int i_peak_mu = 0; i_peak_mu < M2_State.N_peaks_mu; i_peak_mu++) {
        cout << "mu_peak = " << M2_State.mu_peak[i_peak_mu] << endl;
    }
    
    for (int i_peak_phi = 0; i_peak_phi < M2_State.N_peaks_phi; i_peak_phi++) {
        cout << "phi_peak = " << M2_State.phi_peak[i_peak_phi] << endl;
    }
    
    cout << "Quadrature bounds for Phi" << endl;
    for (int i_phi_quad = 0; i_phi_quad < M2_State.N_intervals_phi_quad; i_phi_quad++) {
        phi_start = M2_State.bounds_quad_phi_refin[i_phi_quad];
        phi_end = M2_State.bounds_quad_phi_refin[i_phi_quad+1];
        
        cout << "i_phi_quad = " << i_phi_quad << "  " << "N_intervals_phi_quad = " << M2_State.N_intervals_phi_quad << "  " << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << endl;
    }
    
    cout << endl;
    
    cout << "Quadrature bounds for theta" << endl;
    for (int i_mu_quad = 0; i_mu_quad < M2_State.N_intervals_mu_quad; i_mu_quad++) {
        mu_start = M2_State.bounds_quad_mu_refin[i_mu_quad];
        mu_end = M2_State.bounds_quad_mu_refin[i_mu_quad+1];
                    
        cout << "i_mu_quad = " << i_mu_quad << "  " << "N_intervals_mu_quad = " << M2_State.N_intervals_mu_quad << "  " << "mu_start = " << mu_start << "  " << "mu_end = " << mu_end << "  " << endl;
    }            
}



int DIRECTIONS_COSINES(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata) {
     long double Omega1, Omega2, Omega3;
     
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     
     Omega1 = M2_State->Omega1[Id_angle];
     Omega2 = M2_State->Omega2[Id_angle];
     Omega3 = M2_State->Omega3[Id_angle];
     
     MOMENTS[0] = pow(Omega1,2);
     MOMENTS[1] = pow(Omega2,2);
     MOMENTS[2] = pow(Omega3,2);
     
     MOMENTS[3] = pow(Omega1,3);
     MOMENTS[4] = pow(Omega2,3);
     MOMENTS[5] = pow(Omega3,3);
     
     return 0;
 }
 
 void DIRECTIONS_COSINES_TEST(M2_State_Param &M2_State) {
    int NFUN = 6;
    long double *Moments_Vals;
    Moments_Vals = new long double[NFUN];
    
    Lebedev_Quadrature_Func(Moments_Vals, NFUN, DIRECTIONS_COSINES, &M2_State);
    
    cout << "Omega1^2 = " << Moments_Vals[0] << "    " << "Omega2^2 = " << Moments_Vals[1] << "    " << "Omega3^2 = " << Moments_Vals[2] << "    " << "Omega1^3 = " << Moments_Vals[3] << "    " << "Omega2^3 = " << Moments_Vals[4] << "    " << "Omega3^3 = " << Moments_Vals[5] << "\n\n";
    
    delete[] Moments_Vals;
 }
