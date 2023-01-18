/*******************************************************************
  File: N3_Non_Gray_M2_3D_RT_Uniform.h

  Description:  ...  

  Author:  Joachim A.R. Sarr

  Date:    May 05th, 2020
*******************************************************************/

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

class N3_Non_Gray_M2_3D_RT_Uniform {
public:
    char path_out[256], prefix[256], extension[256];
    
    record_N3 rec_N3;
    
    int N_Points_Triangle_gam1_gam2;
    
    int N_pts_Mob_Scale;
    
    int N_Points_E;
    int N_Points_f;
    int N_Points_Theta;
    int N_Points_Phi;
    int N_Points_gam1;
    
    long double *E_NON_GRAY = NULL;
    long double *N1_1_NON_GRAY = NULL;
    long double *N1_2_NON_GRAY = NULL;
    long double *N1_3_NON_GRAY = NULL;
    long double *gam1_NON_GRAY = NULL;
    long double *gam2_NON_GRAY = NULL;
    
    long double *N3_111_NON_GRAY = NULL;
    long double *N3_122_NON_GRAY = NULL;
    long double *N3_123_NON_GRAY = NULL;
    
    long double *dN3_111_NON_GRAY_dN1_1 = NULL;
    long double *dN3_122_NON_GRAY_dN1_1 = NULL;
    long double *dN3_123_NON_GRAY_dN1_1 = NULL;
    long double *dN3_123_NON_GRAY_dN1_2 = NULL;
    long double *dN3_123_NON_GRAY_dN1_3 = NULL;
    
    long double *f_N3_111_NON_GRAY = NULL;
    long double *f_N3_122_NON_GRAY = NULL;
    long double *f_N3_123_NON_GRAY = NULL;
    
    long double L_inf_Norm_N3;
    long double L2_Norm_N3;
    
    // Constructor
    N3_Non_Gray_M2_3D_RT_Uniform(const int &Num_points_E, const int &Num_points_f, const int &Num_points_Theta, const int &Num_points_Phi, const int &Num_points_gam1, const int &Num_pts_Mobius_Scale) {
        N_Points_E = Num_points_E;
        N_Points_f = Num_points_f;
        N_Points_Theta = Num_points_Theta;
        N_Points_Phi = Num_points_Phi;
        N_Points_gam1 = Num_points_gam1;
        N_pts_Mob_Scale = Num_pts_Mobius_Scale;
        
        allocate();
    }
    
    // Destructor
    ~N3_Non_Gray_M2_3D_RT_Uniform() {
        deallocate();
    }
    
    void allocate();
    
    void deallocate();
    
    void ReadInputData(char *filename);
    
    void SetupInterpolant_Values();
    
    void Write_Data_Matlab();
};

inline void N3_Non_Gray_M2_3D_RT_Uniform :: allocate() {
    N_Points_Triangle_gam1_gam2 = 0;
    for (int i_gam1 = 0; i_gam1 < N_Points_gam1; i_gam1++) {
        for (int i_gam2 = 0; i_gam2 < N_Points_gam1 - i_gam1; i_gam2++) {
            N_Points_Triangle_gam1_gam2++;
        }
    }
    
    // Deallocate first
    deallocate();
    
    // Allocate now
    int N_pts_total = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
    
    E_NON_GRAY = new long double [N_pts_total];
    N1_1_NON_GRAY = new long double [N_pts_total];
    N1_2_NON_GRAY = new long double [N_pts_total];
    N1_3_NON_GRAY = new long double [N_pts_total];
    gam1_NON_GRAY = new long double [N_pts_total];
    gam2_NON_GRAY = new long double [N_pts_total];
    
    N3_111_NON_GRAY = new long double [N_pts_total];
    N3_122_NON_GRAY = new long double [N_pts_total];
    N3_123_NON_GRAY = new long double [N_pts_total];
    
    dN3_111_NON_GRAY_dN1_1 = new long double[N_pts_total];
    dN3_122_NON_GRAY_dN1_1 = new long double[N_pts_total];
    dN3_123_NON_GRAY_dN1_1 = new long double[N_pts_total];
    dN3_123_NON_GRAY_dN1_2 = new long double[N_pts_total];
    dN3_123_NON_GRAY_dN1_3 = new long double[N_pts_total];
    
    f_N3_111_NON_GRAY = new long double [N_pts_total];
    f_N3_122_NON_GRAY = new long double [N_pts_total];
    f_N3_123_NON_GRAY = new long double [N_pts_total];
}

inline void N3_Non_Gray_M2_3D_RT_Uniform :: deallocate() {
    if (E_NON_GRAY != NULL) {
        delete[] E_NON_GRAY; E_NON_GRAY = NULL;
    }
    if (N1_1_NON_GRAY != NULL) {
        delete[] N1_1_NON_GRAY; N1_1_NON_GRAY = NULL;
    }
    if (N1_2_NON_GRAY != NULL) {
        delete[] N1_2_NON_GRAY; N1_2_NON_GRAY = NULL;
    }
    if (N1_3_NON_GRAY != NULL) {
        delete[] N1_3_NON_GRAY; N1_3_NON_GRAY = NULL;
    }
    if (gam1_NON_GRAY != NULL) {
        delete[] gam1_NON_GRAY; gam1_NON_GRAY = NULL;
    }
    if (gam2_NON_GRAY != NULL) {
        delete[] gam2_NON_GRAY; gam2_NON_GRAY = NULL;
    }
    
    if (N3_111_NON_GRAY != NULL) {
        delete[] N3_111_NON_GRAY; N3_111_NON_GRAY = NULL;
    }
    if (N3_122_NON_GRAY != NULL) {
        delete[] N3_122_NON_GRAY; N3_122_NON_GRAY = NULL;
    }
    if (N3_123_NON_GRAY != NULL) {
        delete[] N3_123_NON_GRAY; N3_123_NON_GRAY = NULL;
    }
    
    if (dN3_111_NON_GRAY_dN1_1 != NULL) {
        delete[] dN3_111_NON_GRAY_dN1_1; dN3_111_NON_GRAY_dN1_1 = NULL;
    }
    if (dN3_122_NON_GRAY_dN1_1 != NULL) {
        delete[] dN3_122_NON_GRAY_dN1_1; dN3_122_NON_GRAY_dN1_1 = NULL;
    }
    
    if (dN3_123_NON_GRAY_dN1_1 != NULL) {
        delete[] dN3_123_NON_GRAY_dN1_1; dN3_123_NON_GRAY_dN1_1 = NULL;
    }
    if (dN3_123_NON_GRAY_dN1_2 != NULL) {
        delete[] dN3_123_NON_GRAY_dN1_2; dN3_123_NON_GRAY_dN1_2 = NULL;
    }
    if (dN3_123_NON_GRAY_dN1_3 != NULL) {
        delete[] dN3_123_NON_GRAY_dN1_3; dN3_123_NON_GRAY_dN1_3 = NULL;
    }
    
    if (f_N3_111_NON_GRAY != NULL) {
        delete[] f_N3_111_NON_GRAY; f_N3_111_NON_GRAY = NULL;
    }
    if (f_N3_122_NON_GRAY != NULL) {
        delete[] f_N3_122_NON_GRAY; f_N3_122_NON_GRAY = NULL;
    }
    if (f_N3_123_NON_GRAY != NULL) {
        delete[] f_N3_123_NON_GRAY; f_N3_123_NON_GRAY = NULL;
    }
}

inline void N3_Non_Gray_M2_3D_RT_Uniform :: ReadInputData(char *filename) {
    int index;
    strcpy(path_out, getenv(PATHVAR));
    strcat(path_out, "/M2_Model/Non_Gray_Model/3D_gam1_gam2_gam3/");
    strcat(path_out, filename);
    sprintf(extension, "_%.6d", 0);
    strcat(extension, ".dat");
    strcat(path_out, extension);
        
    fstream in_out;
    in_out.open(path_out, ios::in|ios::binary);
        
    if (!in_out) {
        cout << "N3_M2_3D_gam1_gam2_gam3_RT.dat could not be accessed!" << endl;    
    }
    
    if (in_out.good()) {
        for (int id_Length_Scale_Mobius = 0; id_Length_Scale_Mobius < N_pts_Mob_Scale; id_Length_Scale_Mobius++) {
            for (int index_e = 0; index_e < N_Points_E; index_e++){
                for (int i_f = 0; i_f < N_Points_f; i_f++) {
                    for (int i_phi = 0; i_phi < N_Points_Phi; i_phi++) {
                        for (int i_theta = 0; i_theta < N_Points_Theta; i_theta++) {
                            for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_gam1_gam2++) {
                                index = (id_Length_Scale_Mobius * N_Points_E + index_e)*N_Points_f + i_f;
                                index = (index * N_Points_Phi+ i_phi) * N_Points_Theta + i_theta;
                                index = index*N_Points_Triangle_gam1_gam2 + i_gam1_gam2;
                                    
                                read<record_N3>(in_out, rec_N3, index);
                                
                                // cout << "id_Length_Scale_Mobius = " << id_Length_Scale_Mobius << "   " << "ratio_I0 = " << rec_N3.ratio_I0 << "   " << "N1_1 = " << rec_N3.N1_1 << "   " << "N1_2 = " << rec_N3.N1_2 << "   " << "N1_3 = " << rec_N3.N1_3 << "   " << "gam1 = " << rec_N3.gam1 << "   " << "gam2 = " << rec_N3.gam2 << endl;
                                    
                                E_NON_GRAY[index] = rec_N3.ratio_I0;
                                N1_1_NON_GRAY[index] = rec_N3.N1_1;
                                N1_2_NON_GRAY[index] = rec_N3.N1_2;
                                N1_3_NON_GRAY[index] = rec_N3.N1_3;
                                gam1_NON_GRAY[index] = rec_N3.gam1;
                                gam2_NON_GRAY[index] = rec_N3.gam2;
                                    
                                N3_111_NON_GRAY[index] = rec_N3.N3_111;
                                N3_122_NON_GRAY[index] = rec_N3.N3_122;
                                N3_123_NON_GRAY[index] = rec_N3.N3_123;
                                
                                dN3_111_NON_GRAY_dN1_1[index] = rec_N3.dN3_111_dN1_1;
                                dN3_122_NON_GRAY_dN1_1[index] = rec_N3.dN3_122_dN1_1;
                                dN3_123_NON_GRAY_dN1_1[index] = rec_N3.dN3_123_dN1_1;
                                dN3_123_NON_GRAY_dN1_2[index] = rec_N3.dN3_123_dN1_2;
                                dN3_123_NON_GRAY_dN1_3[index] = rec_N3.dN3_123_dN1_3;
                            }
                        }
                    }
                }
            }
        }
        in_out.close();
    }
}

inline void N3_Non_Gray_M2_3D_RT_Uniform :: SetupInterpolant_Values() {
    int index;
    long double norm_f, gam3;
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        for (int i_E = 0; i_E < N_Points_E; i_E++) {
            for (int i_f = 0; i_f < N_Points_f; i_f++) {
                for (int i_Phi = 0; i_Phi < N_Points_Phi; i_Phi++) {
                    for (int i_Theta = 0; i_Theta < N_Points_Theta; i_Theta++) {
                        for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_gam1_gam2++) {
                            index = id_Mobius*N_Points_E + i_E;
                            index = ((index*N_Points_f + i_f) * N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
                            index = index*N_Points_Triangle_gam1_gam2 + i_gam1_gam2;
                            
                            norm_f = pow(N1_1_NON_GRAY[index],2) + pow(N1_2_NON_GRAY[index],2) + pow(N1_3_NON_GRAY[index],2);
                            f_N3_111_NON_GRAY[index] = (N3_111_NON_GRAY[index] - pow(N1_1_NON_GRAY[index],3))/(N1_1_NON_GRAY[index]*(1.0 - norm_f));
                            
                            f_N3_122_NON_GRAY[index] = (N3_122_NON_GRAY[index] - N1_1_NON_GRAY[index]*pow(N1_2_NON_GRAY[index],2))/(N1_1_NON_GRAY[index]*(1.0 - norm_f));
                            
                            f_N3_123_NON_GRAY[index] = N3_123_NON_GRAY[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);
                            
                            if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
                                f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dN1_1[index]/(1.0 - norm_f);
                                f_N3_122_NON_GRAY[index] = (dN3_122_NON_GRAY_dN1_1[index]-pow(N1_2_NON_GRAY[index],2))/(1.0 - norm_f);
                                f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);
                            } else if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
                                f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dN1_2[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);
                            } else if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
                                f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dN1_3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
                            }
                            
                            if (norm_f < 1.0e-8) {
                                f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dN1_1[index]*dN3_123_NON_GRAY_dN1_2[index]*dN3_123_NON_GRAY_dN1_3[index];
                            } else if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8 && fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
                                f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dN1_2[index]*dN3_123_NON_GRAY_dN1_3[index]/(N1_1_NON_GRAY[index]);
                            }
                            
                            f_N3_111_NON_GRAY[index] /= gam1_NON_GRAY[index];
                            f_N3_111_NON_GRAY[index] = (f_N3_111_NON_GRAY[index] - 1.0)/(1.0 - gam1_NON_GRAY[index]);
                            
                            f_N3_122_NON_GRAY[index] /= gam2_NON_GRAY[index];
                            f_N3_122_NON_GRAY[index] = (f_N3_122_NON_GRAY[index] - 1.0)/gam1_NON_GRAY[index];
                            
                            gam3 = 1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index];
                            f_N3_123_NON_GRAY[index] = (N3_123_NON_GRAY[index] - N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index])/(gam1_NON_GRAY[index]*gam2_NON_GRAY[index]*gam3);
                        } // end for i_gam1_gam2
                    } // end for i_Theta
                } // end for i_Phi
            } // end for i_f
        } // end for i_E
    } // end for id_Mobius
}

inline void N3_Non_Gray_M2_3D_RT_Uniform :: Write_Data_Matlab() {
    int index;
    char path_out[256];
    strcpy(path_out, getenv(PATHVAR));
    strcat(path_out, "/M2_Model/Non_Gray_Model/3D_gam1_gam2_gam3/N3_3D_RT_Datas_For_Matlab.dat");
        
    ofstream in_out;
    in_out.open(path_out);
        
    if (!in_out) {
        cout << "N3_3D_RT_Datas_For_Matlab.dat could not be accessed!" << endl;    
    }
    
    if (in_out.good()) {
        for (int index_e = 0 ; index_e < N_Points_E; index_e++){
            for (int i_f = 0; i_f < N_Points_f; i_f++) {
                for (int i_Phi = 0 ; i_Phi < N_Points_Phi; i_Phi++) {
                    for (int i_Theta = 0 ; i_Theta < N_Points_Theta; i_Theta++) {
                        for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_gam1_gam2++) {
                            index = index_e*N_Points_f + i_f;
                            index = (index * N_Points_Phi+ i_Phi) * N_Points_Theta + i_Theta;
                            index = index*N_Points_Triangle_gam1_gam2 + i_gam1_gam2;
                            
                            in_out << 1.0 << setw(18) << N1_1_NON_GRAY[index] << setw(18) << N1_2_NON_GRAY[index] << setw(18) << N1_3_NON_GRAY[index] << setw(18) << gam1_NON_GRAY[index] << setw(18) << gam2_NON_GRAY[index] << setw(18) << N3_111_NON_GRAY[index] << setw(18) << N3_122_NON_GRAY[index] << setw(18) << N3_123_NON_GRAY[index] << endl;
                        }
                    }
                }
            }
        }
        in_out.close();
    }
}
