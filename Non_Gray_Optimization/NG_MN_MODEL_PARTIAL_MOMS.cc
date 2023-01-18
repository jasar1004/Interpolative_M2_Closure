#include "NG_MN_Model_3D_OPTIM.h"

void Partial_Moments_M2_GRAY_1D(long double *PARTIAL_MOMS, const long double *x, int &NVARS);
 
struct rec_Partial_Moments {
 long double E, f, gam1, gam2, e_plus, fx_plus, pxx_plus, pyy_plus;
};
 
struct rec_Partial_Moments_gam2_0 {
 long double E, f, gam1, e_plus, fx_plus, pxx_plus;
};

// Function to write a record at a specific position
ostream& write(ostream& ios, const rec_Partial_Moments& r, int &pos) {
 ios.clear(); // clear any errors
 ios.seekg(sizeof(rec_Partial_Moments) * pos); // move to record's position
 ios.write(reinterpret_cast<const char*>(&r), sizeof(r)); //write the record to file
 
 return ios; // return the stream (for easy error detection/chaining)
}

// Function to read a record at a specific position
iostream& read(iostream& ios, rec_Partial_Moments_gam2_0 &r, int &pos) {
 ios.clear(); // clear any errors
 ios.seekg(sizeof(rec_Partial_Moments_gam2_0) * pos); // move to record's position
 ios.read(reinterpret_cast<char *>(&r), sizeof(r)); //read record from file
 return ios; // return the stream (for easy error detection/chaining)
}

int main() {
    int num_points, n_procs;
    long double E, f_test, gam1_test, gam2_test, Increment_gam2;
    long double a_temp, b_temp, c_temp, d_temp;
    char path11_out[256], extension[256];
    const char *path12_out, *path11_in;
    long double *PARTIAL_MOMS;
    rec_Partial_Moments r_PM;
    rec_Partial_Moments_gam2_0 r_PM_gam_0;
    int index_gam2_0;
    int p = 1;
    
    NVARS = 4;
    num_points = 50;
    n_procs = 1; //8;
    
    long double *pk_final;
    pk_final = new long double[NVARS*NVARS];
    
    for (int i = 0; i < NVARS; i++) {
        for (int j = 0; j < NVARS; j++) {
            if (i == j) {
                pk_final[i*NVARS + j] = 1.0;
            } else {500.0  500.0  500.0  500.0
                pk_final[i*NVARS + j] = 0.0;
            }
        }
    }
    
    path12_out = "./Partial_Moms_M2_1D.dat";
    path11_in = "./Partial_Moms_M2_1D_gam2_0.dat";
    
    fstream in11_in;
    in11_in.open(path11_in, ios::in|ios::binary);
    
    if (!in11_in) {
        cout << "Partial_Moms_M2_1D_gam2_0.dat could not be accessed!" << endl;
    }
    
    fstream in12_out;
    in12_out.open(path12_out, ios::out|ios::binary);
    
    if (!in12_out) {
        cout << "Partial_Moms_M2_1D.dat could not be accessed!" << endl;
    }
    
    long double x[NVARS];
    PARTIAL_MOMS = new long double[NVARS];
    
    for (int id_vals = 0; id_vals < n_procs; id_vals++) {
        strcpy(path11_out, "./Lag_M2_3D_eigen");
        sprintf(extension, "_%.4d", id_vals);
        strcat(extension, ".dat");
        strcat(path11_out, extension);
        
        ifstream in11_out;
        in11_out.open(path11_out);
        
        if (!in11_out) {
            cout << "Lag_M2_3D_Flux_1D.dat could not be accessed on proc id = " << id_vals << endl;
        }
        
        if (in11_out.good()) {
            E = 1.0;
            
            for (int i = 0; i < 2; k++) {
                for (int j = 0; j < num_points; j++) {
                    f_test = f_skewed(j, num_points);
                    for (int k = 0; k < num_points; k++) {
                        gam1_test = gam1_skewed(k, num_points);
                        Increment_gam2 = ((1.0 - gam1_test) - 0.0)/(num_points-1);
                        for (int l = 0; l < num_points; l++) {
                            gam2_test = l*Increment_gam2;
                        
                            r_PM.E = E;
                            r_PM.f = f_test;
                            r_PM.gam1 = gam1_test;
                            r_PM.gam2 = gam2_test;
                        
                            if (fabs(f_test) < 1.0 - TOLER) {
                                if (gam1_test > TOLER && (1.0 - gam1_test) > TOLER) {
                                    if (gam2_test > TOLER && (1.0 - gam1_test - gam2_test) > TOLER) {
                                        if (i == 0) {
                                            in11_out >> setprecision(12) >> a_temp >> b_temp >> c_temp >> d_temp >> x[0] >> x[1] >> x[2] >> x[3];
                                        }
                                        
                                        if (i == 1) {
                                            x[1] = -x[1];
                                        }
                                        
                                        //                             cout << setprecision(12) << "x[0] = " << x[0] << "    " << "x[1] = " << x[1] << "    " << "x[2] = " << x[2] << endl;
                                        
                                        Partial_Moments_M2_GRAY_1D(PARTIAL_MOMS, x, NVARS);
                                    } else {
                                        index_gam2_0 = j*num_points+k;
                                        if (gam2_test < TOLER) {
                                            read(in11_in, r_PM_gam_0, index_gam2_0);
                                            PARTIAL_MOMS[0] = r_PM_gam_0.e_plus;
                                            PARTIAL_MOMS[1] = r_PM_gam_0.fx_plus;
                                            PARTIAL_MOMS[2] = r_PM_gam_0.pxx_plus;
                                            PARTIAL_MOMS[3] = 0.5*(r_PM_gam_0.e_plus - r_PM_gam_0.pxx_plus);
                                        } else {
                                            read(in11_in, r_PM_gam_0, index_gam2_0);
                                            PARTIAL_MOMS[0] = r_PM_gam_0.e_plus;
                                            PARTIAL_MOMS[1] = r_PM_gam_0.fx_plus;
                                            PARTIAL_MOMS[2] = r_PM_gam_0.pxx_plus;
                                            PARTIAL_MOMS[3] = 0.5*(r_PM_gam_0.e_plus - r_PM_gam_0.pxx_plus);
                                        }
                                    }
                                } else {
                                    if (gam1_test < TOLER) {
                                        PARTIAL_MOMS[0] = E;
                                        PARTIAL_MOMS[1] = f_test;
                                        PARTIAL_MOMS[2] = pow(f_test, 2);
                                        PARTIAL_MOMS[3] = pow(f_test, 2);
                                    } else {
                                        PARTIAL_MOMS[0] = E*(1.0+f_test)/2.0;
                                        PARTIAL_MOMS[1] = E*(1.0+f_test)/2.0;
                                        PARTIAL_MOMS[2] = E*(1.0+f_test)/2.0;
                                        PARTIAL_MOMS[3] = E*(1.0+f_test)/2.0;
                                    }   
                                }
                            } else {
                                PARTIAL_MOMS[0] = E;
                                PARTIAL_MOMS[1] = f_test;
                                PARTIAL_MOMS[2] = pow(f_test, 2); // fix this
                                PARTIAL_MOMS[3] = pow(f_test, 2); // fix this
                            }
                            
                            r_PM.e_plus = PARTIAL_MOMS[0];
                            r_PM.fx_plus = PARTIAL_MOMS[1];
                            r_PM.pxx_plus = PARTIAL_MOMS[2];
                            r_PM.pyy_plus = PARTIAL_MOMS[3];
                            
                            if (i == 1) {
                                index = ((j - num_points - 1)*num_points + k)*num_points + l;
                            } else {
                                index = (j*num_points + k)*num_points + l;
                            }
                            
                            write(in12_out, r_PM, index);
                        }
                    }
                }
            }
            in11_out.close();
        }
    }
    in11_in.close();
    in12_out.close();
//     delete[] pk_final;
}

void Partial_Numerical_Cubature(unsigned fdim, integrand f, unsigned dim, 
                        unsigned maxEval, long double *val, void *fdata) {
    
    long double reqRelError, reqAbsError;
    long double *err;
    int maxevals = 1e5;
    long double xmin[2] = {0.0, 0.0};
    long double xmax[2] = {1.0, 2.0*PI};
    err = new long double[fdim];
    
    reqRelError = 0.0;
    reqAbsError = 1e-12;
     
    hcubature(fdim, f, fdata, dim, xmin, xmax, maxevals, reqAbsError, 
              reqRelError, ERROR_INDIVIDUAL, val, err);

//     for (int i = 0; i < fdim; i++) {
//         val[i] = 2.0*PI*val[i];
//     }
    
    delete[] err;
}
 
 long double sum_Lagrange_M2_1D(const long double *x, const long double &mu, const long double &phi) {
     long double sum;
     long double *poly;
     poly = new long double[4];
     generate_polynomials(poly, mu, phi);
     
     sum = x[0]*poly[0] + x[1]*poly[1] + x[2]*poly[2] + x[3]*poly[3];
     
     delete[] poly;
     return sum;
 }
 
int Partial_Moments(unsigned NDIM, const long double *Int_Vars, void *fdata, unsigned NFUN, long double *MOMENTS) {
     long double mu, phi;
     long double coeff_1;
     M2_State_Param *M2_State = (M2_State_Param *) fdata;
     
     mu = Int_Vars[0];
     phi = Int_Vars[1];
     coeff_1 = 1.0;
     
     MOMENTS[0] = coeff_1/(pow(sum_Lagrange_M2_1D(M2_State->x,mu,phi), 4));
     MOMENTS[1] = coeff_1*mu/(pow(sum_Lagrange_M2_1D(M2_State->x,mu,phi), 4));
     MOMENTS[2] = coeff_1*pow(mu,2)/(pow(sum_Lagrange_M2_1D(M2_State->x,mu,phi), 4));
     MOMENTS[3] = coeff_1*pow(sqrt(1.0-pow(mu,2))*cos(phi),2)/(pow(sum_Lagrange_M2_1D(M2_State->x,mu,phi), 4));
     
     return 0;
 }
 
 void Partial_Moments_M2_GRAY_1D(long double *PARTIAL_MOMS, const long double *x, int &NVARS) {
    M2_State_Param M2_State;
    int NF_PM;
    NF_PM = 4;
    
    M2_State.x = x;
    
    Partial_Numerical_Cubature(NF_PM, Partial_Moments, 2, 0, PARTIAL_MOMS, &M2_State);
    
//     if (id == 0) {
	    cout << "E_plus = " << PARTIAL_MOMS[0] << "    " << "fx_plus = " << PARTIAL_MOMS[1] << "    " << "pxx_plus = " << PARTIAL_MOMS[2] << "    " << "pyy_plus = " << PARTIAL_MOMS[3] << endl;
//     }
 }
