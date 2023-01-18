#ifndef _M2_CLOSURE_FIT_EIGENSTRUCTURE_H_INCLUDED
#include "M2_Closure_Fit_Eigenstructure.h"
#endif // _M2_CLOSURE_FIT_EIGENSTRUCTURE_H_INCLUDED

extern "C"
{
void F77NAME(dgeev) (char *, char *, int *, double *, int *,
                     double *, double *, double *, int *, double *, int *,
                     double *, int *, int *);
}

void Eigenvalues_overwrite( long double *dFdU, long double *REAL, long double *IMAG) {
    char JOBVL('N'), JOBVR('N');
    int N(6);
    int LDA(N), LDVL(N), LDVR(N);
    double dummy;
    double *WORK;
    int LWORK;
    int INFO;
    
    double REAL_temp[6];
    double IMAG_temp[6];
    double dFdU_temp[6*6];
    
    for (int i = 0; i < 6*6; i++) {
        dFdU_temp[i] = double(dFdU[i]);
    }
    
    LWORK = 10*N;   //how big should this thing be?
    WORK = new double[LWORK];
    
    // dgeev only works with double and not long double
    
    /* Call Fortran subroutine */
    F77NAME(dgeev)(&JOBVL, &JOBVR,
                   &N, &dFdU_temp[0], &LDA,
                   &REAL_temp[0], &IMAG_temp[0],
                   &dummy, &LDVL, &dummy, &LDVR,
                   WORK, &LWORK, &INFO);
    
    assert(INFO==0);
    
    for (int i = 0; i < N; i++) {
        REAL[i] = REAL_temp[i];
        IMAG[i] = IMAG_temp[i];
    }
    
    // cout << "INFO = " << INFO << endl;
    
    delete [] WORK; WORK=NULL;
}

//***********************************************************************
// This routine computes ...
//***********************************************************************
long double Closure_RT::q_index(const int &i, const int &j, const int &k) {
    long double q_temp;
    if (i == j && i == k) {
        if (i == 0) {
            q_temp = qxxx_val;
        } else if (i == 1){
            q_temp = qyyy_val;  
        } else {
            cout << "AAAAAAAAAAAAAA" << endl;
            exit(0);
        }
    } else if (i == j || i == k) {
        if (i == 0) {
            q_temp = qxxy_val;
        } else if (i == 1){
            q_temp = qxyy_val;  
        } else {
            cout << "BBBBBBBBBBBBBBB" << endl;
            exit(0);
        }
    } else if (j == k) {
        if (i == 0) {
            q_temp = qxyy_val;
        } else if (i == 1){
            q_temp = qxxy_val;  
        } else {
            cout << "CCCCCCCCCCCCC" << endl;
            exit(0);
        }
    } else {
        cout << "Invalid option !!!!!!" << endl;
        exit(0);
    }
    
    return q_temp;
}

//***********************************************************************
// This routine computes ...
//***********************************************************************
long double Closure_RT::d_q_index(const int &i, const int &j, const int &k, const int &index_U) {
    long double dq_dvar;
    if (i == j && i == k) {
        if (i == 0) {
            dq_dvar = dqxxx_dU[index_U];
        } else if (i == 1){
            dq_dvar = dqyyy_dU[index_U];   
        }
    } else if (i == j || i == k) {
        if (i == 0) {
            dq_dvar = dqxxy_dU[index_U];
        } else if (i == 1){
            dq_dvar = dqxyy_dU[index_U];   
        }
    } else if (j == k) {
        if (i == 0) {
            dq_dvar = dqxyy_dU[index_U];
        } else if (i == 1){
            dq_dvar = dqxxy_dU[index_U];   
        }
    } else {
        cout << "Invalid option !!!!!!" << endl;
        exit(0);
    }
    
    return dq_dvar;
}

long double Closure_RT::q_RT_index(const int &i, const int &j, const int &k) {
    long double q_temp;
    
//     if (flag_compute_closing_fluxes) {
//         cout << "qxxx_RT_val, qxxy_RT_val, qxyy_RT_val, and qyyy_RT_val not been computed yet !!!!!!" << endl;
//         exit(0);
//     }
    
    // cout << "Make sure you compute the closing properly!!!" << endl;
    
    if (i == j && i == k) {
        if (i == 0) {
            q_temp = qxxx_RT_val;
        } else if (i == 1){
            q_temp = qyyy_RT_val;   
        } else {
            cout << "AAAAAAAAAAAAA" << endl;
            exit(0);
        }
    } else if (i == j || i == k) {
        if (i == 0) {
            q_temp = qxxy_RT_val;
        } else if (i == 1){
            q_temp = qxyy_RT_val;   
        } else {
            cout << "BBBBBBBBBBBBB" << endl;
            exit(0);
        }
    } else if (j == k) {
        if (i == 0) {
            q_temp = qxyy_RT_val;
        } else if (i == 1){
            q_temp = qxxy_RT_val;   
        } else {
            cout << "CCCCCCCCCCCCC" << endl;
            exit(0);
        }
    } else {
        cout << "Invalid option here in q_RT_index !!!!!!!!!!!!!!!" << endl;
        exit(0);
    }
    
    return q_temp;
}

long double Closure_RT::d_q_RT_index(const int &i, const int &j, const int &k, const int &index_U) {
    long double dq_RT_dvar;
    if (i == j && i == k) {
        if (i == 0) {
            dq_RT_dvar = dqxxx_RT[index_U];
        } else if (i == 1){
            dq_RT_dvar = dqyyy_RT[index_U];   
        } else {
            cout << "AAAAAAAAAAAAAAA" << endl;
            exit(0);
        }
    } else if (i == j || i == k) {
        if (i == 0) {
            dq_RT_dvar = dqxxy_RT[index_U];
        } else if (i == 1){
            dq_RT_dvar = dqxyy_RT[index_U]; 
        } else {
            cout << "BBBBBBBBBBBBBB" << endl;
            exit(0);
        }
    } else if (j == k) {
        if (i == 0) {
            dq_RT_dvar = dqxyy_RT[index_U];
        } else if (i == 1){
            dq_RT_dvar = dqxxy_RT[index_U]; 
        } else {
            cout << "CCCCCCCCCCCCC" << endl;
            exit(0);
        }
    } else {
        cout << "Invalid option here in d_q_RT_index !!!!!!!!!!!!!!!" << endl;
        exit(0);
    }
    return dq_RT_dvar;
}
