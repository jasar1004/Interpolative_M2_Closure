/*******************************************************************
  File: N3_Non_Gray_M2_3D_RT_Cheby.cc

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
#include <new>

#ifndef _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED
#include "N3_Non_Gray_M2_3D_RT_Cheby.h"
#endif // _N3_NON_GRAY_M2_3D_RT_CHEBY_H_INCLUDED

int N3_ijk_Entry_Index[3] = {N3_111_ENTRY, N3_122_ENTRY, N3_123_ENTRY};

int N3_Non_Gray_M2_3D_RT_Cheby :: Problem_Type = GRAY;
MPI_Datatype N3_Non_Gray_M2_3D_RT_Cheby :: rec_N3_type;
int N3_Non_Gray_M2_3D_RT_Cheby :: closure_type = MOMENT_CLOSURE_M2; // MOMENT_CLOSURE_PROJECTION;
int N3_Non_Gray_M2_3D_RT_Cheby :: flag_basis_type = MONOMIAL_BASIS;

//************************************************************************
// Allocate arrays for data structure N3_Non_Gray_M2_3D_RT_Cheby
//************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: allocate() {
    // Deallocate first
    deallocate();
    
    int N_pts_total;
    
    int sum;
    N_Coeffs_SH = 0;
    for (int l_SH = 0; l_SH <= Order_SH; l_SH++) {
        for (int m_SH = -l_SH; m_SH <= l_SH; m_SH++) {
            sum  = l_SH + fabs(m_SH);
            if (m_SH >= 0 && l_SH % 2 == 0 && m_SH % 2 == 0) {
                // cout << "l_SH = " << l_SH << "   " << "m_SH = " << m_SH << "   " << "sum = " << sum << endl;
                N_Coeffs_SH++;
            }
        }
    }
    
    // Allocate now
    Array_l_SH = new int[N_Coeffs_SH];
    Array_m_SH = new int[N_Coeffs_SH];
    
    Setup_SH_Indexes();
    
    N_Points_Triangle_gam1_gam2 = 0;
    for (int i_gam1 = 0; i_gam1 < N_Points_gam1; i_gam1++) {
        for (int i_gam2 = 0; i_gam2 < N_Points_gam1 - i_gam1; i_gam2++) {
            N_Points_Triangle_gam1_gam2++;
        }
    }
    
//     cout << "N_Points_Triangle_gam1_gam2 = " << N_Points_Triangle_gam1_gam2 << endl;
    
    // Allocate now
    Coefficients_Vander_Matrix_SH = new long double[N_Coeffs_SH*N_Coeffs_SH];
    
    Coefficients_Vander_Matrix_Least_Squares_L_I0_star = new long double[N_Points_E*N_Points_E];
    
    N_pts_total = N_Points_f*N_Points_Triangle_gam1_gam2;
    Coefficients_Vander_Matrix_L_I0_star_N3_ijk = new long double[N_pts_total*N_pts_total];
    
    N_pts_total = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
    Coefficients_Vander_Matrix_NG_N3_ijk = new long double[N_pts_total*N_pts_total];
    
    N_pts_total = N_Points_f*N_Points_Theta*N_Points_Phi*N_Points_Triangle_gam1_gam2;
    Optim_Length_Scale_N3_ijk = new long double[N_pts_total];
    
    N_pts_total = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Theta*N_Points_Phi*N_Points_Triangle_gam1_gam2;
    rec_N3 = new record_N3[N_pts_total];
    
    E_NON_GRAY = new long double [N_pts_total];
    N1_1_NON_GRAY = new long double [N_pts_total];
    N1_2_NON_GRAY = new long double [N_pts_total];
    N1_3_NON_GRAY = new long double [N_pts_total];
    gam1_NON_GRAY = new long double [N_pts_total];
    gam2_NON_GRAY = new long double [N_pts_total];
    
    x_SH = new long double [N_pts_total];
    y_SH = new long double [N_pts_total];
    z_SH = new long double [N_pts_total];
    
    N3_111_NON_GRAY = new long double [N_pts_total];
    N3_112_NON_GRAY = new long double [N_pts_total];
    N3_122_NON_GRAY = new long double [N_pts_total];
    N3_123_NON_GRAY = new long double [N_pts_total];
    N3_222_NON_GRAY = new long double [N_pts_total];
    
    f_N3_111_NON_GRAY = new long double [N_pts_total];
    f_N3_122_NON_GRAY = new long double [N_pts_total];
    f_N3_123_NON_GRAY = new long double [N_pts_total];
    
    // Derivatives for N3_111
    dN3_111_NON_GRAY_dN1_1 = new long double [N_pts_total];
    dN3_111_NON_GRAY_dN1_2 = new long double [N_pts_total];
    dN3_111_NON_GRAY_dN1_3 = new long double [N_pts_total];
    dN3_111_NON_GRAY_dmu = new long double [N_pts_total];
    dN3_111_NON_GRAY_dnorm_f = new long double [N_pts_total];
    dN3_111_NON_GRAY_dgam1 = new long double [N_pts_total];
    d2_N3_111_NON_GRAY_dnorm_f_dN1_1 = new long double [N_pts_total];
    d2_N3_111_NON_GRAY_dnorm_f_dgam1 = new long double [N_pts_total];
    d2_N3_111_NON_GRAY_dgam1_dN1_1 = new long double [N_pts_total];
    d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1 = new long double [N_pts_total];
    
    // Derivatives for N3_122
    dN3_122_NON_GRAY_dN1_1 = new long double [N_pts_total];
    dN3_122_NON_GRAY_dN1_2 = new long double [N_pts_total];
    dN3_122_NON_GRAY_dN1_3 = new long double [N_pts_total];
    dN3_122_NON_GRAY_dmu = new long double [N_pts_total];
    dN3_122_NON_GRAY_dnorm_f = new long double [N_pts_total];
    dN3_122_NON_GRAY_dgam1 = new long double [N_pts_total];
    dN3_122_NON_GRAY_dgam2 = new long double [N_pts_total];
    d2_N3_122_NON_GRAY_dgam1_dgam2 = new long double [N_pts_total];
    d2_N3_122_NON_GRAY_dnorm_f_dN1_1 = new long double [N_pts_total];
    d2_N3_122_NON_GRAY_dnorm_f_dgam1 = new long double [N_pts_total];
    d2_N3_122_NON_GRAY_dnorm_f_dgam2 = new long double [N_pts_total];
    d2_N3_122_NON_GRAY_dgam1_dN1_1 = new long double [N_pts_total];
    d2_N3_122_NON_GRAY_dgam2_dN1_1 = new long double [N_pts_total];
    d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1 = new long double [N_pts_total];
    d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1 = new long double [N_pts_total];
    d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2 = new long double [N_pts_total];
    d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1 = new long double [N_pts_total];
    d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1 = new long double [N_pts_total];
    
    // Derivatives for N3_123
    dN3_123_NON_GRAY_dN1_1 = new long double [N_pts_total];
    dN3_123_NON_GRAY_dN1_2 = new long double [N_pts_total];
    dN3_123_NON_GRAY_dN1_3 = new long double [N_pts_total];
    dN3_123_NON_GRAY_dmu = new long double [N_pts_total];
    dN3_123_NON_GRAY_dgam1 = new long double [N_pts_total];
    dN3_123_NON_GRAY_dgam2 = new long double [N_pts_total];
    dN3_123_NON_GRAY_dgam3 = new long double [N_pts_total];
    
    d2_N3_123_NON_GRAY_dN1_1_dN1_2 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_1_dN1_3 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_2_dN1_3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3 = new long double [N_pts_total];
    
    d2_N3_123_NON_GRAY_dgam1_dN1_1 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dgam2_dN1_1 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_1_dgam3 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_2_dgam1 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_2_dgam2 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_2_dgam3 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_3_dgam1 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_3_dgam2 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dN1_3_dgam3 = new long double [N_pts_total];
    
    d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3 = new long double [N_pts_total];
    
    d2_N3_123_NON_GRAY_dgam1_dgam2 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dgam1_dgam3 = new long double [N_pts_total];
    d2_N3_123_NON_GRAY_dgam2_dgam3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3 = new long double [N_pts_total];
    d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3 = new long double [N_pts_total];
    
    d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3 = new long double [N_pts_total];
    d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3 = new long double [N_pts_total];
    d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2 = new long double [N_pts_total];
    d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3 = new long double [N_pts_total];
    d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3 = new long double [N_pts_total];
    
    N_pts_total = N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
    Coefficient_Mobius_Scale_N3_ijk = new long double[N_pts_total];
    
    Coefficient_Matrix_Fit_N3_111_HL = new long double[N_pts_total];
    Coefficient_Matrix_Fit_N3_111_LL = new long double[N_pts_total];
    
    Coefficient_Matrix_Fit_N3_122_HL = new long double[N_pts_total];
    Coefficient_Matrix_Fit_N3_122_LL = new long double[N_pts_total];
    
    Coefficient_Matrix_Fit_N3_123_HL = new long double[N_pts_total];
    Coefficient_Matrix_Fit_N3_123_LL = new long double[N_pts_total];
    
    N_pts_total = N_Points_E*N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
    
    Coeff_Lebedev_Quadrature_N3_111 = new long double [N_pts_total];
    Coeff_Lebedev_Quadrature_N3_122 = new long double [N_pts_total];
    Coeff_Lebedev_Quadrature_N3_123 = new long double [N_pts_total];
    
    Coefficient_Matrix_Fit_N3_111 = new long double [N_pts_total];
    Coefficient_Matrix_Fit_N3_122 = new long double [N_pts_total];
    Coefficient_Matrix_Fit_N3_123 = new long double [N_pts_total];
    
    
    Coefficient_Matrix_Fit_N3_111_Cheby_Basis = new long double [N_pts_total];
    Coefficient_Matrix_Fit_N3_122_Cheby_Basis = new long double [N_pts_total];
    Coefficient_Matrix_Fit_N3_123_Cheby_Basis = new long double [N_pts_total];
}                

//************************************************************************
// Deallocate arrays for data structure N3_Non_Gray_M2_3D_RT_Cheby
//************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: deallocate() {
    if (Optim_Length_Scale_N3_ijk != NULL) {
        delete[] Optim_Length_Scale_N3_ijk; Optim_Length_Scale_N3_ijk = NULL;
    }
    if (Coefficient_Mobius_Scale_N3_ijk != NULL) {
        delete[] Coefficient_Mobius_Scale_N3_ijk; Coefficient_Mobius_Scale_N3_ijk = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_111_HL!= NULL) {
        delete[] Coefficient_Matrix_Fit_N3_111_HL; Coefficient_Matrix_Fit_N3_111_HL = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_111_LL!= NULL) {
        delete[] Coefficient_Matrix_Fit_N3_111_LL; Coefficient_Matrix_Fit_N3_111_LL = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_122_HL!= NULL) {
        delete[] Coefficient_Matrix_Fit_N3_122_HL; Coefficient_Matrix_Fit_N3_122_HL = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_122_LL!= NULL) {
        delete[] Coefficient_Matrix_Fit_N3_122_LL; Coefficient_Matrix_Fit_N3_122_LL = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_123_HL!= NULL) {
        delete[] Coefficient_Matrix_Fit_N3_123_HL; Coefficient_Matrix_Fit_N3_123_HL = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_123_LL!= NULL) {
        delete[] Coefficient_Matrix_Fit_N3_123_LL; Coefficient_Matrix_Fit_N3_123_LL = NULL;
    }
    
    if (Coefficients_Vander_Matrix_SH != NULL) {
        delete[] Coefficients_Vander_Matrix_SH; Coefficients_Vander_Matrix_SH = NULL;
    }
    
    if (Coefficients_Vander_Matrix_Least_Squares_L_I0_star != NULL) {
        delete[] Coefficients_Vander_Matrix_Least_Squares_L_I0_star; Coefficients_Vander_Matrix_Least_Squares_L_I0_star = NULL;
    }
    
    if (Coefficients_Vander_Matrix_L_I0_star_N3_ijk != NULL) {
        delete[] Coefficients_Vander_Matrix_L_I0_star_N3_ijk; Coefficients_Vander_Matrix_L_I0_star_N3_ijk = NULL;
    }
    if (Coefficients_Vander_Matrix_NG_N3_ijk != NULL) {
        delete[] Coefficients_Vander_Matrix_NG_N3_ijk; Coefficients_Vander_Matrix_NG_N3_ijk = NULL;
    }
    
    if (rec_N3 != NULL) {
        delete[] rec_N3; rec_N3 = NULL;
    }
    
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
    
    if (x_SH != NULL) {
        delete[] x_SH; x_SH = NULL;
    }
    if (y_SH != NULL) {
        delete[] y_SH; y_SH = NULL;
    }
    if (z_SH != NULL) {
        delete[] z_SH; z_SH = NULL;
    }
    
    if (N3_111_NON_GRAY != NULL) {   
        delete[] N3_111_NON_GRAY; N3_111_NON_GRAY = NULL;
    }
    if (N3_112_NON_GRAY != NULL) {
        delete[] N3_112_NON_GRAY; N3_112_NON_GRAY = NULL;
    }
    if (N3_122_NON_GRAY != NULL) {
        delete[] N3_122_NON_GRAY; N3_122_NON_GRAY = NULL;
    }
    if (N3_123_NON_GRAY != NULL) {
        delete[] N3_123_NON_GRAY; N3_123_NON_GRAY = NULL;
    }
    if (N3_222_NON_GRAY != NULL) {
        delete[] N3_222_NON_GRAY; N3_222_NON_GRAY = NULL;
    }
    
    // Derivatives for N3_111
    if (dN3_111_NON_GRAY_dN1_1 != NULL) {   
        delete[] dN3_111_NON_GRAY_dN1_1; dN3_111_NON_GRAY_dN1_1 = NULL;
    }
    if (dN3_111_NON_GRAY_dN1_2 != NULL) {   
        delete[] dN3_111_NON_GRAY_dN1_2; dN3_111_NON_GRAY_dN1_2 = NULL;
    }
    if (dN3_111_NON_GRAY_dN1_3 != NULL) {   
        delete[] dN3_111_NON_GRAY_dN1_3; dN3_111_NON_GRAY_dN1_3 = NULL;
    }
    if (dN3_111_NON_GRAY_dmu != NULL) {   
        delete[] dN3_111_NON_GRAY_dmu; dN3_111_NON_GRAY_dmu = NULL;
    }
    
    if (dN3_111_NON_GRAY_dgam1 != NULL) {   
        delete[] dN3_111_NON_GRAY_dgam1; dN3_111_NON_GRAY_dgam1 = NULL;
    }
    if (dN3_111_NON_GRAY_dnorm_f != NULL) {   
        delete[] dN3_111_NON_GRAY_dnorm_f; dN3_111_NON_GRAY_dnorm_f = NULL;
    }
    if (d2_N3_111_NON_GRAY_dnorm_f_dN1_1 != NULL) {   
        delete[] d2_N3_111_NON_GRAY_dnorm_f_dN1_1; d2_N3_111_NON_GRAY_dnorm_f_dN1_1 = NULL;
    }
    if (d2_N3_111_NON_GRAY_dnorm_f_dgam1 != NULL) {   
        delete[] d2_N3_111_NON_GRAY_dnorm_f_dgam1; d2_N3_111_NON_GRAY_dnorm_f_dgam1 = NULL;
    }
    if (d2_N3_111_NON_GRAY_dgam1_dN1_1 != NULL) {   
        delete[] d2_N3_111_NON_GRAY_dgam1_dN1_1; d2_N3_111_NON_GRAY_dgam1_dN1_1 = NULL;
    }
    if (d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1 != NULL) {   
        delete[] d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1; d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1 = NULL;
    }
    
    // Derivatives for N3_122
    if (dN3_122_NON_GRAY_dN1_1 != NULL) {
        delete[] dN3_122_NON_GRAY_dN1_1; dN3_122_NON_GRAY_dN1_1 = NULL;
    }
    if (dN3_122_NON_GRAY_dN1_2 != NULL) {
        delete[] dN3_122_NON_GRAY_dN1_2; dN3_122_NON_GRAY_dN1_2 = NULL;
    }
    if (dN3_122_NON_GRAY_dN1_3 != NULL) {
        delete[] dN3_122_NON_GRAY_dN1_3; dN3_122_NON_GRAY_dN1_3 = NULL;
    }
    if (dN3_122_NON_GRAY_dmu != NULL) {
        delete[] dN3_122_NON_GRAY_dmu; dN3_122_NON_GRAY_dmu = NULL;
    }
    if (dN3_122_NON_GRAY_dgam1 != NULL) {   
        delete[] dN3_122_NON_GRAY_dgam1; dN3_122_NON_GRAY_dgam1 = NULL;
    }
    if (dN3_122_NON_GRAY_dgam2 != NULL) {   
        delete[] dN3_122_NON_GRAY_dgam2; dN3_122_NON_GRAY_dgam2 = NULL;
    }
    if (dN3_122_NON_GRAY_dnorm_f != NULL) {   
        delete[] dN3_122_NON_GRAY_dnorm_f; dN3_122_NON_GRAY_dnorm_f = NULL;
    }
    if (d2_N3_122_NON_GRAY_dgam1_dgam2 != NULL) {   
        delete[] d2_N3_122_NON_GRAY_dgam1_dgam2; d2_N3_122_NON_GRAY_dgam1_dgam2 = NULL;
    }
    if (d2_N3_122_NON_GRAY_dnorm_f_dN1_1 != NULL) {   
        delete[] d2_N3_122_NON_GRAY_dnorm_f_dN1_1; d2_N3_122_NON_GRAY_dnorm_f_dN1_1 = NULL;
    }
    if (d2_N3_122_NON_GRAY_dnorm_f_dgam1 != NULL) {   
        delete[] d2_N3_122_NON_GRAY_dnorm_f_dgam1; d2_N3_122_NON_GRAY_dnorm_f_dgam1 = NULL;
    }
    if (d2_N3_122_NON_GRAY_dnorm_f_dgam2 != NULL) {   
        delete[] d2_N3_122_NON_GRAY_dnorm_f_dgam2; d2_N3_122_NON_GRAY_dnorm_f_dgam2 = NULL;
    }
    if (d2_N3_122_NON_GRAY_dgam1_dN1_1 != NULL) {   
        delete[] d2_N3_122_NON_GRAY_dgam1_dN1_1; d2_N3_122_NON_GRAY_dgam1_dN1_1 = NULL;
    }
    if (d2_N3_122_NON_GRAY_dgam2_dN1_1 != NULL) {   
        delete[] d2_N3_122_NON_GRAY_dgam2_dN1_1; d2_N3_122_NON_GRAY_dgam2_dN1_1 = NULL;
    }
    if (d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1 != NULL) {   
        delete[] d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1; d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1 = NULL;
    }
    if (d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1 != NULL) {   
        delete[] d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1; d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1 = NULL;
    }
    if (d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2 != NULL) {   
        delete[] d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2; d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2 = NULL;
    }
    if (d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1 != NULL) {   
        delete[] d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1; d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1 = NULL;
    }
    if (d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1 != NULL) {   
        delete[] d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1; d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1 = NULL;
    }
    
    // Derivatives for N3_123
    if (dN3_123_NON_GRAY_dN1_1 != NULL) {
        delete[] dN3_123_NON_GRAY_dN1_1; dN3_123_NON_GRAY_dN1_1 = NULL;
    }
    if (dN3_123_NON_GRAY_dN1_2 != NULL) {
        delete[] dN3_123_NON_GRAY_dN1_2; dN3_123_NON_GRAY_dN1_2 = NULL;
    }
    if (dN3_123_NON_GRAY_dN1_3 != NULL) {
        delete[] dN3_123_NON_GRAY_dN1_3; dN3_123_NON_GRAY_dN1_3 = NULL;
    }
    if (dN3_123_NON_GRAY_dmu != NULL) {
        delete[] dN3_123_NON_GRAY_dmu; dN3_123_NON_GRAY_dmu = NULL;
    }
    if (dN3_123_NON_GRAY_dgam1 != NULL) {   
        delete[] dN3_123_NON_GRAY_dgam1; dN3_123_NON_GRAY_dgam1 = NULL;
    }
    if (dN3_123_NON_GRAY_dgam2 != NULL) {   
        delete[] dN3_123_NON_GRAY_dgam2; dN3_123_NON_GRAY_dgam2 = NULL;
    }
    if (dN3_123_NON_GRAY_dgam3 != NULL) {   
        delete[] dN3_123_NON_GRAY_dgam3; dN3_123_NON_GRAY_dgam3 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_1_dN1_2 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_1_dN1_2; d2_N3_123_NON_GRAY_dN1_1_dN1_2 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_1_dN1_3 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_1_dN1_3; d2_N3_123_NON_GRAY_dN1_1_dN1_3 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_2_dN1_3 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_2_dN1_3; d2_N3_123_NON_GRAY_dN1_2_dN1_3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3; d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dgam1_dN1_1 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dgam1_dN1_1; d2_N3_123_NON_GRAY_dgam1_dN1_1 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dgam2_dN1_1 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dgam2_dN1_1; d2_N3_123_NON_GRAY_dgam2_dN1_1 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_1_dgam3 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_1_dgam3; d2_N3_123_NON_GRAY_dN1_1_dgam3 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_2_dgam1 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_2_dgam1; d2_N3_123_NON_GRAY_dN1_2_dgam1 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_2_dgam2 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_2_dgam2; d2_N3_123_NON_GRAY_dN1_2_dgam2 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_2_dgam3 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_2_dgam3; d2_N3_123_NON_GRAY_dN1_2_dgam3 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_3_dgam1 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_3_dgam1; d2_N3_123_NON_GRAY_dN1_3_dgam1 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_3_dgam2 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_3_dgam2; d2_N3_123_NON_GRAY_dN1_3_dgam2 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dN1_3_dgam3 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dN1_3_dgam3; d2_N3_123_NON_GRAY_dN1_3_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1; d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1; d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1; d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1; d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2; d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2; d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2; d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2; d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3; d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3; d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3; d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3; d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dgam1_dgam2 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dgam1_dgam2; d2_N3_123_NON_GRAY_dgam1_dgam2 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dgam1_dgam3 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dgam1_dgam3; d2_N3_123_NON_GRAY_dgam1_dgam3 = NULL;
    }
    if (d2_N3_123_NON_GRAY_dgam2_dgam3 != NULL) {   
        delete[] d2_N3_123_NON_GRAY_dgam2_dgam3; d2_N3_123_NON_GRAY_dgam2_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1; d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3; d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3; d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2; d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3; d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3; d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2; d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3; d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3 = NULL;
    }
    if (d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3 != NULL) {   
        delete[] d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3; d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2; d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3; d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3; d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2; d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3; d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3; d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2; d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3; d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3 = NULL;
    }
    if (d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3 != NULL) {   
        delete[] d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3; d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3 = NULL;
    }
    if (d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2 != NULL) {   
        delete[] d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2; d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2 = NULL;
    }
    if (d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3 != NULL) {   
        delete[] d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3; d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3 = NULL;
    }
    if (d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3 != NULL) {   
        delete[] d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3; d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3 = NULL;
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
    
    //
    if (Coefficient_Matrix_Fit_N3_111 != NULL) {
        delete[] Coefficient_Matrix_Fit_N3_111; Coefficient_Matrix_Fit_N3_111 = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_122 != NULL) {
        delete[] Coefficient_Matrix_Fit_N3_122; Coefficient_Matrix_Fit_N3_122 = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_123 != NULL) {
        delete[] Coefficient_Matrix_Fit_N3_123; Coefficient_Matrix_Fit_N3_123 = NULL;
    }
    
    //
    if (Coefficient_Matrix_Fit_N3_111_Cheby_Basis != NULL) {
        delete[] Coefficient_Matrix_Fit_N3_111_Cheby_Basis; Coefficient_Matrix_Fit_N3_111_Cheby_Basis = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_122_Cheby_Basis != NULL) {
        delete[] Coefficient_Matrix_Fit_N3_122_Cheby_Basis; Coefficient_Matrix_Fit_N3_122_Cheby_Basis = NULL;
    }
    if (Coefficient_Matrix_Fit_N3_123_Cheby_Basis != NULL) {
        delete[] Coefficient_Matrix_Fit_N3_123_Cheby_Basis; Coefficient_Matrix_Fit_N3_123_Cheby_Basis = NULL;
    }
    
    //
    if (Coeff_Lebedev_Quadrature_N3_111 != NULL) {
        delete[] Coeff_Lebedev_Quadrature_N3_111; Coeff_Lebedev_Quadrature_N3_111 = NULL;
    }
    if (Coeff_Lebedev_Quadrature_N3_122 != NULL) {
        delete[] Coeff_Lebedev_Quadrature_N3_122; Coeff_Lebedev_Quadrature_N3_122 = NULL;
    }
    if (Coeff_Lebedev_Quadrature_N3_123 != NULL) {
        delete[] Coeff_Lebedev_Quadrature_N3_123; Coeff_Lebedev_Quadrature_N3_123 = NULL;
    }
}

//******************************************************************************************
// This routines sets up the indexes for the spherical harmonics used for our interpolative-
// based approximation of the third-order closing fluxes for the gray or non-gray M2 closure
//******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Setup_SH_Indexes() {
    int sum;
    int index_SH = 0;
    for (int l_SH = 0; l_SH <= Order_SH; l_SH++) {
        for (int m_SH = -l_SH; m_SH <= l_SH; m_SH++) {
            sum  = l_SH + fabs(m_SH);
            if (m_SH >= 0 && l_SH % 2 == 0 && m_SH % 2 == 0) {
                Array_l_SH[index_SH] = l_SH;
                Array_m_SH[index_SH] = m_SH;
                index_SH++;
            }
        }
    }
    
    if (index_SH != N_Coeffs_SH) {
        cout << "Isssue with Setup_SH_Indexes; index_SH = " << index_SH << "  " << "N_Coeffs_SH = " << N_Coeffs_SH << endl;
        exit(0);
    }
}

//******************************************************************************************
// Setup parameters required for the purpose of parallel computing
//******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Setup_MPI_Processes() {
    int num_proc_per_var_E, num_proc_per_var_f, num_proc_per_var_phi, num_proc_per_var_theta, num_proc_per_var_Triangle_gam1_gam2;
    int id_proc_E, id_proc_f, id_proc_theta, id_proc_phi, id_proc_Triangle_gam1_gam2;
    int N_pts_total;
    int size_rec;
    int num_proc_temp;
    
    N_pts_total = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Theta*N_Points_Phi*N_Points_Triangle_gam1_gam2;
    
    if ((Problem_Type == NON_GRAY) // && 
        // (num_points->Maximum_Entropy_Solution_Regime != HYPERBOLIC_LIMIT) && 
        // (num_points->Maximum_Entropy_Solution_Regime != LOGARITHMIC_LIMIT)
       ) {
        // Compute number of processors per for each variable
        if (num_proc > 1) {
            num_proc_per_var_E = num_proc_E;
            num_proc_per_var_f = num_proc_f;
            num_proc_per_var_phi = num_proc_Phi;
            num_proc_per_var_theta = num_proc_Theta;
            num_proc_per_var_Triangle_gam1_gam2 = num_proc_Triangle_gam1_gam2;
        } else {
            num_proc_per_var_E = 1;
            num_proc_per_var_f = 1;
            num_proc_per_var_phi = 1;
            num_proc_per_var_theta = 1;
            num_proc_per_var_Triangle_gam1_gam2 = 1;
        }
        
        num_proc_used = num_proc_per_var_E*num_proc_per_var_f*num_proc_per_var_phi*num_proc_per_var_theta*num_proc_per_var_Triangle_gam1_gam2;
    
        // Determine id for each variable based on the current id_proc
        if (id_proc < num_proc_used) {
            id_proc_E = floor(id_proc/(num_proc_per_var_f*num_proc_per_var_theta*num_proc_per_var_phi*num_proc_per_var_Triangle_gam1_gam2));
            id_proc_f = floor(id_proc/(num_proc_per_var_theta*num_proc_per_var_phi*num_proc_per_var_Triangle_gam1_gam2)) - id_proc_E*num_proc_per_var_f;
            id_proc_phi = floor(id_proc/(num_proc_per_var_theta*num_proc_per_var_Triangle_gam1_gam2)) - ((id_proc_E*num_proc_per_var_f+id_proc_f)*num_proc_per_var_phi);
            id_proc_theta = floor(id_proc/num_proc_per_var_Triangle_gam1_gam2) - ((id_proc_E*num_proc_per_var_f+id_proc_f)*num_proc_per_var_phi+id_proc_phi)*num_proc_per_var_theta;
            id_proc_Triangle_gam1_gam2 = id_proc - (((id_proc_E*num_proc_per_var_f + id_proc_f)*num_proc_per_var_phi + id_proc_phi)*num_proc_per_var_theta+id_proc_theta)*num_proc_per_var_Triangle_gam1_gam2;
        } else {
            id_proc_E = 0;
            id_proc_f = 0;
            id_proc_phi = 0;
            id_proc_theta = 0;
            id_proc_Triangle_gam1_gam2 = 0;
        }
    } else if ((Problem_Type == GRAY) // ||
               // (num_points->Maximum_Entropy_Solution_Regime != HYPERBOLIC_LIMIT) ||
               // (num_points->Maximum_Entropy_Solution_Regime != LOGARITHMIC_LIMIT)
              ) {
        // Compute number of processors per for each variable
        if (num_proc > 1) {
            num_proc_per_var_E = 1;
            num_proc_per_var_f = num_proc_f;
            num_proc_per_var_phi = num_proc_Phi;
            num_proc_per_var_theta = num_proc_Theta;
            num_proc_per_var_Triangle_gam1_gam2 = num_proc_Triangle_gam1_gam2;
        } else {
            num_proc_per_var_E = 1;
            num_proc_per_var_f = 1;
            num_proc_per_var_phi = 1;
            num_proc_per_var_theta = 1;
            num_proc_per_var_Triangle_gam1_gam2 = 1;
        }
        
        num_proc_used = num_proc_per_var_E*num_proc_per_var_f*num_proc_per_var_phi*num_proc_per_var_theta*num_proc_per_var_Triangle_gam1_gam2;
    
        // Determine id for each variable based on the current id_proc
        if (id_proc < num_proc_used) {
            id_proc_E = 0;
            id_proc_f = floor(id_proc/(num_proc_per_var_theta*num_proc_per_var_phi*num_proc_per_var_Triangle_gam1_gam2));   
            id_proc_phi = floor(id_proc/(num_proc_per_var_theta*num_proc_per_var_Triangle_gam1_gam2)) - id_proc_f*num_proc_per_var_phi;
            id_proc_theta = floor(id_proc/num_proc_per_var_Triangle_gam1_gam2) - (id_proc_f*num_proc_per_var_phi+id_proc_phi)*num_proc_per_var_theta;
            id_proc_Triangle_gam1_gam2 = id_proc - ((id_proc_f*num_proc_per_var_phi+id_proc_phi)*num_proc_per_var_theta+id_proc_theta)*num_proc_per_var_Triangle_gam1_gam2;
        } else {
            id_proc_E = 0;
            id_proc_f = 0;
            id_proc_phi = 0;
            id_proc_theta = 0;
            id_proc_Triangle_gam1_gam2 = 0;
        }
    } else {
        cout << "Invalid value for Problem_Type" << endl;
        exit(0);
    }
    
    size_rec = N_pts_total;
    size_rec /= num_proc_used;
    
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
    
    MPI_proc_parameters.size_rec = size_rec;
    
    if (id_proc == PRIMARY_ID) {
        cout << "num_proc_per_var_E = " << num_proc_per_var_E << "  " << "num_proc_per_var_f = " << num_proc_per_var_f << "  " << "num_proc_per_var_phi = " << num_proc_per_var_phi << "  " << "num_proc_per_var_theta = " << num_proc_per_var_theta << "  " << "num_proc_per_var_Triangle_gam1_gam2 = " << num_proc_per_var_Triangle_gam1_gam2 << endl;
    }
    
//     if (id_proc == 1) {
//         cout << "id_proc = " << id_proc << "  " << "id_proc_E = " << id_proc_E << "  " << "id_proc_f = " << id_proc_f << "  " << "id_proc_phi = " << id_proc_phi << "  " << "id_proc_theta = " << id_proc_theta << "  " << "id_proc_Triangle_gam1_gam2 = " << id_proc_Triangle_gam1_gam2 << endl;
//         cout << "num_proc_per_var_E = " << num_proc_per_var_E << "  " << "num_proc_per_var_f = " << num_proc_per_var_f << "  " << "num_proc_per_var_phi = " << num_proc_per_var_phi << "  " << "num_proc_per_var_theta = " << num_proc_per_var_theta << "  " << "num_proc_per_var_Triangle_gam1_gam2 = " << num_proc_per_var_Triangle_gam1_gam2 << endl;
//     }
    
    if (N_Points_E % num_proc_per_var_E != 0) {
        if (id_proc == PRIMARY_ID) {
            cout << "The number N_Points_E = " << N_Points_E << "is not divisible by num_proc_per_var_E = " << num_proc_per_var_E << "the quotient is " << N_Points_E/num_proc_per_var_E << endl;
        }
        exit(0);
    }
    
    if (N_Points_f % num_proc_per_var_f != 0) {
        if (id_proc == PRIMARY_ID) {
            cout << "The number N_Points_f = " << N_Points_f << "is not divisible by num_proc_per_var_f = " << num_proc_per_var_f << "the quotient is " << N_Points_f/num_proc_per_var_f << endl;
        }
        exit(0);
    }
    
    if (N_Points_Phi % num_proc_per_var_phi != 0) {
        if (id_proc == PRIMARY_ID) {
            cout << "The number N_Points_Phi = " << N_Points_Phi << "is not divisible by num_proc_per_var_phi = " << num_proc_per_var_phi << "the quotient is " << N_Points_Phi/num_proc_per_var_phi << endl;
        }
        exit(0);
    }
    
    if (N_Points_Theta % num_proc_per_var_theta != 0) {
        if (id_proc == PRIMARY_ID) {
            cout << "The number N_Points_Theta = " << N_Points_Theta << "is not divisible by num_proc_per_var_theta = " << num_proc_per_var_theta << "the quotient is " << N_Points_Theta/num_proc_per_var_theta << endl;
        }
        exit(0);
    }
    
    if (N_Points_Triangle_gam1_gam2 % num_proc_per_var_Triangle_gam1_gam2 != 0) {
        if (id_proc == PRIMARY_ID) {
            cout << "The number N_Points_Triangle_gam1_gam2 = " << N_Points_Triangle_gam1_gam2 << "is not divisible by num_proc_per_var_Triangle_gam1_gam2 = " << num_proc_per_var_Triangle_gam1_gam2 << "the quotient is " << N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2 << endl;
        }
        exit(0);
    }
    
    MPI_proc_parameters.N_Points_E_per_proc = N_Points_E/num_proc_per_var_E;
    MPI_proc_parameters.N_Points_f_per_proc = N_Points_f/num_proc_per_var_f;
    MPI_proc_parameters.N_Points_Phi_per_proc = N_Points_Phi/num_proc_per_var_phi;
    MPI_proc_parameters.N_Points_Theta_per_proc = N_Points_Theta/num_proc_per_var_theta;
    MPI_proc_parameters.N_Points_Triangle_gam1_gam2_Per_Proc = N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2;
    
    if (id_proc == PRIMARY_ID) {
        cout << endl;
        cout << "id_proc_E = " << id_proc_E << "  " << "id_proc_f = " << id_proc_f << "  " << "id_proc_phi = " << id_proc_phi << "  " << "id_proc_theta = " << id_proc_theta << "  " << "id_proc_Triangle_gam1_gam2 = " << id_proc_Triangle_gam1_gam2 << endl;
//         exit(0);
    }
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Compute_MPI_Processes_Max_Min_Indexes(int &id_E_min, int &id_E_max, int &id_f_min, int &id_f_max, int &id_phi_min, int &id_phi_max, int &id_theta_min, int &id_theta_max, int &id_Triangle_gam1_gam2_min, int &id_Triangle_gam1_gam2_max) {
    int num_proc_per_var_E, num_proc_per_var_f, num_proc_per_var_phi, num_proc_per_var_theta, num_proc_per_var_Triangle_gam1_gam2;
    int id_proc_E, id_proc_f, id_proc_theta, id_proc_phi, id_proc_Triangle_gam1_gam2;
    
    num_proc_per_var_E = MPI_proc_parameters.num_proc_per_var_E;
    num_proc_per_var_f = MPI_proc_parameters.num_proc_per_var_f;
    num_proc_per_var_phi = MPI_proc_parameters.num_proc_per_var_phi;
    num_proc_per_var_theta = MPI_proc_parameters.num_proc_per_var_theta;
    num_proc_per_var_Triangle_gam1_gam2 = MPI_proc_parameters.num_proc_per_var_Triangle_gam1_gam2;
    
    id_proc_E = MPI_proc_parameters.id_proc_E;
    id_proc_f = MPI_proc_parameters.id_proc_f;
    id_proc_theta = MPI_proc_parameters.id_proc_theta;
    id_proc_phi = MPI_proc_parameters.id_proc_phi;
    id_proc_Triangle_gam1_gam2 = MPI_proc_parameters.id_proc_Triangle_gam1_gam2;
    
    id_E_min = id_proc_E*N_Points_E/num_proc_per_var_E;
    id_E_max = (id_proc_E+1)*N_Points_E/num_proc_per_var_E;
    
    id_f_min = id_proc_f*N_Points_f/num_proc_per_var_f;
    id_f_max = (id_proc_f+1)*N_Points_f/num_proc_per_var_f;
    
    id_phi_min = id_proc_phi*N_Points_Phi/num_proc_per_var_phi;
    id_phi_max = (id_proc_phi+1)*N_Points_Phi/num_proc_per_var_phi;
    
    id_theta_min = id_proc_theta*N_Points_Theta/num_proc_per_var_theta;
    id_theta_max = (id_proc_theta+1)*N_Points_Theta/num_proc_per_var_theta;
    
    id_Triangle_gam1_gam2_min = id_proc_Triangle_gam1_gam2*N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2;
    id_Triangle_gam1_gam2_max = (id_proc_Triangle_gam1_gam2+1)*N_Points_Triangle_gam1_gam2/num_proc_per_var_Triangle_gam1_gam2;
    
    // cout << "id_proc = " << id_proc << "  " << "id_proc_f = " << id_proc_f << "  " << "N_Points_f = " << N_Points_f << "  " << "num_proc_per_var_f = " << num_proc_per_var_f << endl;
    
    // cout << "id_E_min = " << id_E_min << "  " << "id_E_max = " << id_E_max << "  " << "id_f_min = " << id_f_min << "  " << "id_f_max = " << id_f_max << "  " << "id_phi_min = " << id_phi_min << "  " << "id_phi_max = " << id_phi_max << "  " << "id_theta_min = " << id_theta_min << "  " << "id_theta_max = " << id_theta_max << "  " << "id_Triangle_gam1_gam2_min = " << id_Triangle_gam1_gam2_min << "  " << "id_Triangle_gam1_gam2_max = " << id_Triangle_gam1_gam2_max << endl;
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Create_MPI_Data_Type_rec_N3() {
    // Create the rec_N3 datatype for MPI
    int lengths[1] = { 98 };
    const MPI_Aint displacements[1] = { 0 };
    MPI_Datatype types[1] = { MPI_LONG_DOUBLE };
    MPI_Type_create_struct(1, lengths, displacements, types, &rec_N3_type);
    MPI_Type_commit(&rec_N3_type);
}

//******************************************************************************************
// This routines sets up the path to the file in which precomputed maximum entropy solutions
// are stored and then open the latter
//******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: OpenInputFile(char *filename) {
    strcpy(path_out, getenv(PATHVAR));
    strcat(path_out, "/M2_Model/");
    strcat(path_out, filename);
    sprintf(extension, "_%.6d", 0);
    strcat(extension, ".dat");
    strcat(path_out, extension);
        
    in_out.open(path_out, ios::in|ios::binary);
        
    if (!in_out) {
        cout << filename << ".dat could not be accessed!" << endl;  
        exit(0); 
    }
}

//*****************************************************************************
// This routines closes the file in which precomputed maximum entropy solutions
//*****************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: CloseInputFile() {
    if (in_out.good()) {
        in_out.close();
    }
}

//******************************************************************************************
// This routines reads the file containing precomputed maximum entropy solutions and stores
// them in the data structure N3_Non_Gray_M2_3D_RT_Cheby
//******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: ReadInputData() {
    int index;
    
    if (id_proc == PRIMARY_ID) {
        if (in_out.good()) {
            for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
                // cout << "id_Mobius = " << id_Mobius << endl;
                for (int index_e = 0 ; index_e < N_Points_E; index_e++){
                    for (int i_f = 0; i_f < N_Points_f; i_f++) {
                        for (int i_Phi = 0 ; i_Phi < N_Points_Phi; i_Phi++) {
                            for (int i_Theta = 0 ; i_Theta < N_Points_Theta; i_Theta++) {
                                for (int i_gam1_gam2 = 0 ; i_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_gam1_gam2++) {
                                    index = Compute_Full_Id(id_Mobius, index_e, i_f, i_Phi, i_Theta, i_gam1_gam2);
                                    
                                    read<record_N3, record_Ncoeffs>(in_out, rec_N3[index], index_rec_N3, index_Ncoeffs);
                                    
                                    // cout << "id_Mobius = " << id_Mobius << "   " << "ratio_I0 = " << rec_N3[index].ratio_I0 << "   " << "I0 = " << rec_N3[index].I0 << "   " << "N1_1 = " << rec_N3[index].N1_1 << "   " << "N1_2 = " << rec_N3[index].N1_2 << "   " << "N1_3 = " << rec_N3[index].N1_3 << "   " << "gam1 = " << rec_N3[index].gam1 << "   " << "gam2 = " << rec_N3[index].gam2 << endl;
                                    
                                    index_rec_N3 = index_rec_N3 + 1;
                                } // end for i_gam1_gam2
                            } // end for i_Theta
                        } // end for i_Phi
                    } // end for i_f
                } // end for index_e
            } // end for id_Mobius
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    int Npts_Total = Compute_Num_Pts_Full();
    // Broadcast rec_N3 to other processors
    MPI_Bcast(rec_N3, Npts_Total, rec_N3_type, PRIMARY_ID, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    Read_Maximum_Entropy_Data(rec_N3);
    
    MPI_Barrier(MPI_COMM_WORLD);
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Write_Maximum_Entropy_Data(const record_N3 *rec_N3_global) {
    int index;
    
    strcpy(path_out, getenv(PATHVAR));
    strcat(path_out, "/M2_Model/");
    if (Problem_Type == GRAY) {
        strcat(path_out, "Gray_M2_Model/Cheby_N3_Gray_M2_3D");
    } else if (Problem_Type == NON_GRAY) {
        strcat(path_out, "Non_Gray_M2_Model/Cheby_N3_Non_Gray_M2_3D");
    } else {
        cout << "Problem type not specified !!!!!!!!!!!!!!!!!" << endl;
        exit(0);
    }
    sprintf(extension, "_%.6d", 0);
    strcat(extension, ".dat");
    strcat(path_out, extension);
        
    in_out.open(path_out, ios::out|ios::binary);
        
    if (!in_out) {
        cout << "Cheby_N3_Gray_M2_3D.dat could not be accessed!" << endl;  
        exit(0); 
    }
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        for (int id_E = 0; id_E < N_Points_E; id_E++) {
            for (int id_f = 0; id_f < N_Points_f; id_f++) {
                for (int i_Phi = 0; i_Phi < N_Points_Phi; i_Phi++) {
                    for (int i_Theta = 0; i_Theta < N_Points_Theta; i_Theta++) {
                        for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_gam1_gam2++) {
                            index = Compute_Full_Id(id_Mobius, id_E, id_f, i_Phi, i_Theta, i_gam1_gam2);
                            
                            // cout << "id_Mobius = " << id_Mobius << "   " << "ratio_I0 = " << rec_N3_global[index].ratio_I0 << "   " << "N1_1 = " << rec_N3_global[index].N1_1 << "   " << "N1_2 = " << rec_N3_global[index].N1_2 << "   " << "N1_3 = " << rec_N3_global[index].N1_3 << "   " << "gam1 = " << rec_N3_global[index].gam1 << "   " << "gam2 = " << rec_N3_global[index].gam2 << endl;
                            
                            write<record_N3>(in_out, rec_N3_global[index]);
                        } // end for i_gam1_gam2
                    } // end for i_Theta
                } // end for i_Phi
            } // end for id_f
        } // end for id_E
    } // end for id_Mobius
    
    in_out.close();
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Read_Maximum_Entropy_Data(const record_N3 *rec_N3_global) {
    int index;
    
    // double N3_111_LB, N3_111_UB;
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        for (int id_E = 0; id_E < N_Points_E; id_E++) {
            for (int id_f = 0; id_f < N_Points_f; id_f++) {
                for (int i_Phi = 0; i_Phi < N_Points_Phi; i_Phi++) {
                    for (int i_Theta = 0; i_Theta < N_Points_Theta; i_Theta++) {
                        for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_gam1_gam2++) {
                            index = Compute_Full_Id(id_Mobius, id_E, id_f, i_Phi, i_Theta, i_gam1_gam2);
                            
                            E_NON_GRAY[index] = rec_N3_global[index].ratio_I0;
                            N1_1_NON_GRAY[index] = rec_N3_global[index].N1_1;
                            N1_2_NON_GRAY[index] = rec_N3_global[index].N1_2;
                            N1_3_NON_GRAY[index] = rec_N3_global[index].N1_3;
                            gam1_NON_GRAY[index] = rec_N3_global[index].gam1;
                            gam2_NON_GRAY[index] = rec_N3_global[index].gam2;
                            
                            x_SH[index] = rec_N3_global[index].x_SH;
                            y_SH[index] = rec_N3_global[index].y_SH;
                            z_SH[index] = rec_N3_global[index].z_SH;
                            
                            N3_111_NON_GRAY[index] = rec_N3_global[index].N3_111;
                            N3_112_NON_GRAY[index] = rec_N3_global[index].N3_112;
                            N3_122_NON_GRAY[index] = rec_N3_global[index].N3_122;
                            N3_123_NON_GRAY[index] = rec_N3_global[index].N3_123;
                            N3_222_NON_GRAY[index] = rec_N3_global[index].N3_222;
                            
//                             N3_111_LB = pow(N1_1_NON_GRAY[index], 3);
//                             N3_111_UB = N1_1_NON_GRAY[index]*(pow(N1_1_NON_GRAY[index], 2) + (1.0 - pow(N1_1_NON_GRAY[index], 2) - pow(N1_2_NON_GRAY[index], 2) - pow(N1_3_NON_GRAY[index], 2)));
//                             if (N3_111_NON_GRAY[index] < N3_111_LB || N3_111_NON_GRAY[index] > N3_111_UB) {
//                                 cout << "N3_111_LB = " << N3_111_LB << "  " << "N3_111_UB = " << N3_111_UB << "  " << "N3_111_NON_GRAY = " << N3_111_NON_GRAY[index] << endl;
//                             }
                            
                            dN3_111_NON_GRAY_dN1_1[index] = rec_N3_global[index].dN3_111_dN1_1;
                            dN3_111_NON_GRAY_dN1_2[index] = rec_N3_global[index].dN3_111_dN1_2;
                            dN3_111_NON_GRAY_dN1_3[index] = rec_N3_global[index].dN3_111_dN1_3;
                            dN3_111_NON_GRAY_dgam1[index] = rec_N3_global[index].dN3_111_dgam1;
                            dN3_111_NON_GRAY_dnorm_f[index] = rec_N3_global[index].dN3_111_dnorm_f;
                            d2_N3_111_NON_GRAY_dnorm_f_dN1_1[index] = rec_N3_global[index].d2_N3_111_dnorm_f_dN1_1;
                            d2_N3_111_NON_GRAY_dnorm_f_dgam1[index] = rec_N3_global[index].d2_N3_111_dnorm_f_dgam1;
                            d2_N3_111_NON_GRAY_dgam1_dN1_1[index] = rec_N3_global[index].d2_N3_111_dgam1_dN1_1;
                            d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index] = rec_N3_global[index].d3_N3_111_dnorm_f_dgam1_dN1_1;
                            
                            dN3_122_NON_GRAY_dN1_1[index] = rec_N3_global[index].dN3_122_dN1_1;
                            dN3_122_NON_GRAY_dN1_2[index] = rec_N3_global[index].dN3_122_dN1_2;
                            dN3_122_NON_GRAY_dN1_3[index] = rec_N3_global[index].dN3_122_dN1_3;
                            dN3_122_NON_GRAY_dgam1[index] = rec_N3_global[index].dN3_122_dgam1;
                            dN3_122_NON_GRAY_dgam2[index] = rec_N3_global[index].dN3_122_dgam2;
                            dN3_122_NON_GRAY_dnorm_f[index] = rec_N3_global[index].dN3_122_dnorm_f;
                            d2_N3_122_NON_GRAY_dgam1_dgam2[index] = rec_N3_global[index].d2_N3_122_dgam1_dgam2;
                            d2_N3_122_NON_GRAY_dnorm_f_dN1_1[index] = rec_N3_global[index].d2_N3_122_dnorm_f_dN1_1;
                            d2_N3_122_NON_GRAY_dnorm_f_dgam1[index] = rec_N3_global[index].d2_N3_122_dnorm_f_dgam1;
                            d2_N3_122_NON_GRAY_dnorm_f_dgam2[index] = rec_N3_global[index].d2_N3_122_dnorm_f_dgam2;
                            d2_N3_122_NON_GRAY_dgam1_dN1_1[index] = rec_N3_global[index].d2_N3_122_dgam1_dN1_1;
                            d2_N3_122_NON_GRAY_dgam2_dN1_1[index] = rec_N3_global[index].d2_N3_122_dgam2_dN1_1;
                            d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1[index] = rec_N3_global[index].d3_N3_122_dnorm_f_dgam1_dN1_1;
                            d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1[index] = rec_N3_global[index].d3_N3_122_dnorm_f_dgam2_dN1_1;
                            d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2[index] = rec_N3_global[index].d3_N3_122_dnorm_f_dgam1_dgam2;
                            d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1[index] = rec_N3_global[index].d3_N3_122_dgam1_dgam2_dN1_1;
                            d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1[index] = rec_N3_global[index].d4_N3_122_dnorm_f_dgam1_dgam2_dN1_1;
                            
                            dN3_123_NON_GRAY_dN1_1[index] = rec_N3_global[index].dN3_123_dN1_1;
                            dN3_123_NON_GRAY_dN1_2[index] = rec_N3_global[index].dN3_123_dN1_2;
                            dN3_123_NON_GRAY_dN1_3[index] = rec_N3_global[index].dN3_123_dN1_3;
                            dN3_123_NON_GRAY_dgam1[index] = rec_N3_global[index].dN3_123_dgam1;
                            dN3_123_NON_GRAY_dgam2[index] = rec_N3_global[index].dN3_123_dgam2;
                            dN3_123_NON_GRAY_dgam3[index] = rec_N3_global[index].dN3_123_dgam3;
                            d2_N3_123_NON_GRAY_dN1_1_dN1_2[index] = rec_N3_global[index].d2_N3_123_dN1_1_dN1_2;
                            d2_N3_123_NON_GRAY_dN1_1_dN1_3[index] = rec_N3_global[index].d2_N3_123_dN1_1_dN1_3;
                            d2_N3_123_NON_GRAY_dN1_2_dN1_3[index] = rec_N3_global[index].d2_N3_123_dN1_2_dN1_3;
                            d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3[index] = rec_N3_global[index].d3_N3_123_dN1_1_dN1_2_dN1_3;
                            d2_N3_123_NON_GRAY_dgam1_dN1_1[index] = rec_N3_global[index].d2_N3_123_dgam1_dN1_1;
                            d2_N3_123_NON_GRAY_dgam2_dN1_1[index] = rec_N3_global[index].d2_N3_123_dgam2_dN1_1;
                            d2_N3_123_NON_GRAY_dN1_1_dgam3[index] = rec_N3_global[index].d2_N3_123_dN1_1_dgam3;
                            d2_N3_123_NON_GRAY_dN1_2_dgam1[index] = rec_N3_global[index].d2_N3_123_dN1_2_dgam1;
                            d2_N3_123_NON_GRAY_dN1_2_dgam2[index] = rec_N3_global[index].d2_N3_123_dN1_2_dgam2;
                            d2_N3_123_NON_GRAY_dN1_2_dgam3[index] = rec_N3_global[index].d2_N3_123_dN1_2_dgam3;
                            d2_N3_123_NON_GRAY_dN1_3_dgam1[index] = rec_N3_global[index].d2_N3_123_dN1_3_dgam1;
                            d2_N3_123_NON_GRAY_dN1_3_dgam2[index] = rec_N3_global[index].d2_N3_123_dN1_3_dgam2;
                            d2_N3_123_NON_GRAY_dN1_3_dgam3[index] = rec_N3_global[index].d2_N3_123_dN1_3_dgam3;
                            d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1[index] = rec_N3_global[index].d3_N3_123_dN1_1_dN1_2_dgam1;
                            d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1[index] = rec_N3_global[index].d3_N3_123_dN1_1_dN1_3_dgam1;
                            d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1[index] = rec_N3_global[index].d3_N3_123_dN1_2_dN1_3_dgam1;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_2_dN1_3_dgam1;
                            d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2[index] = rec_N3_global[index].d3_N3_123_dN1_1_dN1_2_dgam2;
                            d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2[index] = rec_N3_global[index].d3_N3_123_dN1_1_dN1_3_dgam2;
                            d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2[index] = rec_N3_global[index].d3_N3_123_dN1_2_dN1_3_dgam2;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_2_dN1_3_dgam2;
                            d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3[index] = rec_N3_global[index].d3_N3_123_dN1_1_dN1_2_dgam3;
                            d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3[index] = rec_N3_global[index].d3_N3_123_dN1_1_dN1_3_dgam3;
                            d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3[index] = rec_N3_global[index].d3_N3_123_dN1_2_dN1_3_dgam3;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_2_dN1_3_dgam3;
                            d2_N3_123_NON_GRAY_dgam1_dgam2[index] = rec_N3_global[index].d2_N3_123_dgam1_dgam2;
                            d2_N3_123_NON_GRAY_dgam1_dgam3[index] = rec_N3_global[index].d2_N3_123_dgam1_dgam3;
                            d2_N3_123_NON_GRAY_dgam2_dgam3[index] = rec_N3_global[index].d2_N3_123_dgam2_dgam3;
                            d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1[index] = rec_N3_global[index].d3_N3_123_dgam1_dgam2_dN1_1;
                            d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3[index] = rec_N3_global[index].d3_N3_123_dgam1_dN1_1_dgam3;
                            d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3[index] = rec_N3_global[index].d3_N3_123_dgam2_dN1_1_dgam3;
                            d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2[index] = rec_N3_global[index].d3_N3_123_dN1_2_dgam1_dgam2;
                            d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3[index] = rec_N3_global[index].d3_N3_123_dN1_2_dgam1_dgam2;
                            d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3[index] = rec_N3_global[index].d3_N3_123_dN1_2_dgam2_dgam3;
                            d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2[index] = rec_N3_global[index].d3_N3_123_dN1_3_dgam1_dgam2;
                            d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3[index] = rec_N3_global[index].d3_N3_123_dN1_3_dgam1_dgam3;
                            d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3[index] = rec_N3_global[index].d3_N3_123_dN1_3_dgam2_dgam3;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_2_dgam1_dgam2;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_2_dgam1_dgam3;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_2_dgam2_dgam3;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_3_dgam1_dgam2;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_3_dgam1_dgam3;
                            d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3[index] = rec_N3_global[index].d4_N3_123_dN1_1_dN1_3_dgam2_dgam3;
                            d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2[index] = rec_N3_global[index].d4_N3_123_dN1_2_dN1_3_dgam1_dgam2;
                            d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3[index] = rec_N3_global[index].d4_N3_123_dN1_2_dN1_3_dgam1_dgam3;
                            d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3[index] = rec_N3_global[index].d4_N3_123_dN1_2_dN1_3_dgam2_dgam3;
                            d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2[index] = rec_N3_global[index].d5_N3_123_dN1_1_dN1_2_dN1_3_dgam1_dgam2;
                            d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3[index] = rec_N3_global[index].d5_N3_123_dN1_1_dN1_2_dN1_3_dgam1_dgam3;
                            d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3[index] = rec_N3_global[index].d5_N3_123_dN1_1_dN1_2_dN1_3_dgam2_dgam3;
                            
                            dN3_111_NON_GRAY_dmu[index] = Cartesian_to_spherical_Coordinates_Jacobian(rec_N3_global[index].N1,
                                                                                                      rec_N3_global[index].x_SH,
                                                                                                      rec_N3_global[index].y_SH,
                                                                                                      rec_N3_global[index].z_SH,
                                                                                                      rec_N3_global[index].dN3_111_dN1_1, 
                                                                                                      rec_N3_global[index].dN3_111_dN1_2,
                                                                                                      rec_N3_global[index].dN3_111_dN1_3,
                                                                                                      1);
                            
                            // cout << "N1 = " << rec_N3_global[index].N1 << " " << "x_SH = " << rec_N3_global[index].x_SH << " " << "y_SH = " << rec_N3_global[index].y_SH << " " << "z_SH = " << rec_N3_global[index].z_SH << " " << "dN3_111_dN1_1 = " << rec_N3_global[index].dN3_111_dN1_1 << " " << "dN3_111_dN1_2 = " << rec_N3_global[index].dN3_111_dN1_2 << " " << "dN3_111_dN1_3 = " << rec_N3_global[index].dN3_111_dN1_3 << "  " << "dN3_111_NON_GRAY_dmu = " << dN3_111_NON_GRAY_dmu[index] << endl;
                            
                            dN3_122_NON_GRAY_dmu[index] = Cartesian_to_spherical_Coordinates_Jacobian(rec_N3_global[index].N1,
                                                                                                      rec_N3_global[index].x_SH,
                                                                                                      rec_N3_global[index].y_SH,
                                                                                                      rec_N3_global[index].z_SH,
                                                                                                      rec_N3_global[index].dN3_122_dN1_1, 
                                                                                                      rec_N3_global[index].dN3_122_dN1_2,
                                                                                                      rec_N3_global[index].dN3_122_dN1_3,
                                                                                                      1);
                            
                            dN3_123_NON_GRAY_dmu[index] = Cartesian_to_spherical_Coordinates_Jacobian(rec_N3_global[index].N1,
                                                                                                      rec_N3_global[index].x_SH,
                                                                                                      rec_N3_global[index].y_SH,
                                                                                                      rec_N3_global[index].z_SH,
                                                                                                      rec_N3_global[index].dN3_123_dN1_1, 
                                                                                                      rec_N3_global[index].dN3_123_dN1_2,
                                                                                                      rec_N3_global[index].dN3_123_dN1_3,
                                                                                                      1);
                            
//                             if (id_proc == PRIMARY_ID) {
//                                 cout << "index = " << index << "   " << "E = " << E_NON_GRAY[index] << "   " << "x_SH = " << x_SH[index] << "   " << "y_SH = " << y_SH[index] << "   " << "z_SH = " << z_SH[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "   " << "N3_111 = " << N3_111_NON_GRAY[index] << "  " << "N3_122 = " << N3_122_NON_GRAY[index] << endl;
//                             } // end if    
                        } // end for i_gam1_gam2
                    } // end for i_Theta
                } // end for i_Phi
            } // end for id_f
        } // end for id_E
    } // end for id_Mobius
}

int N3_Non_Gray_M2_3D_RT_Cheby :: Compute_Num_Pts_Full() {
    int Npts_total;
    
    Npts_total = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Theta*N_Points_Phi*N_Points_Triangle_gam1_gam2;
    
    return Npts_total;
}

int N3_Non_Gray_M2_3D_RT_Cheby :: Compute_Num_Coeffs_Full_Without_SH() {
    int N_Coeffs;
    
    N_Coeffs = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
    
    return N_Coeffs;
}

int N3_Non_Gray_M2_3D_RT_Cheby :: Compute_Num_Coeffs_Full() {
    int N_Coeffs;
    
    N_Coeffs = N_Points_E*N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
    
    return N_Coeffs;
}

int N3_Non_Gray_M2_3D_RT_Cheby :: Compute_Full_Id_Triangle(const int &i_gam1, const int &i_gam2) {
    int index;
    index = ((i_gam1+1)*(i_gam1+2))/2 + i_gam2;
    index = N_Points_Triangle_gam1_gam2 - index;
    return index;
}

int N3_Non_Gray_M2_3D_RT_Cheby :: Compute_Id_Triangle_Single_Block(const int &i_gam1, const int &i_gam2) {
    int index;
    
//     index = ((i_gam1+1)*(i_gam1+2))/2 + i_gam2;
//     index = N_Points_Triangle_gam1_gam2 - index;
//     return index;
    cout << "Incorrect !!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    exit(0);
    
    return index;
}

int N3_Non_Gray_M2_3D_RT_Cheby :: Compute_N_pts_Triangle_Single_Block(const int &id_Triangle_gam1_gam2_min, const int &id_Triangle_gam1_gam2_max) {
    int N_Points_Triangle_Per_Proc;
    N_Points_Triangle_Per_Proc = 0;
    for (int i_gam1 = id_Triangle_gam1_gam2_min; i_gam1 < id_Triangle_gam1_gam2_max; i_gam1++) {
        for (int i_gam2 = 0; i_gam2 < N_Points_gam1 - i_gam1; i_gam2++) {
            N_Points_Triangle_Per_Proc++;
        }
    }
    return N_Points_Triangle_Per_Proc;
}

int N3_Non_Gray_M2_3D_RT_Cheby :: Compute_Full_Id(const int &id_Mobius, const int &i_Cheby_E, const int &i_Cheby_f, const int &i_Phi, const int &i_Theta, const int &i_Cheby_gam1_gam2) {
    int index;
    index = id_Mobius*N_Points_E + i_Cheby_E;
    index = index*N_Points_f + i_Cheby_f;
    index = (index*N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
    index = index*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2;
    return index;
}

int N3_Non_Gray_M2_3D_RT_Cheby :: Compute_Full_Id_Single_proc(const int &id_Mobius, const int &i_Cheby_E, const int &i_Cheby_f, const int &i_Phi, const int &i_Theta, const int &i_Cheby_gam1_gam2) {
    int index;
    int id_E_min, id_E_max;
    int id_f_min, id_f_max;
    int id_phi_min, id_phi_max;
    int id_theta_min, id_theta_max;
    int id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max;
    int i_Cheby_E_temp, i_Cheby_f_temp, i_Phi_temp, i_Theta_Temp, i_Triangle_gam1_gam2_Temp;
    int N_Points_E_per_proc, N_Points_f_per_proc, N_Points_Phi_per_proc, N_Points_Theta_per_proc, N_Points_Triangle_gam1_gam2_Per_Proc;
    
    N_Points_E_per_proc = MPI_proc_parameters.N_Points_E_per_proc;
    N_Points_f_per_proc = MPI_proc_parameters.N_Points_f_per_proc;
    N_Points_Phi_per_proc = MPI_proc_parameters.N_Points_Phi_per_proc;
    N_Points_Theta_per_proc = MPI_proc_parameters.N_Points_Theta_per_proc;
    N_Points_Triangle_gam1_gam2_Per_Proc = MPI_proc_parameters.N_Points_Triangle_gam1_gam2_Per_Proc;
    
    Compute_MPI_Processes_Max_Min_Indexes(id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, id_theta_min, id_theta_max, id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max);
    
    // Compute ralative Id first
    i_Cheby_E_temp = i_Cheby_E - id_E_min;
    i_Cheby_f_temp = i_Cheby_f - id_f_min;
    i_Phi_temp = i_Phi - id_phi_min;
    i_Theta_Temp = i_Theta - id_theta_min;
    i_Triangle_gam1_gam2_Temp = i_Cheby_gam1_gam2 - id_Triangle_gam1_gam2_min;
    
    index = id_Mobius*N_Points_E_per_proc + i_Cheby_E_temp;
    index = index*N_Points_f_per_proc + i_Cheby_f_temp;
    index = (index*N_Points_Phi_per_proc + i_Phi_temp)*N_Points_Theta_per_proc + i_Theta_Temp;
    index = index*N_Points_Triangle_gam1_gam2_Per_Proc + i_Triangle_gam1_gam2_Temp;
    
//     if (id_proc == 0) {
//         cout << "id_Mobius = " << id_Mobius << endl;
//         
//         cout << "N_Points_E_per_proc = " << N_Points_E_per_proc << "  " << "N_Points_f_per_proc = " << N_Points_f_per_proc << "  " << "N_Points_Phi_per_proc = " << N_Points_Phi_per_proc << "  " << "N_Points_Theta_per_proc = " << N_Points_Theta_per_proc << endl;
//         
//         cout << "i_Cheby_E_temp = " << i_Cheby_E_temp << "  " << "i_Cheby_f_temp = " << i_Cheby_f_temp << "  " << "i_Phi_temp = " << i_Phi_temp << "  " << "i_Theta_Temp = " << i_Theta_Temp << endl;
//     }
    
    return index;
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Print_Data(const int &index) {
    cout << "index = " << index << "   " << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "  " << "N3_111 = " << N3_111_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY[index] << "  " << "N3_122 = " << N3_122_NON_GRAY[index] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index] << endl;
}

//********************************************************************************************************
// This routine sets up values of the weighting function, g_N3_ijk, of the polynomial interpolant, 
// f_N3_ijk, for the purpose of our interpolative-based approximations of the closing fluxes for either
// the gray or the non-gray M2 closure
//********************************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: SetupInterpolant_Values_BE() {
    int index;
    long double mu;
    long double norm_f, norm_f_2, gam3;
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
            for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
                for (int i_Phi = 0; i_Phi < N_Points_Phi; i_Phi++) {
                    for (int i_Theta = 0; i_Theta < N_Points_Theta; i_Theta++) {
                        for (int i_Cheby_gam1_gam2 = 0; i_Cheby_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_Cheby_gam1_gam2++) {
                            
                            index = Compute_Full_Id(id_Mobius, i_Cheby_E, i_Cheby_f, i_Phi, i_Theta, i_Cheby_gam1_gam2);
                            
                            norm_f_2 = pow(N1_1_NON_GRAY[index], 2) + pow(N1_2_NON_GRAY[index], 2) + pow(N1_3_NON_GRAY[index], 2);
                            
                            norm_f = sqrt(norm_f_2);
                            
                            mu = x_SH[index];
                            
                            // Compute f_N3_111
                            if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(fabs(1.0 - gam1_NON_GRAY[index]) < 1.0e-8)) {
                                // gam1 > 0 and gam1 < 1
                                if (norm_f < 1.0e-8) {
                                    // norm_f = 0
                                    f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dN1_1[index];
                                } else if (fabs(1.0 - norm_f) < 1.0e-8) {
                                    // norm_f = 1
                                    if (fabs(mu) < 1.0e-8) {
                                        // mu = 0
                                        cout << "Fix this 00000000000000 !!!!!!!!!!!!!!!!" << endl;
                                        Print_Data(index);
                                        // f_N3_111_NON_GRAY[index] = -(1.0/2.0)*d2_N3_111_NON_GRAY_dnorm_f_dmu[index];
                                    } else {
                                        // mu ~= 0
                                        f_N3_111_NON_GRAY[index] = (3.0/2.0)*pow(mu, 2) - (1.0/2.0)*dN3_111_NON_GRAY_dnorm_f[index]/mu;   
                                    }
                                } else {
                                    // norm_f > 0 && norm_f < 1
                                    if (fabs(mu) < 1.0e-8) {
                                        // mu = 0
                                        f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dmu[index];
                                        f_N3_111_NON_GRAY[index] /= norm_f*(1.0 - norm_f_2);
                                        
                                        // cout << "mu = " << mu << "   " << "N1_1_NON_GRAY[index] = " << N1_1_NON_GRAY[index] << "   " << "norm_f = " << norm_f << "   " << "dN3_111_NON_GRAY_dmu = " << dN3_111_NON_GRAY_dmu[index] << endl;
                                    } else {
                                        // mu ~= 0
                                        f_N3_111_NON_GRAY[index] = (N3_111_NON_GRAY[index] - pow(N1_1_NON_GRAY[index], 3));
                                        f_N3_111_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);  
                                    }
                                }
                                f_N3_111_NON_GRAY[index] = f_N3_111_NON_GRAY[index] - gam1_NON_GRAY[index];
                                f_N3_111_NON_GRAY[index] /= gam1_NON_GRAY[index]*(1.0 - gam1_NON_GRAY[index]);
                            } else if (gam1_NON_GRAY[index] < 1.0e-8) {
                                cout << "Fix this first 11111111111111111" << endl;
                                ///////////////////////////
                                Print_Data(index);
                                exit(0);
                                ///////////////////////////////////
                                // gam1 = 0
                                if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
                                    // N1_1 = 0
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_111_NON_GRAY[index] = -d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;
                                    } else {
                                        // norm_f < 1
                                        f_N3_111_NON_GRAY[index] = d2_N3_111_NON_GRAY_dgam1_dN1_1[index]/(1.0 - norm_f_2);
                                    }
                                } else {
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_111_NON_GRAY[index] = -d2_N3_111_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
                                    } else {
                                        // norm_f < 1
                                        f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dgam1[index];
                                        f_N3_111_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
                                    }
                                }
                                f_N3_111_NON_GRAY[index] = (f_N3_111_NON_GRAY[index] - 1.0)/(1.0 - gam1_NON_GRAY[index]);
                            } else if (fabs(1.0 - gam1_NON_GRAY[index]) < 1.0e-8) {
                                cout << "Fix this first 222222222222222222222" << endl;
                                ///////////////////////////
                                Print_Data(index);
                                exit(0);
                                ///////////////////////////////////
                                // gam1 = 1
                                if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
                                    // N1_1 = 0
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_111_NON_GRAY[index] = -d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;  
                                    } else {
                                        // norm_f < 1
                                        f_N3_111_NON_GRAY[index] = d2_N3_111_NON_GRAY_dgam1_dN1_1[index]/(1.0 - norm_f_2);
                                    }
                                } else {
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_111_NON_GRAY[index] = -d2_N3_111_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
                                    } else {
                                        // norm_f < 1
                                        f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dgam1[index];
                                        f_N3_111_NON_GRAY[index] = N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
                                    }
                                }
                                f_N3_111_NON_GRAY[index] = (1.0 - f_N3_111_NON_GRAY[index])/gam1_NON_GRAY[index];
                            } else {
                                cout << "Error in Setup Interpolant Values for N3 111" << endl;
                                exit(0);
                            }
                            
                            // Compute f_N3_122
                            if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8)) {
                                // gam1 > 0 and gam2 > 0
                                if (norm_f < 1.0e-8) {
                                    // norm_f = 0
                                    f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dN1_1[index];
                                } else if (fabs(1.0 - norm_f) < 1.0e-8) {
                                    // norm_f = 1
                                    if (fabs(mu) < 1.0e-8) {
                                        // mu = 0
                                        cout << "Fix this !!!!!!!!!!!!!!!!" << endl;
                                        Print_Data(index);
                                        // f_N3_122_NON_GRAY[index] = (3.0/2.0)*pow(N1_2_NON_GRAY[index], 2) - (1.0/2.0)*d2_N3_122_NON_GRAY_dnorm_f_dmu[index];
                                    } else {
                                        // mu ~= 0
                                        f_N3_122_NON_GRAY[index] = (3.0/2.0)*pow(N1_2_NON_GRAY[index], 2) - (1.0/2.0)*dN3_122_NON_GRAY_dnorm_f[index]/mu;   
                                    }
                                } else {
                                    // norm_f > 0 && norm_f < 1
                                    if (fabs(mu) < 1.0e-8) {
                                        // mu = 0
                                        f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dmu[index]/norm_f - pow(N1_2_NON_GRAY[index], 2);
                                        f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
                                    } else {
                                        // mu ~= 0
                                        f_N3_122_NON_GRAY[index] = (N3_122_NON_GRAY[index] - N1_1_NON_GRAY[index]*pow(N1_2_NON_GRAY[index], 2));
                                        f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
                                    }
                                }
                                f_N3_122_NON_GRAY[index] = f_N3_122_NON_GRAY[index] - gam2_NON_GRAY[index];
                                f_N3_122_NON_GRAY[index] /= gam2_NON_GRAY[index]*(1.0 - gam2_NON_GRAY[index]);
                                // f_N3_122_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam2_NON_GRAY[index];
                            } else if ((gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8)) {
                                cout << "Fix this first 3333333333333333" << endl;
                                Print_Data(index);
                                // gam1 = 0 and gam2 > 0 and gam2 < 1 - gam1
                                if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
                                    // N1_1 = 0
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;
                                    } else {
                                        // norm_f < 1
                                        f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam1_dN1_1[index];
                                        f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
                                    }
                                } else {
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_122_NON_GRAY[index] = -d2_N3_122_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
                                    } else {
                                        // norm_f < 1
                                        f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dgam1[index];
                                        f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
                                    }
                                }
                                f_N3_122_NON_GRAY[index] /= gam2_NON_GRAY[index];
                            } else if (!(gam1_NON_GRAY[index] < 1.0e-8) && (gam2_NON_GRAY[index] < 1.0e-8)) {
                                cout << "Fix this first 4444444444444444444" << endl;
                                Print_Data(index);
                                // gam1 > 0 and gam1 < 1 and gam2 = 0
                                if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
                                    // N1_1 = 0
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1[index]/2.0;
                                    } else {
                                        // norm_f < 1
                                        f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam2_dN1_1[index];
                                        f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
                                    }
                                } else {
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_122_NON_GRAY[index] = -d2_N3_122_NON_GRAY_dnorm_f_dgam2[index]/(2.0*N1_1_NON_GRAY[index]);
                                    } else {
                                        // norm_f < 1
                                        f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dgam2[index];
                                        f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
                                    }
                                }
                                f_N3_122_NON_GRAY[index] = (f_N3_122_NON_GRAY[index] - 1.0)/gam1_NON_GRAY[index];
                            } else if ((gam1_NON_GRAY[index] < 1.0e-8) && (gam2_NON_GRAY[index] < 1.0e-8)) {
                                cout << "Fix this first 555555555555555555" << endl;
                                Print_Data(index);
                                // gam1 = 0 && gam2 = 0
                                if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
                                    // N1_1 = 0
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_122_NON_GRAY[index] = -d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1[index]/2.0;
                                    } else {
                                        // norm_f < 1
                                        f_N3_122_NON_GRAY[index] = d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1[index];
                                        f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
                                    }
                                } else {
                                    if (fabs(1.0 - norm_f) < 1.0e-8) {
                                        // norm_f = 1
                                        f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2[index]/(2.0*N1_1_NON_GRAY[index]);
                                    } else {
                                        // norm_f < 1
                                        f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam1_dgam2[index];
                                        f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
                                    }
                                }
                            } else {
                                cout << "Error in Setup Interpolant Values for N3 122" << endl;
                                exit(0);
                            }
                            
                            
//                             f_N3_111_NON_GRAY[index] = (N3_111_NON_GRAY[index] - pow(N1_1_NON_GRAY[index], 3));
//                             f_N3_111_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2); 
//                             
//                             f_N3_122_NON_GRAY[index] = (N3_122_NON_GRAY[index] - N1_1_NON_GRAY[index]*pow(N1_2_NON_GRAY[index], 2));
//                             f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                             
//                             f_N3_123_NON_GRAY[index] = N3_123_NON_GRAY[index];
                            
                            // cout /*<< "index = " << index << "   "*/ << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "   " << "N3_111 = " << N3_111_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY[index] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index] << "  " << "dN3_122 = " << d2_N3_122_NON_GRAY_dgam2_dN1_1[index] << endl;
                            
                            //                                                     if (fabs(f_N3_111_NON_GRAY[index]) > 1.0 || fabs(f_N3_122_NON_GRAY[index]) > 1.0) {
//                             if (f_N3_111_NON_GRAY[index] != f_N3_111_NON_GRAY[index]) {
//                                 cout /*<< "index = " << index << "   "*/ << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "  " << "N3_111 = " << N3_111_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY[index] << "  " << "N3_122 = " << N3_122_NON_GRAY[index] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index] /*<< "  " << "dN3_111_dN1_1 = " << dN3_111_NON_GRAY_dN1_1[index]*/ /*<< "  " << "dN3_122_dN1_1 = " << dN3_122_NON_GRAY_dN1_1[index]*/ << endl;
//                                 exit(0);
//                             }
                        } // end for i_Cheby_gam1_gam2
                    } // end for i_Theta
                } // end for i_Phi
            } // end for i_Cheby_f
        } // end for i_Cheby_E
    } // end for id_Mobius
}

// void N3_Non_Gray_M2_3D_RT_Cheby :: SetupInterpolant_Values_BE_MPI() {
//     int index, index_temp;
//     long double norm_f, norm_f_2, gam3;
//     int id_E_min, id_E_max;
//     int id_f_min, id_f_max;
//     int id_theta_min, id_theta_max;
//     int id_phi_min, id_phi_max;
//     int size_rec = MPI_proc_parameters.size_rec;
//     
//     long double *f_N3_111_NON_GRAY_temp, *f_N3_122_NON_GRAY_temp, *f_N3_123_NON_GRAY_temp;
//     f_N3_111_NON_GRAY_temp = new long double[size_rec];
//     f_N3_122_NON_GRAY_temp = new long double[size_rec];
//     f_N3_123_NON_GRAY_temp = new long double[size_rec];
//     
//     Compute_MPI_Processes_Max_Min_Indexes(id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, id_theta_min, id_theta_max);
//     
//     for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
//         for (int i_Cheby_E = id_E_min; i_Cheby_E < id_E_max; i_Cheby_E++) {
//             for (int i_Cheby_f = id_f_min; i_Cheby_f < id_f_max; i_Cheby_f++) {
//                 for (int i_Phi = id_phi_min ; i_Phi < id_phi_max; i_Phi++) {
//                     for (int i_Theta = id_theta_min ; i_Theta < id_theta_max; i_Theta++) {
//                         for (int i_Cheby_gam1_gam2 = 0; i_Cheby_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_Cheby_gam1_gam2++) {
//                             
//                             index = Compute_Full_Id(id_Mobius, i_Cheby_E, i_Cheby_f, i_Phi, i_Theta, i_Cheby_gam1_gam2);
//                             
//                             index_temp = Compute_Full_Id_Single_proc(id_Mobius, i_Cheby_E, i_Cheby_f, i_Phi, i_Theta, i_Cheby_gam1_gam2);
//                             
// //                             if (id_proc == 0) {
// //                                 cout << "index = " << index << "   " << "index_temp = " << index_temp  << "   " << "size_rec = " << size_rec << endl;
// //                             }
//                             
//                             norm_f_2 = pow(N1_1_NON_GRAY[index], 2) + pow(N1_2_NON_GRAY[index], 2) + pow(N1_3_NON_GRAY[index], 2);
//                             
//                             // Compute f_N3_111
//                             if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(fabs(1.0 - gam1_NON_GRAY[index]) < 1.0e-8)) {
//                                 // gam1 > 0 and gam1 < 1
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 ///////////////////////////
//                                 exit(0);
//                                 ///////////////////////////////////
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = -d2_N3_111_NON_GRAY_dnorm_f_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = dN3_111_NON_GRAY_dN1_1[index]/(1.0 - norm_f_2);
//                                         
//                                         // cout << "dN3_111_NON_GRAY_dN1_1[index] = " << dN3_111_NON_GRAY_dN1_1[index] << endl;
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 ///////////////////////////
//                                 exit(0);
//                                 ///////////////////////////////////
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = (3.0/2.0)*pow(N1_1_NON_GRAY[index], 2) - dN3_111_NON_GRAY_dnorm_f[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = (N3_111_NON_GRAY[index] - pow(N1_1_NON_GRAY[index], 3));
//                                         f_N3_111_NON_GRAY_temp[index_temp] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_111_NON_GRAY_temp[index_temp] = f_N3_111_NON_GRAY_temp[index_temp] - gam1_NON_GRAY[index];
//                                 f_N3_111_NON_GRAY_temp[index_temp] /= gam1_NON_GRAY[index]*(1.0 - gam1_NON_GRAY[index]);
//                             } else if (gam1_NON_GRAY[index] < 1.0e-8) {
//                                 ///////////////////////////
//                                 exit(0);
//                                 ///////////////////////////////////
//                                 // gam1 = 0
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = -d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = d2_N3_111_NON_GRAY_dgam1_dN1_1[index]/(1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = -d2_N3_111_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = dN3_111_NON_GRAY_dgam1[index];
//                                         f_N3_111_NON_GRAY_temp[index_temp] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_111_NON_GRAY_temp[index_temp] = (f_N3_111_NON_GRAY_temp[index_temp] - 1.0)/(1.0 - gam1_NON_GRAY[index]);
//                             } else if (fabs(1.0 - gam1_NON_GRAY[index]) < 1.0e-8) {
//                                 ///////////////////////////
//                                 exit(0);
//                                 ///////////////////////////////////
//                                 // gam1 = 1
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = -d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;  
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = d2_N3_111_NON_GRAY_dgam1_dN1_1[index]/(1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = -d2_N3_111_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY_temp[index_temp] = dN3_111_NON_GRAY_dgam1[index];
//                                         f_N3_111_NON_GRAY_temp[index_temp] = N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_111_NON_GRAY_temp[index_temp] = (1.0 - f_N3_111_NON_GRAY_temp[index_temp])/gam1_NON_GRAY[index];
//                             } else {
//                                 cout << "Error in Setup Interpolant Values for N3 111" << endl;
//                                 exit(0);
//                             }
//                             
//                             // Compute f_N3_122
//                             if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 > 0 and gam1 < 1 and gam2 > 0 and gam2 < 1 - gam1
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = pow(N1_2_NON_GRAY[index], 2) - d2_N3_122_NON_GRAY_dnorm_f_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = dN3_122_NON_GRAY_dN1_1[index] - pow(N1_2_NON_GRAY[index], 2);
//                                         f_N3_122_NON_GRAY_temp[index_temp] /= (1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     // N1_1 != 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = (3.0/2.0)*pow(N1_2_NON_GRAY[index], 2) - dN3_122_NON_GRAY_dnorm_f[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = (N3_122_NON_GRAY[index] - N1_1_NON_GRAY[index]*pow(N1_2_NON_GRAY[index], 2));
//                                         f_N3_122_NON_GRAY_temp[index_temp] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_122_NON_GRAY_temp[index_temp] = f_N3_122_NON_GRAY_temp[index_temp] - gam2_NON_GRAY[index];
//                                 f_N3_122_NON_GRAY_temp[index_temp] /= gam1_NON_GRAY[index]*gam2_NON_GRAY[index];
//                             } else if ((gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 = 0 and gam2 > 0 and gam2 < 1 - gam1
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = -d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = d2_N3_122_NON_GRAY_dgam1_dN1_1[index];
//                                         f_N3_122_NON_GRAY_temp[index_temp] /= (1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = -d2_N3_122_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = dN3_122_NON_GRAY_dgam1[index];
//                                         f_N3_122_NON_GRAY_temp[index_temp] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_122_NON_GRAY_temp[index_temp] /= gam2_NON_GRAY[index];
//                             } else if (!(gam1_NON_GRAY[index] < 1.0e-8) && (gam2_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 > 0 and gam1 < 1 and gam2 = 0
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = -d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = d2_N3_122_NON_GRAY_dgam2_dN1_1[index];
//                                         f_N3_122_NON_GRAY_temp[index_temp] /= (1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = -d2_N3_122_NON_GRAY_dnorm_f_dgam2[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = dN3_122_NON_GRAY_dgam2[index];
//                                         f_N3_122_NON_GRAY_temp[index_temp] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_122_NON_GRAY_temp[index_temp] = (f_N3_122_NON_GRAY_temp[index_temp] - 1.0)/gam1_NON_GRAY[index];
//                             } else if ((gam1_NON_GRAY[index] < 1.0e-8) && (gam2_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 = 0 && gam2 = 0
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = -d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1[index];
//                                         f_N3_122_NON_GRAY_temp[index_temp] /= (1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = -d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY_temp[index_temp] = d2_N3_122_NON_GRAY_dgam1_dgam2[index];
//                                         f_N3_122_NON_GRAY_temp[index_temp] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                             } else {
//                                 cout << "Error in Setup Interpolant Values for N3 122" << endl;
//                                 exit(0);
//                             }
//                             
//                             // cout /*<< "index = " << index << "   "*/ << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "   " << "N3_111 = " << N3_111_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY_temp[index_temp] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index] << "  " << "dN3_122 = " << d2_N3_122_NON_GRAY_dgam2_dN1_1[index] << endl;
//                             
//                             //                                                     if (fabs(f_N3_111_NON_GRAY_temp[index_temp]) > 1.0 || fabs(f_N3_122_NON_GRAY[index]) > 1.0) {
//                             if (f_N3_111_NON_GRAY_temp[index_temp] != f_N3_111_NON_GRAY_temp[index_temp]) {
//                                 cout /*<< "index = " << index << "   "*/ << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY_temp[index_temp] << "  " << "N3_122 = " << N3_122_NON_GRAY[index] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index] << "  " << "dN3_111_dN1_1 = " << dN3_111_NON_GRAY_dN1_1[index] << "  " << "dN3_122_dN1_1 = " << dN3_122_NON_GRAY_dN1_1[index] << endl;
//                                 exit(0);
//                             }
//                         } // end for i_Cheby_gam1_gam2
//                     } // end for i_Theta
//                 } // end for i_Phi
//             } // end for i_Cheby_f
//         } // end for i_Cheby_E
//     } // end for id_Mobius
//     
// //     if (id_proc == 0) {
// //         cout  << id_E_min << "  " << id_E_max << "  " << id_f_min << "  " << id_f_max << "  " << id_theta_min << "  " << id_theta_max << "  " << id_phi_min << "  " << id_phi_max << endl;
// //         for (int i = 0; i < size_rec; i++) {
// //             cout << "i = " << i << "  " << "index_temp = " << index_temp << "  " << "N3_111_NON_GRAY = " << N3_111_NON_GRAY[i] << "  " << "f_N3_111_NON_GRAY_temp = " << f_N3_111_NON_GRAY_temp[i] << endl;
// //         }
// //     }
//     
//     MPI_Barrier(MPI_COMM_WORLD); 
//     
//     // Gather all maximum entropy solutions to primary processor 
//     MPI_Gather(f_N3_111_NON_GRAY_temp, size_rec, MPI_LONG_DOUBLE, f_N3_111_NON_GRAY, size_rec, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
//     MPI_Gather(f_N3_122_NON_GRAY_temp, size_rec, MPI_LONG_DOUBLE, f_N3_122_NON_GRAY, size_rec, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
//     // MPI_Gather(f_N3_123_NON_GRAY_temp, size_rec, MPI_LONG_DOUBLE, f_N3_123_NON_GRAY, size_rec, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
//     
//     MPI_Barrier(MPI_COMM_WORLD);
//     
//     if (id_proc == PRIMARY_ID) {
//         Reorder_Maximum_Entropy_Data_f_N3_ijk_MPI(f_N3_111_NON_GRAY, f_N3_122_NON_GRAY, f_N3_123_NON_GRAY);
//     }
//     
//     delete[] f_N3_111_NON_GRAY_temp; f_N3_111_NON_GRAY_temp = NULL;
//     delete[] f_N3_122_NON_GRAY_temp; f_N3_122_NON_GRAY_temp = NULL;
//     // delete[] f_N3_123_NON_GRAY_temp; f_N3_123_NON_GRAY_temp = NULL;
// }


void N3_Non_Gray_M2_3D_RT_Cheby :: Setup_Vandermonde_Matrix_Spherical_Harmonic() {
    long double *VanderMonde_Matrix, *VanderMonde_Vector;
    long double poly_SH;
    int index_full, i_Cheby_gam1_gam2;
    VanderMonde_Matrix = new long double[N_Coeffs_SH*N_Coeffs_SH];
    VanderMonde_Vector = new long double[N_Coeffs_SH];
    
    for (int index_SH_lin = 0; index_SH_lin < N_Coeffs_SH; index_SH_lin++) {
        index_full = index_SH_lin*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2;
        
        for (int index_SH_col = 0; index_SH_col < N_Coeffs_SH; index_SH_col++) {
            poly_SH = Cartesian_Spherical_harmonics(x_SH[index_full], y_SH[index_full], z_SH[index_full], Array_l_SH[index_SH_col], Array_m_SH[index_SH_col]);
            VanderMonde_Matrix[index_SH_lin*N_Coeffs_SH+index_SH_col] = poly_SH;     
        } // end for index_SH_col    
    } // end for i_lin_Cheby_E
    
    for (int index_SH_lin = 0; index_SH_lin < N_Coeffs_SH; index_SH_lin++) {
        for (int index_SH_col = 0; index_SH_col < N_Coeffs_SH; index_SH_col++) {
            VanderMonde_Vector[index_SH_col] = VanderMonde_Vector_N_vars(index_SH_col, index_SH_lin);
            //                                                     cout << "VanderMonde_Vector[i_Coeff_Cheby] = " << VanderMonde_Vector[i_Coeff_Cheby] << endl;
        } // end for index_SH_col
        
        //                                                 cout  << "*************************** solving the Vandermonde System ************************"<< endl;
        Solve_A_x_b(VanderMonde_Matrix, Coefficients_Vander_Matrix_SH, VanderMonde_Vector, N_Coeffs_SH, index_SH_lin);
    } // end for index_SH_lin
    
    cout << "Checking Vandermonde System for Spherical Harmonic !!!!!!!!!!!!!!!!!!!!" << endl;
    Check_A_x_b(VanderMonde_Matrix, Coefficients_Vander_Matrix_SH, N_Coeffs_SH);
    
    delete[] VanderMonde_Matrix;
    delete[] VanderMonde_Vector;
}


void N3_Non_Gray_M2_3D_RT_Cheby :: Vandermonde_Interpolation_Spherical_Harmonic() {
    
    int N_Points_NON_GRAY = Compute_Num_Coeffs_Full_Without_SH();
    
    for (int i_Cheby = 0; i_Cheby < N_Points_NON_GRAY; i_Cheby++) {
        for (int i_SH = 0; i_SH < N_Coeffs_SH; i_SH++) {
            Coeff_Lebedev_Quadrature_N3_111[i_SH*N_Points_NON_GRAY + i_Cheby] = 0.0;
            Coeff_Lebedev_Quadrature_N3_122[i_SH*N_Points_NON_GRAY + i_Cheby] = 0.0;
            Coeff_Lebedev_Quadrature_N3_123[i_SH*N_Points_NON_GRAY + i_Cheby] = 0.0;
            
            for (int j_SH = 0; j_SH < N_Coeffs_SH; j_SH++) {
                Coeff_Lebedev_Quadrature_N3_111[i_SH*N_Points_NON_GRAY + i_Cheby] += Coefficients_Vander_Matrix_SH[i_SH*N_Points_NON_GRAY+j_SH]*f_N3_111_NON_GRAY[i_SH*N_Points_NON_GRAY + j_SH];
                
                Coeff_Lebedev_Quadrature_N3_122[i_SH*N_Points_NON_GRAY + i_Cheby] += Coefficients_Vander_Matrix_SH[i_SH*N_Points_NON_GRAY+j_SH]*f_N3_122_NON_GRAY[i_SH*N_Points_NON_GRAY + j_SH];
                
                Coeff_Lebedev_Quadrature_N3_123[i_SH*N_Points_NON_GRAY + i_Cheby] += Coefficients_Vander_Matrix_SH[i_SH*N_Points_NON_GRAY+j_SH]*f_N3_123_NON_GRAY[i_SH*N_Points_NON_GRAY + j_SH];
            }
        }
    }
    
    
//     for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
//         for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
//             for (int i_Cheby_gam1_gam2 = 0; i_Cheby_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_Cheby_gam1_gam2++) {
//                 index_Cheby = (i_Cheby_E*N_Points_f + i_Cheby_f)*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2;
//                 
//                 for (int i_Phi = 0; i_Phi < N_Points_Phi; i_Phi++) {
//                     for (int i_Theta = 0; i_Theta < N_Points_Theta; i_Theta++) {
//                         index = i_Cheby_E*N_Points_f + i_Cheby_f;
//                         index = (index*N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
//                         index = index*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2; 
//                         
//                         id_Cheby_Lebed_rule = i_Phi*N_Points_Theta + i_Theta;
//                         
//                         Coeff_Vec_Lebedev_N3_111[id_Cheby_Lebed_rule] = f_N3_111_NON_GRAY[index];
//                         Coeff_Vec_Lebedev_N3_122[id_Cheby_Lebed_rule] = f_N3_122_NON_GRAY[index];
//                         Coeff_Vec_Lebedev_N3_123[id_Cheby_Lebed_rule] = f_N3_123_NON_GRAY[index];
//                     }
//                 }
//             }
//         }
//     }
}


void N3_Non_Gray_M2_3D_RT_Cheby :: Lebedev_Quadrature_Interpolation() {
    int index_Cheby, index;
    long double *Coeff_Vec_Lebedev_N3_111;
    long double *Coeff_Vec_Lebedev_N3_122;
    long double *Coeff_Vec_Lebedev_N3_123;
    long double *x_SH_temp, *y_SH_temp, *z_SH_temp;
    int N_Points_NON_GRAY_without_Lebed;
    int index_SH, Full_Id_SH;
    long double temp_val;
    long double error_SH = 0.0;
    int id_Cheby_Lebed_rule;
    
    Coeff_Vec_Lebedev_N3_111 = new long double [N_Points_Phi*N_Points_Theta];
    Coeff_Vec_Lebedev_N3_122 = new long double [N_Points_Phi*N_Points_Theta];
    Coeff_Vec_Lebedev_N3_123 = new long double [N_Points_Phi*N_Points_Theta];
    
    x_SH_temp = new long double [N_Points_Phi*N_Points_Theta];
    y_SH_temp = new long double [N_Points_Phi*N_Points_Theta];
    z_SH_temp = new long double [N_Points_Phi*N_Points_Theta];
    
    N_Points_NON_GRAY_without_Lebed = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
    
    // Orthogonality_Test(2);
    
    // Lebedev quadrature
    for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
        for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
            for (int i_Cheby_gam1_gam2 = 0; i_Cheby_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_Cheby_gam1_gam2++) {
                index_Cheby = (i_Cheby_E*N_Points_f + i_Cheby_f)*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2;
                    
                for (int i_Phi = 0; i_Phi < N_Points_Phi; i_Phi++) {
                    for (int i_Theta = 0; i_Theta < N_Points_Theta; i_Theta++) {
                        index = i_Cheby_E*N_Points_f + i_Cheby_f;
                        index = (index*N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
                        index = index*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2; 
                        
                        id_Cheby_Lebed_rule = i_Phi*N_Points_Theta + i_Theta;
                        
//                         // cout << "index = " << index << "  " << "x_SH = " << x_SH[index]  << "  " << "y_SH = " << y_SH[index]  << "  " << "z_SH = " << z_SH[index] << endl;
//                         
//                         cout << "index = " << index << "  " << "I0 = " << E_NON_GRAY[index]  << "  " << "N1_1 = " << N1_1_NON_GRAY[index]  << "  " << "N1_2 = " << N1_2_NON_GRAY[index]  << "  " << "N1_3 = " << N1_3_NON_GRAY[index] << "  " << "gam1 = " << gam1_NON_GRAY[index] << "  " << "gam2 = " << gam2_NON_GRAY[index]  << "  " << "N3_111 = " << N3_111_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY[index] /*<< "  " << "d_N3_111_d_N1_1 = " << dN3_111_NON_GRAY_dN1_1[index]*/ << endl;
//                         
// //                         Coeff_Vec_Lebedev_N3_111[id_Cheby_Lebed_rule] = N3_111_NON_GRAY[index];
// //                         Coeff_Vec_Lebedev_N3_122[id_Cheby_Lebed_rule] = N3_122_NON_GRAY[index];
// //                         Coeff_Vec_Lebedev_N3_123[id_Cheby_Lebed_rule] = N3_123_NON_GRAY[index];
                        
                        Coeff_Vec_Lebedev_N3_111[id_Cheby_Lebed_rule] = f_N3_111_NON_GRAY[index];
                        Coeff_Vec_Lebedev_N3_122[id_Cheby_Lebed_rule] = f_N3_122_NON_GRAY[index];
                        Coeff_Vec_Lebedev_N3_123[id_Cheby_Lebed_rule] = f_N3_123_NON_GRAY[index];
                        
                        x_SH_temp[id_Cheby_Lebed_rule] = x_SH[index];
                        y_SH_temp[id_Cheby_Lebed_rule] = y_SH[index];
                        z_SH_temp[id_Cheby_Lebed_rule] = z_SH[index];
                    }
                }
                
                for (int index_SH = 0; index_SH < N_Coeffs_SH; index_SH++) {
                    Full_Id_SH = index_SH*N_Points_NON_GRAY_without_Lebed+index_Cheby;
                    
                    Coeff_Lebedev_Quadrature_N3_111[Full_Id_SH] = Lebedev_Quadrature_Matrix_SH_Temp(Coeff_Vec_Lebedev_N3_111, Array_l_SH[index_SH], Array_m_SH[index_SH], N_Points_Phi, N_Points_Theta, x_SH_temp, y_SH_temp, z_SH_temp);
                     
                    Coeff_Lebedev_Quadrature_N3_122[Full_Id_SH] = Lebedev_Quadrature_Matrix_SH_Temp(Coeff_Vec_Lebedev_N3_122, Array_l_SH[index_SH], Array_m_SH[index_SH], N_Points_Phi, N_Points_Theta, x_SH_temp, y_SH_temp, z_SH_temp);

//                     Coeff_Lebedev_Quadrature_N3_123[Full_Id_SH] = Lebedev_Quadrature_Matrix_SH_Temp(Coeff_Vec_Lebedev_N3_123, Array_l_SH[index_SH], Array_m_SH[index_SH], N_Points_Phi, N_Points_Theta, x_SH, y_SH, z_SH);
                            
//                     cout << "Array_l_SH[index_SH] = " << Array_l_SH[index_SH] << "    " << "Array_m_SH[index_SH] = " << Array_m_SH[index_SH] << "    " << "Coeff_Lebedev_Quadrature_N3_111 = " << Coeff_Lebedev_Quadrature_N3_111[Full_Id_SH] << endl;
                }
                
                error_SH = 0.0;
                for (int i_Phi = 0; i_Phi < N_Points_Phi; i_Phi++) {
                    for (int i_Theta = 0; i_Theta < N_Points_Theta; i_Theta++) {
                        index = i_Cheby_E*N_Points_f + i_Cheby_f;
                        index = (index*N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
                        index = index*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2; 
                        
                        id_Cheby_Lebed_rule = i_Phi*N_Points_Theta + i_Theta;
                        
                        // temp_val = Test_Spherical_Harmonics_Even_Odd(Coeff_Lebedev_Quadrature_N3_111, index_Cheby, N_Points_NON_GRAY_without_Lebed, x_SH[index], y_SH[index], z_SH[index], Order_SH);
                        temp_val = Test_Spherical_Harmonics(Coeff_Lebedev_Quadrature_N3_111, index_Cheby, N_Points_NON_GRAY_without_Lebed, x_SH[index], y_SH[index], z_SH[index], Order_SH);
                        
                        error_SH = max(error_SH, fabs(temp_val - Coeff_Vec_Lebedev_N3_111[id_Cheby_Lebed_rule]));
                        
//                         if (fabs(temp_val - Coeff_Vec_Lebedev_N3_111[id_Cheby_Lebed_rule]) > 1.0e0) {
//                             cout << "id_Cheby_Lebed_rule = " << id_Cheby_Lebed_rule << "   " << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "    " << "Coeff_Vec_Lebedev_N3_111 = " << Coeff_Vec_Lebedev_N3_111[id_Cheby_Lebed_rule] << "     " << "temp_val = " << temp_val << "     " << "diff = " << temp_val - Coeff_Vec_Lebedev_N3_111[id_Cheby_Lebed_rule] << endl;
// //                             exit(0);
//                         }
                    }
                }
                // cout << "error_SH = " << error_SH << endl;
            } // end for i_Cheby_gam1_gam2
        } // end for i_Cheby_f 
    } // end for i_Cheby_E 
    
    delete[] Coeff_Vec_Lebedev_N3_111;
    delete[] Coeff_Vec_Lebedev_N3_122;
    delete[] Coeff_Vec_Lebedev_N3_123;
    
    delete[] x_SH_temp;
    delete[] y_SH_temp;
    delete[] z_SH_temp;
}

//********************************************************************************************************
// This routine sets up and solves the Vandermonde system of equations for the interpolative-based 
// approximation of the closing fluxes for either the gray or the non-gray M2 closure
//********************************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Setup_Vandermonde_Matrix_NG_N3_ijk() {
    int N_Points_NON_GRAY;
    long double ratio_E, norm_f;
    int index_full;
    long double *VanderMonde_Matrix, *VanderMonde_Vector;
    int iter_triangle_points;
    
    long double poly_map_I0_star, poly_norm_f, poly_gam1_gam2;
    
    N_Points_NON_GRAY = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
    VanderMonde_Matrix = new long double[N_Points_NON_GRAY*N_Points_NON_GRAY];
    VanderMonde_Vector = new long double[N_Points_NON_GRAY];
    
    int index_line, index_column;
    
    // Chebyshev quadrature
    for (int i_lin_Cheby_E = 0; i_lin_Cheby_E < N_Points_E; i_lin_Cheby_E++) {
        for (int i_lin_Cheby_f = 0; i_lin_Cheby_f < N_Points_f; i_lin_Cheby_f++) {
            for (int i_lin_Cheby_gam1_gam2 = 0; i_lin_Cheby_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_lin_Cheby_gam1_gam2++) {
                index_line = i_lin_Cheby_E*N_Points_f + i_lin_Cheby_f;
                index_line = index_line*N_Points_Triangle_gam1_gam2 + i_lin_Cheby_gam1_gam2; 
                    
                index_full = Compute_Full_Id(0, i_lin_Cheby_E, i_lin_Cheby_f, 0, 0, i_lin_Cheby_gam1_gam2);
                // index_full = (i_lin_Cheby_E*N_Points_f + i_lin_Cheby_f)*N_Points_Phi*N_Points_Theta;
                // index_full = index_full*N_Points_Triangle_gam1_gam2 + i_lin_Cheby_gam1_gam2; 
                
                ratio_E = E_NON_GRAY[index_full];
                
                norm_f = pow(N1_1_NON_GRAY[index_full],2) + pow(N1_2_NON_GRAY[index_full],2) + pow(N1_3_NON_GRAY[index_full],2);
                
                norm_f = sqrt(norm_f);
                
                // cout << "index_full = " << index_full << "  " << "ratio_E = " << ratio_E << "  " << "norm_f = " << norm_f << "  " << "gam1 = " << gam1_NON_GRAY[index_full] << "  " << "gam2 = " << gam2_NON_GRAY[index_full] << endl;
                
                for (int i_col_Cheby_E = 0; i_col_Cheby_E < N_Points_E; i_col_Cheby_E++) {
                    for (int i_col_Cheby_f = 0; i_col_Cheby_f < N_Points_f; i_col_Cheby_f++) {
                        iter_triangle_points = 0;
                        for (int i_col_Cheby_gam1 = 0; i_col_Cheby_gam1 < N_Points_gam1; i_col_Cheby_gam1++) {
                            for (int i_col_Cheby_gam2 = 0; i_col_Cheby_gam2 < N_Points_gam1 - i_col_Cheby_gam1; i_col_Cheby_gam2++) {
                                index_column = i_col_Cheby_E*N_Points_f + i_col_Cheby_f;
                                index_column = index_column*N_Points_Triangle_gam1_gam2 + iter_triangle_points; 
                                iter_triangle_points++;
                                
                                poly_map_I0_star = Chebyshev_Polynomial_Basis(ratio_E, i_col_Cheby_E);
                                poly_norm_f = Chebyshev_Polynomial_Basis(norm_f, 2*i_col_Cheby_f);
                                poly_gam1_gam2 = Proriol_Polynomials(gam1_NON_GRAY[index_full], gam2_NON_GRAY[index_full], i_col_Cheby_gam1, i_col_Cheby_gam2);
                                
                                // cout << "poly_map_I0_star = " << poly_map_I0_star << endl;
                                
                                VanderMonde_Matrix[index_line*N_Points_NON_GRAY+index_column] = poly_map_I0_star*poly_norm_f*poly_gam1_gam2;
                            } // end for i_col_Cheby_gam2
                        } // end for i_col_Cheby_gam1
                        
                        if (iter_triangle_points != N_Points_Triangle_gam1_gam2) {
                            cout << "iter_triangle_points = " << iter_triangle_points << "  "  << "N_Points_Triangle_gam1_gam2 = " << N_Points_Triangle_gam1_gam2 << endl;
                            exit(0);
                        }
                        
                    } // end for i_col_Cheby_f    
                } // end for i_col_Cheby_E
                
            } // end for i_lin_Cheby_gam1_gam2
        } // end for i_lin_Cheby_f
    } // end for i_lin_Cheby_E
    
    for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
        for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
            for (int i_Cheby_gam1_gam2 = 0; i_Cheby_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_Cheby_gam1_gam2++) {
                index_full = (i_Cheby_E*N_Points_f + i_Cheby_f)*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2; 
                
                for (int i_Coeff_Cheby = 0; i_Coeff_Cheby < N_Points_NON_GRAY; i_Coeff_Cheby++) {
                    VanderMonde_Vector[i_Coeff_Cheby] = VanderMonde_Vector_N_vars(i_Coeff_Cheby, index_full);
                    //                                                     cout << "VanderMonde_Vector[i_Coeff_Cheby] = " << VanderMonde_Vector[i_Coeff_Cheby] << endl;
                } // end for i_Coeff_Cheby
                
                //                                                 cout  << "*************************** solving the Vandermonde System ************************"<< endl;
                Solve_A_x_b(VanderMonde_Matrix, Coefficients_Vander_Matrix_NG_N3_ijk, VanderMonde_Vector, N_Points_NON_GRAY, index_full);
            } // end for i_Cheby_gam1_gam2
        } // end for i_Cheby_f
    } // end for i_Cheby_E
    
    // cout << "Checking Vandermonde System !!!!!!!!!!!!!!!!!!!!" << endl;
    // Check_A_x_b(VanderMonde_Matrix, Coefficients_Vander_Matrix_NG_N3_ijk, N_Points_NON_GRAY);
    
    delete[] VanderMonde_Matrix;
    delete[] VanderMonde_Vector;
}

//********************************************************************************************************
// This routine performs the Vandermonde interpolation for the closing fluxes for either the gray or 
// the non-gray M2 closure
//********************************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Vandermonde_Interpolation_NG_N3_ijk() {
    int N_Points_NON_GRAY = Compute_Num_Coeffs_Full_Without_SH();
    
    for (int i_SH = 0; i_SH < N_Coeffs_SH; i_SH++) {
        for (int i_Coeff_Cheby = 0; i_Coeff_Cheby < N_Points_NON_GRAY; i_Coeff_Cheby++) {
            Coefficient_Matrix_Fit_N3_111[i_SH*N_Points_NON_GRAY + i_Coeff_Cheby] = 0.0;
            Coefficient_Matrix_Fit_N3_122[i_SH*N_Points_NON_GRAY + i_Coeff_Cheby] = 0.0;
            Coefficient_Matrix_Fit_N3_123[i_SH*N_Points_NON_GRAY + i_Coeff_Cheby] = 0.0;
            
            for (int j_Coeff_Cheby = 0; j_Coeff_Cheby < N_Points_NON_GRAY; j_Coeff_Cheby++) {
//                 if (i_SH == 0) {
//                     cout << "Coefficients_Vander_Matrix_NG_N3_ijk = " << Coefficients_Vander_Matrix_NG_N3_ijk[i_Coeff_Cheby*N_Points_NON_GRAY+j_Coeff_Cheby] << endl;
//                 }
                
                Coefficient_Matrix_Fit_N3_111[i_SH*N_Points_NON_GRAY + i_Coeff_Cheby] += Coefficients_Vander_Matrix_NG_N3_ijk[i_Coeff_Cheby*N_Points_NON_GRAY+j_Coeff_Cheby]*Coeff_Lebedev_Quadrature_N3_111[i_SH*N_Points_NON_GRAY + j_Coeff_Cheby];
                
                Coefficient_Matrix_Fit_N3_122[i_SH*N_Points_NON_GRAY + i_Coeff_Cheby] += Coefficients_Vander_Matrix_NG_N3_ijk[i_Coeff_Cheby*N_Points_NON_GRAY+j_Coeff_Cheby]*Coeff_Lebedev_Quadrature_N3_122[i_SH*N_Points_NON_GRAY + j_Coeff_Cheby];
                
                Coefficient_Matrix_Fit_N3_123[i_SH*N_Points_NON_GRAY + i_Coeff_Cheby] += Coefficients_Vander_Matrix_NG_N3_ijk[i_Coeff_Cheby*N_Points_NON_GRAY+j_Coeff_Cheby]*Coeff_Lebedev_Quadrature_N3_123[i_SH*N_Points_NON_GRAY + j_Coeff_Cheby];
            }
            
//                 if (i_SH == 0) {
//                     cout << "i_Coeff_Cheby = " << i_Coeff_Cheby << "  " << "Coefficient_Matrix_Fit_N3_111[i_Coeff_Cheby] = " << Coefficient_Matrix_Fit_N3_111[i_SH*N_Points_NON_GRAY + i_Coeff_Cheby] << "  " << "Coeff_Lebedev_Quadrature_N3_111[i_Coeff_Cheby] = " << Coeff_Lebedev_Quadrature_N3_111[i_Coeff_Cheby] << endl;
//                 }
        }
    }
    
    
    int id_temp;
    for (int i_SH = 0; i_SH < N_Coeffs_SH; i_SH++) {
        for (int i_Coeff_Cheby = 0; i_Coeff_Cheby < N_Points_NON_GRAY; i_Coeff_Cheby++) {
            id_temp = i_SH*N_Points_NON_GRAY + i_Coeff_Cheby;
            Coefficient_Matrix_Fit_N3_111_Cheby_Basis[id_temp] = Coefficient_Matrix_Fit_N3_111[id_temp];
            Coefficient_Matrix_Fit_N3_122_Cheby_Basis[id_temp] = Coefficient_Matrix_Fit_N3_122[id_temp];
            Coefficient_Matrix_Fit_N3_123_Cheby_Basis[id_temp] = Coefficient_Matrix_Fit_N3_123[id_temp];
        }
    }
    
    // Recompute Coefficients for N3_111 in monomial basis
    Chebyshev_First_Kind_to_Monomial_Basis_ratio_E(Coefficient_Matrix_Fit_N3_111, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2);
    Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(Coefficient_Matrix_Fit_N3_111, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2);
    Proriol_to_Monomial_Basis(Coefficient_Matrix_Fit_N3_111, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2, N_Points_gam1);
    Spherical_Harmonics_to_Monomial_Basis(Coefficient_Matrix_Fit_N3_111, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2, Order_SH, Array_l_SH, Array_m_SH);
    
    // Recompute Coefficients for N3_122 in monomial basis
    Chebyshev_First_Kind_to_Monomial_Basis_ratio_E(Coefficient_Matrix_Fit_N3_122, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2);
    Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(Coefficient_Matrix_Fit_N3_122, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2);
    Proriol_to_Monomial_Basis(Coefficient_Matrix_Fit_N3_122, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2, N_Points_gam1);
    Spherical_Harmonics_to_Monomial_Basis(Coefficient_Matrix_Fit_N3_122, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2, Order_SH, Array_l_SH, Array_m_SH);
    
    // Recompute Coefficients for N3_123 in monomial basis
    Chebyshev_First_Kind_to_Monomial_Basis_ratio_E(Coefficient_Matrix_Fit_N3_123, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2);
    Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(Coefficient_Matrix_Fit_N3_123, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2);
    Proriol_to_Monomial_Basis(Coefficient_Matrix_Fit_N3_123, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2, N_Points_gam1);
    Spherical_Harmonics_to_Monomial_Basis(Coefficient_Matrix_Fit_N3_123, N_Points_E, N_Points_f, N_Coeffs_SH, N_Points_Triangle_gam1_gam2, Order_SH, Array_l_SH, Array_m_SH);
}

// ******************************************************************************************
// This routine computes the coefficients resulting from the Vandermonde interpolation of the 
// maximum-entropy solutions in the case of the gray M2 closure
// ******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Polynomial_Interpolation_Gray_M2_Closure(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, ofstream &output_Opt_Coefficients) {
    int size_rec = MPI_proc_parameters.size_rec;
    ofstream out_L_inf_Norm;
    
    Problem_Type = GRAY;
    
    Create_MPI_Data_Type_rec_N3();
    
    // Setup parameters for MPI-based parallel calculations
    Setup_MPI_Processes();
    
    // Precompute maximum-entropy solutions at the interpolations nodes for the purpose of the Vandermonde
    // interpolation
    iteration = LAST_ITERATION;
    Precompute_Final_Max_Ent_Solution();
    
//     OpenInputFile("Gray_M2_Model/Cheby_N3_Gray_M2_3D");
//     ReadInputData();
//     CloseInputFile();
    
    cout << "Done Mat entropy solution !!!!!!!!!!!!" << endl;
    
    // Set up the interpolant values based on the maximum-entropy solutions for the purpose
    // of the Vandermonde interpolation
    SetupInterpolant_Values_BE();
    
    cout << "Done Interpolant Setup!!!!!!!!!!!!" << endl;
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        // Perform intergration over the full solid angle to obtain coefficients for the spherical harmonics
        Lebedev_Quadrature_Interpolation();
        
        cout << "Done Spherical Harmonic Integration!!!!!!!!!!!!" << endl;
        
        // Set up and solve the Vandermonde system of equations for the purpose of our interpolation
        // procedure
        Setup_Vandermonde_Matrix_NG_N3_ijk();
        
        cout << "Done Setting up Vandermonde matrix!!!!!!!!!!!!" << endl;
        
        // Now perform the Vandermonde Interpolation based on the precomputed maximum-entropy solutions
        Vandermonde_Interpolation_NG_N3_ijk();
        
        cout << "Done Performing Vandermonde Interpolation !!!!!!!!!!!!" << endl;
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Broadcast Interpolative coefficients on all the other processors
    int N_Coeffs_Total = Compute_Num_Coeffs_Full();
    MPI_Bcast(Coefficient_Matrix_Fit_N3_111, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(Coefficient_Matrix_Fit_N3_122, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(Coefficient_Matrix_Fit_N3_123, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    MPI_Bcast(Coefficient_Matrix_Fit_N3_111_Cheby_Basis, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(Coefficient_Matrix_Fit_N3_122_Cheby_Basis, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(Coefficient_Matrix_Fit_N3_123_Cheby_Basis, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ Testing fit in 3D based on selected node for interpolation "<< endl;
    }
    
    // Ensure the Vandermonde interpolation was performed correctly by computing the error of our
    // approximation at the interpolation nodes
    // Note that the error at such points should be zero since the Vandermonde interpolation is exact at the 
    // interpolation nodes
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_111 "<< endl;
    }
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_111_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_112 "<< endl;
    }
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_112_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_122 "<< endl;
    }
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_122_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_123 "<< endl;
    }
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_123_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_222 "<< endl;
    }
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_222_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "****************** Testing fit in 3D based on selected node for interpolation Completed "<< endl;
        cout << endl;
    }
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ Testing accuracy of interpolation in 3D "<< endl;
    }
    
    // Now assess the accuracy of our interpolative-based approximation of the closing fluxes for the 
    // gray M2 closure based on a given set of evaluation points uniformly distributed throughout the 
    // realizable space for angular moments up to second-order
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ For N3_111 "<< endl;
//     }
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_111_ENTRY);
//     
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ For N3_112 "<< endl;
//     }
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_112_ENTRY);
//     
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ For N3_122 "<< endl;
//     }
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_122_ENTRY);
// //     if (id_proc == PRIMARY_ID) {
// //         cout  << "************************ For N3_123 "<< endl;
// //     }
//     // Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_123_ENTRY);
//     
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ For N3_222 "<< endl;
//     }
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_222_ENTRY);
//     
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ Testing accuracy of interpolation in 3D Completed "<< endl;
//     }
//     
//     // N3_3D_RT_Unif.L2_Norm_N3;
    
    
    if (id_proc == PRIMARY_ID) {
        Write_Coefficients_M2_Closure_Interp(output_Opt_Coefficients);
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
}

// ******************************************************************************
// This routine performs the polynomial interpolation for the third-order closing
// fluxes in either the Hyperbolic or the Logarithmic limit
// ******************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Polynomial_Interpolation_HL_LL(const int &Maximum_Entropy_Solution_Regime) {
    // First Precompute maximum entropy solutions
    Precompute_Final_Max_Ent_Solution(Maximum_Entropy_Solution_Regime);
    
    // Setup interpolant values for the purpose of polynomial interpolation
    // SetupInterpolant_Values_HL_LL();
    
    // Setup the Vandermonde matrix used to compute the interpolant
    Setup_Vandermonde_Matrix_NG_N3_ijk();
    
    // solve the Vandermonde system for the coefficients of the interpolant
    Vandermonde_Interpolation_NG_N3_ijk();
}

// ******************************************************************************
// This routine performs the polynomial interpolation for the third-order closing
// fluxes in the Bose Einstein Regime
// ******************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Polynomial_Interpolation_BE(N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_HL, N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_LL) {
    int size_rec = MPI_proc_parameters.size_rec;
    
//     OpenInputFile("Non_Gray_M2_Model/Cheby_N3_Non_Gray_M2_3D");
//     ReadInputData();
//     CloseInputFile();
    
    // First Precompute maximum entropy solutions
    Precompute_Final_Max_Ent_Solution();
    
    // Setup coefficients for the interpolation in both the Hyperbolic and Logarithmic limits
    // Setup_Coefficients_HL_LL(M2_3D_Data_N3_HL, M2_3D_Data_N3_LL);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        cout << "Done max entropy solution !!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    
    if (id_proc == PRIMARY_ID) {
        // Setup interpolant values for the purpose of polynomial interpolation
        SetupInterpolant_Values_BE();
        
        cout << "Done Setting up interpolant values !!!!!!!!!!!!!!!!!!!!!" << endl;
    
        // Perform intergration over the full solid angle to obtain coefficients for the spherical harmonics
        Lebedev_Quadrature_Interpolation();
        
        cout << "Done performing Lebedev quadrature !!!!!!!!!!!!!!!!!!!!!" << endl;
        
        // Setup the Vandermonde matrix used to compute the interpolant
        Setup_Vandermonde_Matrix_NG_N3_ijk();
        
        cout << "Done setting up Vandermonde matrix !!!!!!!!!!!!!!!!!!!!!" << endl;
        
        // solve the Vandermonde system for the coefficients of the interpolant
        Vandermonde_Interpolation_NG_N3_ijk();
        
        cout << "Done performing Vandermonde interpolation !!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Broadcast Interpolative coefficients on all the other processors
    int N_Coeffs_Total = Compute_Num_Coeffs_Full();
    MPI_Bcast(Coefficient_Matrix_Fit_N3_111, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(Coefficient_Matrix_Fit_N3_122, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(Coefficient_Matrix_Fit_N3_123, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    MPI_Bcast(Coefficient_Matrix_Fit_N3_111_Cheby_Basis, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(Coefficient_Matrix_Fit_N3_122_Cheby_Basis, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(Coefficient_Matrix_Fit_N3_123_Cheby_Basis, N_Coeffs_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Broadcast Interpolation data on all the other processors
    int N_Points_Total = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Theta*N_Points_Phi*N_Points_Triangle_gam1_gam2;
    MPI_Bcast(f_N3_111_NON_GRAY, N_Points_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(f_N3_122_NON_GRAY, N_Points_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Bcast(f_N3_123_NON_GRAY, N_Points_Total, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
}



/*inline*/ void N3_Non_Gray_M2_3D_RT_Cheby :: Copy_to(N3_Non_Gray_M2_3D_RT_Cheby *New_N3_M2_3D_RT_Cheby) {
    New_N3_M2_3D_RT_Cheby->Rec_Moms_Test.Copy(Rec_Moms_Test);
    
    int N_pts_total;
    New_N3_M2_3D_RT_Cheby->id_proc = id_proc;
    New_N3_M2_3D_RT_Cheby->num_proc = num_proc;
    New_N3_M2_3D_RT_Cheby->num_proc_used = num_proc_used;
    
    New_N3_M2_3D_RT_Cheby->num_proc_E = num_proc_E;
    New_N3_M2_3D_RT_Cheby->num_proc_f = num_proc_f;
    New_N3_M2_3D_RT_Cheby->num_proc_Phi = num_proc_Phi;
    New_N3_M2_3D_RT_Cheby->num_proc_Theta = num_proc_Theta;
    New_N3_M2_3D_RT_Cheby->num_proc_Triangle_gam1_gam2 = num_proc_Triangle_gam1_gam2;
    
    New_N3_M2_3D_RT_Cheby->MPI_proc_parameters.Copy(MPI_proc_parameters);
    
    New_N3_M2_3D_RT_Cheby->Max_Ent_Data_Type = Max_Ent_Data_Type;
    
    New_N3_M2_3D_RT_Cheby->index_f = index_f;
    New_N3_M2_3D_RT_Cheby->index_Phi = index_Phi;
    New_N3_M2_3D_RT_Cheby->index_Theta = index_Theta;
    New_N3_M2_3D_RT_Cheby->index_gam1 = index_gam1;
    New_N3_M2_3D_RT_Cheby->index_gam2 = index_gam2;
    New_N3_M2_3D_RT_Cheby->index_triangle = index_triangle;
    
    New_N3_M2_3D_RT_Cheby->N_pts_Mob_Scale = N_pts_Mob_Scale;
    New_N3_M2_3D_RT_Cheby->N_Points_E = N_Points_E;
    New_N3_M2_3D_RT_Cheby->N_Points_f = N_Points_f;
    New_N3_M2_3D_RT_Cheby->N_Points_gam1 = N_Points_gam1;
    New_N3_M2_3D_RT_Cheby->N_Points_Phi = N_Points_Phi;
    New_N3_M2_3D_RT_Cheby->N_Points_Theta = N_Points_Theta;
    New_N3_M2_3D_RT_Cheby->N_Points_Triangle_gam1_gam2 = N_Points_Triangle_gam1_gam2;
    
    New_N3_M2_3D_RT_Cheby->Order_SH = Order_SH;
    New_N3_M2_3D_RT_Cheby->N_Coeffs_SH = N_Coeffs_SH;
    
    New_N3_M2_3D_RT_Cheby->index_Ncoeffs = index_Ncoeffs;
    New_N3_M2_3D_RT_Cheby->index_rec_N3 = index_rec_N3;
    
    New_N3_M2_3D_RT_Cheby->Node_Distribution_E = Node_Distribution_E;
    New_N3_M2_3D_RT_Cheby->Node_Distribution_f = Node_Distribution_f;
    New_N3_M2_3D_RT_Cheby->Node_Distribution_Theta_Phi = Node_Distribution_Theta_Phi;
    New_N3_M2_3D_RT_Cheby->Node_Distribution_gam1 = Node_Distribution_gam1;
    
    // cout << "What was the issue when I was not reallocating???????????" << endl;
    // cout << "Most likely because I did not specify Order_SH before" << endl;
    New_N3_M2_3D_RT_Cheby->allocate();
    
    for (int i = 0; i < N_Coeffs_SH; i++) {
        New_N3_M2_3D_RT_Cheby->Array_l_SH[i] = Array_l_SH[i];
        New_N3_M2_3D_RT_Cheby->Array_m_SH[i] = Array_m_SH[i];
    }
    
    for (int i = 0; i < N_Coeffs_SH*N_Coeffs_SH; i++) {
        New_N3_M2_3D_RT_Cheby->Coefficients_Vander_Matrix_SH[i] = Coefficients_Vander_Matrix_SH[i];
    }
    
    for (int i = 0; i < N_Points_E*N_Points_E; i++) {
        New_N3_M2_3D_RT_Cheby->Coefficients_Vander_Matrix_Least_Squares_L_I0_star[i] = Coefficients_Vander_Matrix_Least_Squares_L_I0_star[i];
    }
    
    N_pts_total = N_Points_f*N_Points_Triangle_gam1_gam2;
    for (int i = 0; i < N_pts_total*N_pts_total; i++) {
        New_N3_M2_3D_RT_Cheby->Coefficients_Vander_Matrix_L_I0_star_N3_ijk[i] = Coefficients_Vander_Matrix_L_I0_star_N3_ijk[i];
    }
    
    N_pts_total = N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2;
    for (int i = 0; i < N_pts_total*N_pts_total; i++) {
        New_N3_M2_3D_RT_Cheby->Coefficients_Vander_Matrix_NG_N3_ijk[i] = Coefficients_Vander_Matrix_NG_N3_ijk[i];
    }
    
    N_pts_total = N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
    for (int i = 0; i < N_pts_total; i++) {
        New_N3_M2_3D_RT_Cheby->Coefficient_Mobius_Scale_N3_ijk[i] = Coefficient_Mobius_Scale_N3_ijk[i];
    }
    
//     N_pts_total = N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
//     for (int i = 0; i < N_pts_total; i++) {
//         New_N3_M2_3D_RT_Cheby->Optim_Length_Scale_N3_ijk[i] = Optim_Length_Scale_N3_ijk[i];
//     }
    
    N_pts_total = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
    for (int i = 0; i < N_pts_total; i++) {
        New_N3_M2_3D_RT_Cheby->E_NON_GRAY[i] = E_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->N1_1_NON_GRAY[i] = N1_1_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->N1_2_NON_GRAY[i] = N1_2_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->N1_3_NON_GRAY[i] = N1_3_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->gam1_NON_GRAY[i] = gam1_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->gam2_NON_GRAY[i] = gam2_NON_GRAY[i];
        
        New_N3_M2_3D_RT_Cheby->x_SH[i] = x_SH[i];
        New_N3_M2_3D_RT_Cheby->y_SH[i] = y_SH[i];
        New_N3_M2_3D_RT_Cheby->z_SH[i] = z_SH[i];
        
        New_N3_M2_3D_RT_Cheby->N3_111_NON_GRAY[i] = N3_111_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->N3_112_NON_GRAY[i] = N3_112_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->N3_122_NON_GRAY[i] = N3_122_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->N3_123_NON_GRAY[i] = N3_123_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->N3_222_NON_GRAY[i] = N3_222_NON_GRAY[i];
        
        New_N3_M2_3D_RT_Cheby->f_N3_111_NON_GRAY[i] = f_N3_111_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->f_N3_122_NON_GRAY[i] = f_N3_122_NON_GRAY[i];
        New_N3_M2_3D_RT_Cheby->f_N3_123_NON_GRAY[i] = f_N3_123_NON_GRAY[i];
        
        // Derivatives for N3_111
        New_N3_M2_3D_RT_Cheby->dN3_111_NON_GRAY_dN1_1[i] = dN3_111_NON_GRAY_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->dN3_111_NON_GRAY_dN1_2[i] = dN3_111_NON_GRAY_dN1_2[i];
        New_N3_M2_3D_RT_Cheby->dN3_111_NON_GRAY_dN1_3[i] = dN3_111_NON_GRAY_dN1_3[i];
        New_N3_M2_3D_RT_Cheby->dN3_111_NON_GRAY_dmu[i] = dN3_111_NON_GRAY_dmu[i];
        New_N3_M2_3D_RT_Cheby->dN3_111_NON_GRAY_dnorm_f[i] = dN3_111_NON_GRAY_dnorm_f[i];
        New_N3_M2_3D_RT_Cheby->dN3_111_NON_GRAY_dgam1[i] = dN3_111_NON_GRAY_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_111_NON_GRAY_dnorm_f_dN1_1[i] = d2_N3_111_NON_GRAY_dnorm_f_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_111_NON_GRAY_dnorm_f_dgam1[i] = d2_N3_111_NON_GRAY_dnorm_f_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_111_NON_GRAY_dgam1_dN1_1[i] = d2_N3_111_NON_GRAY_dgam1_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[i] = d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[i];
        
        // Derivatives for N3_122
        New_N3_M2_3D_RT_Cheby->dN3_122_NON_GRAY_dN1_1[i] = dN3_122_NON_GRAY_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->dN3_122_NON_GRAY_dN1_2[i] = dN3_122_NON_GRAY_dN1_2[i];
        New_N3_M2_3D_RT_Cheby->dN3_122_NON_GRAY_dN1_3[i] = dN3_122_NON_GRAY_dN1_3[i];
        New_N3_M2_3D_RT_Cheby->dN3_122_NON_GRAY_dmu[i] = dN3_122_NON_GRAY_dmu[i];
        New_N3_M2_3D_RT_Cheby->dN3_122_NON_GRAY_dnorm_f[i] = dN3_122_NON_GRAY_dnorm_f[i];
        New_N3_M2_3D_RT_Cheby->dN3_122_NON_GRAY_dgam1[i] = dN3_122_NON_GRAY_dgam1[i];
        New_N3_M2_3D_RT_Cheby->dN3_122_NON_GRAY_dgam2[i] = dN3_122_NON_GRAY_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_122_NON_GRAY_dgam1_dgam2[i] = d2_N3_122_NON_GRAY_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_122_NON_GRAY_dnorm_f_dN1_1[i] = d2_N3_122_NON_GRAY_dnorm_f_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_122_NON_GRAY_dnorm_f_dgam1[i] = d2_N3_122_NON_GRAY_dnorm_f_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_122_NON_GRAY_dnorm_f_dgam2[i] = d2_N3_122_NON_GRAY_dnorm_f_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_122_NON_GRAY_dgam1_dN1_1[i] = d2_N3_122_NON_GRAY_dgam1_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_122_NON_GRAY_dgam2_dN1_1[i] = d2_N3_122_NON_GRAY_dgam2_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1[i] = d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1[i] = d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2[i] = d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1[i] = d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1[i] = d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1[i];
        
        // Derivatives for N3_123
        New_N3_M2_3D_RT_Cheby->dN3_123_NON_GRAY_dN1_1[i] = dN3_123_NON_GRAY_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->dN3_123_NON_GRAY_dN1_2[i] = dN3_123_NON_GRAY_dN1_2[i];
        New_N3_M2_3D_RT_Cheby->dN3_123_NON_GRAY_dN1_3[i] = dN3_123_NON_GRAY_dN1_3[i];
        New_N3_M2_3D_RT_Cheby->dN3_123_NON_GRAY_dmu[i] = dN3_123_NON_GRAY_dmu[i];
        New_N3_M2_3D_RT_Cheby->dN3_123_NON_GRAY_dgam1[i] = dN3_123_NON_GRAY_dgam1[i];
        New_N3_M2_3D_RT_Cheby->dN3_123_NON_GRAY_dgam2[i] = dN3_123_NON_GRAY_dgam2[i];
        New_N3_M2_3D_RT_Cheby->dN3_123_NON_GRAY_dgam3[i] = dN3_123_NON_GRAY_dgam3[i];
        
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_1_dN1_2[i] = d2_N3_123_NON_GRAY_dN1_1_dN1_2[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_1_dN1_3[i] = d2_N3_123_NON_GRAY_dN1_1_dN1_3[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_2_dN1_3[i] = d2_N3_123_NON_GRAY_dN1_2_dN1_3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3[i] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3[i];
        
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dgam1_dN1_1[i] = d2_N3_123_NON_GRAY_dgam1_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dgam2_dN1_1[i] = d2_N3_123_NON_GRAY_dgam2_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_1_dgam3[i] = d2_N3_123_NON_GRAY_dN1_1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_2_dgam1[i] = d2_N3_123_NON_GRAY_dN1_2_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_2_dgam2[i] = d2_N3_123_NON_GRAY_dN1_2_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_2_dgam3[i] = d2_N3_123_NON_GRAY_dN1_2_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_3_dgam1[i] = d2_N3_123_NON_GRAY_dN1_3_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_3_dgam2[i] = d2_N3_123_NON_GRAY_dN1_3_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dN1_3_dgam3[i] = d2_N3_123_NON_GRAY_dN1_3_dgam3[i];
        
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1[i] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1[i] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1[i] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2[i] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2[i] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2[i] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3[i] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3[i] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3[i] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3[i];
        
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dgam1_dgam2[i] = d2_N3_123_NON_GRAY_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dgam1_dgam3[i] = d2_N3_123_NON_GRAY_dgam1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d2_N3_123_NON_GRAY_dgam2_dgam3[i] = d2_N3_123_NON_GRAY_dgam2_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1[i] = d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3[i] = d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3[i] = d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2[i] = d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3[i] = d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3[i] = d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2[i] = d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3[i] = d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3[i] = d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3[i];
        
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3[i] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2[i] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3[i] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3[i] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2[i] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2[i];
        New_N3_M2_3D_RT_Cheby->d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3[i] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3[i];
        New_N3_M2_3D_RT_Cheby->d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3[i] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3[i];
    }
    
//     N_pts_total = N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
//     for (int i = 0; i < N_pts_total; i++) {
//         New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_111_HL[i] = Coefficient_Matrix_Fit_N3_111_HL[i];
//         New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_111_LL[i] = Coefficient_Matrix_Fit_N3_111_LL[i];
        
//         New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_122_HL[i] = Coefficient_Matrix_Fit_N3_122_HL[i];
//         New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_122_LL[i] = Coefficient_Matrix_Fit_N3_122_LL[i];
        
//         New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_123_HL[i] = Coefficient_Matrix_Fit_N3_123_HL[i];
//         New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_123_LL[i] = Coefficient_Matrix_Fit_N3_123_LL[i];
//     }
    
    N_pts_total = N_Points_E*N_Points_f*N_Coeffs_SH*N_Points_Triangle_gam1_gam2;
    for (int i = 0; i < N_pts_total; i++) {
        New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_111[i] = Coefficient_Matrix_Fit_N3_111[i];
        New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_122[i] = Coefficient_Matrix_Fit_N3_122[i];
        New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_123[i] = Coefficient_Matrix_Fit_N3_123[i];
        
        New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_111_Cheby_Basis[i] = Coefficient_Matrix_Fit_N3_111_Cheby_Basis[i];
        New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_122_Cheby_Basis[i] = Coefficient_Matrix_Fit_N3_122_Cheby_Basis[i];
        New_N3_M2_3D_RT_Cheby->Coefficient_Matrix_Fit_N3_123_Cheby_Basis[i] = Coefficient_Matrix_Fit_N3_123_Cheby_Basis[i];
    }
}

// ******************************************************************************************
// This routine computes optimal values of the length scale f_L_N3_ijk of the algebraic 
// mapping of the first-order normalized angular moment for the purpose of our interpolative-
// based approximation of the closing fluxes of the non-gray M2 closure. The associated 
// coefficients required for the polynomial interpolation of both the optimal length scale
// L_N3_ijk and the closing fluxes N3_ijk are also computed alongside
// ******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Polynomial_Interpolation_Non_Gray_M2_N3_ijk(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, ofstream &output_Opt_Coefficients) {
    ofstream out_L_inf_Norm;
    
    Problem_Type = NON_GRAY;
    
    Create_MPI_Data_Type_rec_N3();
    
    // Setup parameters for MPI-based parallel calculations
    Setup_MPI_Processes();
    
    // First perform polynomial interpolation in the Hyperbolic and Free-Streaming limits
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_HL(1, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    N3_Non_Gray_M2_3D_RT_Cheby N3_M2_3D_RT_LL(1, N_Points_f, N_Points_Phi, N_Points_Theta, N_Points_gam1, N_pts_Mob_Scale, Order_SH);
    
    // Hyperbolic Limit
    // N3_M2_3D_RT_HL.Polynomial_Interpolation_HL_LL(HYPERBOLIC_LIMIT);
    // Now Logarithmic Limit
    // N3_M2_3D_RT_LL.Polynomial_Interpolation_HL_LL(LOGARITHMIC_LIMIT);
    
    // Compute optimal set of coefficients for the length scale of the exponential 
    // mapping of the radiative energy density
    Least_Squares_Optimization_Coefficients_L_I0_star(N3_3D_RT_Unif, N3_M2_3D_RT_HL, N3_M2_3D_RT_LL);
    
    // Coefficient_Mobius_Scale_N3_ijk[0] = 3.0;
    
    // Perform polynomial interpolation for the non-gray M2 closure, based on the optimal set
    // of coefficients for the length scale of the exponential mapping of the radiative energy
    // density, computed previously
    iteration = LAST_ITERATION;
    Polynomial_Interpolation_BE(N3_M2_3D_RT_HL, N3_M2_3D_RT_LL);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ Testing fit in 3D based on selected node for interpolation "<< endl;
    }
    
    // Ensure the Vandermonde interpolation was performed correctly by computing the error of our
    // approximation at the interpolation nodes
    // Note that the error at such points should be zero since the Vandermonde interpolation is exact at the 
    // interpolation nodes
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_111 "<< endl;
    }
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ In monomial basis "<< endl;
    }
    flag_basis_type = MONOMIAL_BASIS;
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_111_ENTRY);
    
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ In Chebyshev basis "<< endl;
//     }
//     flag_basis_type = CHEBYSHEV_BASIS;
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_111_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_122 "<< endl;
    }
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ In monomial basis "<< endl;
    }
    flag_basis_type = MONOMIAL_BASIS;
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_122_ENTRY);
    
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ In Chebyshev basis "<< endl;
//     }
//     flag_basis_type = CHEBYSHEV_BASIS;
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_122_ENTRY);
    
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ For N3_123 "<< endl;
//     }
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_123_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "****************** Testing fit in 3D based on selected node for interpolation Completed "<< endl;
        cout << endl;
    }
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ Testing accuracy of interpolation in 3D "<< endl;
    }
    
    // Now assess the accuracy of our interpolative-based approximation of the closing fluxes for the 
    // gray M2 closure based on a given set of evaluation points uniformly distributed throughout the 
    // realizable space for angular moments up to second-order
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_111 "<< endl;
    }
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ In monomial basis "<< endl;
    }
    flag_basis_type = MONOMIAL_BASIS;
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_111_ENTRY);
    
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ In Chebyshev basis "<< endl;
//     }
//     flag_basis_type = CHEBYSHEV_BASIS;
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_111_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ For N3_122 "<< endl;
    }
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ In monomial basis "<< endl;
    }
    flag_basis_type = MONOMIAL_BASIS;
    Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_122_ENTRY);
    
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ In Chebyshev basis "<< endl;
//     }
//     flag_basis_type = CHEBYSHEV_BASIS;
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_122_ENTRY);
    
//     if (id_proc == PRIMARY_ID) {
//         cout  << "************************ For N3_123 "<< endl;
//     }
//     Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_3D_RT_Unif, 1, out_L_inf_Norm, N3_123_ENTRY);
    
    if (id_proc == PRIMARY_ID) {
        cout  << "************************ Testing accuracy of interpolation in 3D Completed "<< endl;
    }
        
    if (id_proc == PRIMARY_ID) {
        Write_Coefficients_M2_Closure_Interp(output_Opt_Coefficients);
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Write_Coefficients_M2_Closure_Interp(ofstream &output_Opt_Coefficients) {
    int index_Coeffs_M2_Closure, index_Coeffs_M2_Closure_max;
    int VAR_NUM;
    
    if (Problem_Type == NON_GRAY) {
        // output_Opt_Coefficients << Coefficient_Mobius_Scale_N3_ijk[0];
        // Output the coefficients for the polynomial interpolation approximation of the
        // optimal length scale in the case of non-gray M2
        index_Coeffs_M2_Closure_max = N_Coeffs_SH*N_Points_f*N_Points_Triangle_gam1_gam2 - 1;
        for (int i_SH = 0; i_SH < 1/*N_Coeffs_SH*/; i_SH++) {
            for (int i_fit_f = 0; i_fit_f < N_Points_f; i_fit_f++) {
                for (int i_fit_gam1_gam2 = 0; i_fit_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_fit_gam1_gam2++) {
                    index_Coeffs_M2_Closure = i_SH*N_Points_f + i_fit_f;
                    index_Coeffs_M2_Closure = index_Coeffs_M2_Closure*N_Points_Triangle_gam1_gam2 + i_fit_gam1_gam2;
                    
                    output_Opt_Coefficients << Coefficient_Mobius_Scale_N3_ijk[index_Coeffs_M2_Closure];
                    
                    if (index_Coeffs_M2_Closure < index_Coeffs_M2_Closure_max) {
                        output_Opt_Coefficients << setw(16);
                    } // end if
                } // end for i_fit_gam1_gam2
            } // end for i_fit_f
        } // end for i_SH
        output_Opt_Coefficients << endl;
    } // end if
        
    // Output the coefficients for the polynomial interpolation approximation of the
    // closing fluxes for the gray or non-gray M2 closure
    index_Coeffs_M2_Closure_max = N_Coeffs_SH*N_Points_E*N_Points_f*N_Points_Triangle_gam1_gam2 - 1;
    
    for (int id_N3_ijk = 0; id_N3_ijk < 3; id_N3_ijk++) {
        VAR_NUM = N3_ijk_Entry_Index[id_N3_ijk];
        for (int i_SH = 0; i_SH < N_Coeffs_SH; i_SH++) {
            for (int i_fit_E = 0; i_fit_E < N_Points_E; i_fit_E++) {
                for (int i_fit_f = 0; i_fit_f < N_Points_f; i_fit_f++) {
                    for (int i_fit_gam1_gam2 = 0; i_fit_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_fit_gam1_gam2++) {
                        index_Coeffs_M2_Closure = (i_SH*N_Points_E + i_fit_E)*N_Points_f + i_fit_f;
                        index_Coeffs_M2_Closure = index_Coeffs_M2_Closure*N_Points_Triangle_gam1_gam2 + i_fit_gam1_gam2;
                        
                        switch (VAR_NUM) {
                            case N3_111_ENTRY:
                                output_Opt_Coefficients << setprecision(12) << Coefficient_Matrix_Fit_N3_111[index_Coeffs_M2_Closure];
                                break;
                            case N3_122_ENTRY:
                                output_Opt_Coefficients << setprecision(12) << Coefficient_Matrix_Fit_N3_122[index_Coeffs_M2_Closure];
                                break;
                            case N3_123_ENTRY:
                                output_Opt_Coefficients << setprecision(12) << Coefficient_Matrix_Fit_N3_123[index_Coeffs_M2_Closure];
                                break;
                            default:
                                cout << "Invalid value for VAR_NUM" << endl;
                                exit(0);
                                break;
                        }
                        
                        if (index_Coeffs_M2_Closure < index_Coeffs_M2_Closure_max) {
                            output_Opt_Coefficients << setw(20);
                        } else {
                            output_Opt_Coefficients << endl;
                        }// end for if
                    } // end for i_fit_gam1_gam2
                } // end for i_fit_f
            } // end for i_fit_E
        } // end for i_SH
    } // end for id_N3_ijk
    
    // cout << "index_Coeffs_M2_Closure = " << index_Coeffs_M2_Closure << "  " << "index_Coeffs_M2_Closure_max = " << index_Coeffs_M2_Closure_max << endl;
}

// // ******************************************************************************************
// // Vandermonde Interpolation Check
// // ******************************************************************************************
// void N3_Non_Gray_M2_3D_RT_Cheby :: Vandermonde_Interpolation_Check(const int &VAR_NUM) {
//     long double max_err_N3_N1_1, max_err_N3_N1_2, max_err_N3_N1_3, max_err_N3_gam1, max_err_N3_gam2;
//     long double L2_Norm_N3, L_inf_Norm_N3;
//     long double ratio_E, N1_1, N1_2, N1_3, gam1, gam2;
//     long double N3_Fit, N3_Numerical;
//     int index;
//     int id_E_min, id_E_max;
//     int id_f_min, id_f_max;
//     int id_theta_min, id_theta_max;
//     int id_phi_min, id_phi_max;
//     long double *L2_Norm_N3_array = NULL, *L_inf_Norm_N3_array = NULL;
//     
//     if (id_proc == PRIMARY_ID) {
//         L2_Norm_N3_array = new long double[num_proc];
//         L_inf_Norm_N3_array = new long double[num_proc];
//     }
//     
//     L2_Norm_N3 = 0.0;
//     L_inf_Norm_N3 = 0.0;
//     
//     Compute_MPI_Processes_Max_Min_Indexes(id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, id_theta_min, id_theta_max);
//     
//     for (int i_E = id_E_min; i_E < id_E_max; i_E++) {
//         for (int i_f = id_f_min; i_f < id_f_max; i_f++) {
//             for (int i_Phi = id_phi_min; i_Phi < id_phi_max; i_Phi++) {
//                 for (int i_Theta = id_theta_min; i_Theta < id_theta_max; i_Theta++) {
//                     for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_gam1_gam2++) {
//                         index = Compute_Full_Id(0, i_E, i_f, i_Phi, i_Theta, i_gam1_gam2);
//                         
//                         // index_temp = Compute_Full_Id_Single_proc(0, i_E, i_f, i_Phi, i_Theta, i_gam1_gam2);
//                         
//                         ratio_E = E_NON_GRAY[index];
//                         N1_1 = N1_1_NON_GRAY[index];
//                         N1_2 = N1_2_NON_GRAY[index];
//                         N1_3 = N1_3_NON_GRAY[index];
//                         gam1 = gam1_NON_GRAY[index];
//                         gam2 = gam2_NON_GRAY[index];
//                         
//                         if (VAR_NUM == N3_111_ENTRY) {
//                             // N3_Numerical = f_N3_111_NON_GRAY[index];
//                             // N3_Fit = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
//                             N3_Numerical = Coeff_Lebedev_Quadrature_N3_111[index];
//                             N3_Fit = Evaluate_N3_111(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
//                         } else if (VAR_NUM == N3_122_ENTRY) {
//                             // N3_Numerical = f_N3_122_NON_GRAY[index];
//                             // N3_Fit = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
//                             N3_Numerical = N3_122_NON_GRAY[index];
//                             N3_Fit = Evaluate_N3_122(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
//                         } else if (VAR_NUM == N3_123_ENTRY) {
//                             N3_Numerical = N3_123_NON_GRAY[index];
//                             N3_Fit = Evaluate_N3_123(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
//                         } else {
//                             cout << "********************* VAR_NUM not defined !!!!! ************************" << endl;
//                             exit(0);
//                         }
//                         
//                         L2_Norm_N3 += pow(N3_Fit - N3_Numerical, 2);
//                         L_inf_Norm_N3 = max(L_inf_Norm_N3, fabs(N3_Fit - N3_Numerical));
//                         
// //                         if (id_proc == PRIMARY_ID) {
// //                             cout << "index = " << index << "   " << "ratio_E = " << ratio_E << "   " << "N1_1 = " << N1_1 << "   " << "N1_2 = " << N1_2 << "   " << "N1_3 = " << N1_3 << "   " << "gam1 = " << gam1 << "   " << "gam2 = " << gam2 << "   " << "N3_Fit = " << N3_Fit << "    " << "N3_Numerical = " << N3_Numerical << "    " << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << endl;
// //                         }
//                         
//                         if (L_inf_Norm_N3 == fabs(N3_Fit - N3_Numerical)) {
//                             max_err_N3_N1_1 = N1_1;
//                             max_err_N3_N1_2 = N1_2;
//                             max_err_N3_N1_3 = N1_3;
//                             max_err_N3_gam1 = gam1;
//                             max_err_N3_gam2 = gam2;
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     
//     // cout << "id_proc = " << id_proc << "  " << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << endl; 
//     
//     // Wait for all the processors running to get to this point before continuing
//     MPI_Barrier(MPI_COMM_WORLD);
//     
//     // Gather all maximum entropy solutions to primary processor 
//     MPI_Gather(&L2_Norm_N3, 1, MPI_LONG_DOUBLE, L2_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
//     MPI_Gather(&L_inf_Norm_N3, 1, MPI_LONG_DOUBLE, L_inf_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
//     
//     // Wait for all the processors running to get to this point before continuing
//     MPI_Barrier(MPI_COMM_WORLD);
//     
//     // Compute L2 and L_inf norm among processors
//     L2_Norm_N3 = 0.0;
//     L_inf_Norm_N3 = 0.0;
//     if (id_proc == PRIMARY_ID) {
//         for (int id_prc = 0; id_prc < num_proc; id_prc++) {
//             L2_Norm_N3 += L2_Norm_N3_array[id_prc];
//             L_inf_Norm_N3 = max(L_inf_Norm_N3, fabs(L_inf_Norm_N3_array[id_prc]));
//         }
//         
//         L2_Norm_N3 /= N_Points_E*N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
//         L2_Norm_N3 = sqrt(L2_Norm_N3);
//         
//         cout << "Interpolation Stats ........." << "L2_Norm_N3 = " << L2_Norm_N3 << "     "  << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << "    " << "max_err_N3_N1_1 = " << max_err_N3_N1_1 << "    " << "max_err_N3_N1_2 = " << max_err_N3_N1_2 << "    " << "max_err_N3_N1_3 = " << max_err_N3_N1_3 << "    " << "max_err_N3_gam1 = " << max_err_N3_gam1 << "    " << "max_err_N3_gam2 = " << max_err_N3_gam2 << endl;
//     }
//     
//     // Wait for all the processors running to get to this point before continuing
//     MPI_Barrier(MPI_COMM_WORLD);
//     
//     if (id_proc == PRIMARY_ID) {
//         delete[] L2_Norm_N3_array; L2_Norm_N3_array = NULL;
//         delete[] L_inf_Norm_N3_array; L_inf_Norm_N3_array = NULL;
//     }
// }

// ******************************************************************************************
// This routine computes the L1 and L2 errors of our interpolative-based approximations of the
// the closing fluxes for the gray or the non-gray M2 closure, at the interpolation nodes
// Note that the error of our proposed approximations at such points should be zero since
// the Vandermonde interpolation is exact at such points
// ******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Compute_L_ONE_L_TWO_Errors_N3_ijk(const int &VAR_NUM) {
    long double max_err_N3_N1_1, max_err_N3_N1_2, max_err_N3_N1_3, max_err_N3_gam1, max_err_N3_gam2;
    long double L2_Norm_N3, L_inf_Norm_N3;
    long double ratio_E, N1_1, N1_2, N1_3, gam1, gam2;
    long double N3_Fit, N3_Numerical;
    int index;
    int id_E_min, id_E_max;
    int id_f_min, id_f_max;
    int id_theta_min, id_theta_max;
    int id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max;
    int id_phi_min, id_phi_max;
    long double *L2_Norm_N3_array = NULL, *L_inf_Norm_N3_array = NULL;
    
    if (id_proc == PRIMARY_ID) {
        L2_Norm_N3_array = new long double[num_proc];
        L_inf_Norm_N3_array = new long double[num_proc];
    }
    
    L2_Norm_N3 = 0.0;
    L_inf_Norm_N3 = 0.0;
    
    Compute_MPI_Processes_Max_Min_Indexes(id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, id_theta_min, id_theta_max, id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max);
    
    for (int i_E = id_E_min; i_E < id_E_max; i_E++) {
        for (int i_f = id_f_min; i_f < id_f_max; i_f++) {
            for (int i_Phi = id_phi_min; i_Phi < id_phi_max; i_Phi++) {
                for (int i_Theta = id_theta_min; i_Theta < id_theta_max; i_Theta++) {
                    for (int i_Triangle_gam1_gam2 = id_Triangle_gam1_gam2_min; i_Triangle_gam1_gam2 < id_Triangle_gam1_gam2_max; i_Triangle_gam1_gam2++) {
                        index = Compute_Full_Id(0, i_E, i_f, i_Phi, i_Theta, i_Triangle_gam1_gam2);
                        
                        // index_temp = Compute_Full_Id_Single_proc(0, i_E, i_f, i_Phi, i_Theta, i_gam1_gam2);
                        
                        ratio_E = E_NON_GRAY[index];
                        N1_1 = N1_1_NON_GRAY[index];
                        N1_2 = N1_2_NON_GRAY[index];
                        N1_3 = N1_3_NON_GRAY[index];
                        gam1 = gam1_NON_GRAY[index];
                        gam2 = gam2_NON_GRAY[index];
                        
                        if (VAR_NUM == N3_111_ENTRY) {
//                             N3_Numerical = f_N3_111_NON_GRAY[index];
//                             N3_Fit = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
                            N3_Numerical = N3_111_NON_GRAY[index];
                            N3_Fit = Evaluate_N3_111(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                        } else if (VAR_NUM == N3_112_ENTRY) {
                            N3_Numerical = N3_112_NON_GRAY[index];
                            N3_Fit = Evaluate_N3_112(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            
//                             cout << "index = " << index << "   " << "ratio_E = " << ratio_E << "   " << "N1_1 = " << N1_1 << "   " << "N1_2 = " << N1_2 << "   " << "N1_3 = " << N1_3 << "   " << "gam1 = " << gam1 << "   " << "gam2 = " << gam2 << "   " << "N3_Fit = " << N3_Fit << "    " << "N3_Numerical = " << N3_Numerical << endl;
                        } else if (VAR_NUM == N3_122_ENTRY) {
//                             N3_Numerical = f_N3_122_NON_GRAY[index];
//                             N3_Fit = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
                            N3_Numerical = N3_122_NON_GRAY[index];
                            N3_Fit = Evaluate_N3_122(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                        } else if (VAR_NUM == N3_123_ENTRY) {
                            N3_Numerical = N3_123_NON_GRAY[index];
                            N3_Fit = Evaluate_N3_123(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                        } else if (VAR_NUM == N3_222_ENTRY) {
                            N3_Numerical = N3_222_NON_GRAY[index];
                            N3_Fit = Evaluate_N3_222(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                        } else {
                            cout << "********************* VAR_NUM not defined !!!!! ************************" << endl;
                            exit(0);
                        }
                        
                        L2_Norm_N3 += pow(N3_Fit - N3_Numerical, 2);
                        L_inf_Norm_N3 = max(L_inf_Norm_N3, fabs(N3_Fit - N3_Numerical));
                        
//                         if (id_proc == PRIMARY_ID) {
//                             cout << "index = " << index << "   " << "ratio_E = " << ratio_E << "   " << "N1_1 = " << N1_1 << "   " << "N1_2 = " << N1_2 << "   " << "N1_3 = " << N1_3 << "   " << "gam1 = " << gam1 << "   " << "gam2 = " << gam2 << "   " << "N3_Fit = " << N3_Fit << "    " << "N3_Numerical = " << N3_Numerical << "    " << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << endl;
//                         }
                        
                        if (L_inf_Norm_N3 == fabs(N3_Fit - N3_Numerical)) {
                            max_err_N3_N1_1 = N1_1;
                            max_err_N3_N1_2 = N1_2;
                            max_err_N3_N1_3 = N1_3;
                            max_err_N3_gam1 = gam1;
                            max_err_N3_gam2 = gam2;
                        }
                    }
                }
            }
        }
    }
    
    // cout << "id_proc = " << id_proc << "  " << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << endl; 
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Gather all maximum entropy solutions to primary processor 
    MPI_Gather(&L2_Norm_N3, 1, MPI_LONG_DOUBLE, L2_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Gather(&L_inf_Norm_N3, 1, MPI_LONG_DOUBLE, L_inf_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Compute L2 and L_inf norm among processors
    L2_Norm_N3 = 0.0;
    L_inf_Norm_N3 = 0.0;
    if (id_proc == PRIMARY_ID) {
        for (int id_prc = 0; id_prc < num_proc; id_prc++) {
            L2_Norm_N3 += L2_Norm_N3_array[id_prc];
            L_inf_Norm_N3 = max(L_inf_Norm_N3, fabs(L_inf_Norm_N3_array[id_prc]));
        }
        
        L2_Norm_N3 /= N_Points_E*N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
        L2_Norm_N3 = sqrt(L2_Norm_N3);
        
        cout << "Interpolation Stats ........." << "L2_Norm_N3 = " << L2_Norm_N3 << "     "  << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << "    " << "max_err_N3_N1_1 = " << max_err_N3_N1_1 << "    " << "max_err_N3_N1_2 = " << max_err_N3_N1_2 << "    " << "max_err_N3_N1_3 = " << max_err_N3_N1_3 << "    " << "max_err_N3_gam1 = " << max_err_N3_gam1 << "    " << "max_err_N3_gam2 = " << max_err_N3_gam2 << endl;
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        delete[] L2_Norm_N3_array; L2_Norm_N3_array = NULL;
        delete[] L_inf_Norm_N3_array; L_inf_Norm_N3_array = NULL;
    }
}


void N3_Non_Gray_M2_3D_RT_Cheby :: Compute_L_ONE_L_TWO_Errors_Derivatives_N3_ijk(const int &VAR_NUM, const int &var_index) {
    long double max_err_N3_N1_1, max_err_N3_N1_2, max_err_N3_N1_3, max_err_N3_gam1, max_err_N3_gam2;
    long double L2_Norm_N3, L_inf_Norm_N3;
    long double ratio_E, N1_1, N1_2, N1_3, gam1, gam2;
    long double N3_Fit, N3_Numerical;
    int index;
    int id_E_min, id_E_max;
    int id_f_min, id_f_max;
    int id_phi_min, id_phi_max;
    int id_theta_min, id_theta_max;
    int id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max;
    long double *L2_Norm_N3_array = NULL, *L_inf_Norm_N3_array = NULL;
    long double I0_star_val = 1.0;
    
    if (id_proc == PRIMARY_ID) {
        L2_Norm_N3_array = new long double[num_proc];
        L_inf_Norm_N3_array = new long double[num_proc];
    }
    
    L2_Norm_N3 = 0.0;
    L_inf_Norm_N3 = 0.0;
    
    Compute_MPI_Processes_Max_Min_Indexes(id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, id_theta_min, id_theta_max, id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max);
    
    for (int i_E = id_E_min; i_E < id_E_max; i_E++) {
        for (int i_f = id_f_min; i_f < id_f_max; i_f++) {
            for (int i_Phi = id_phi_min; i_Phi < id_phi_max; i_Phi++) {
                for (int i_Theta = id_theta_min; i_Theta < id_theta_max; i_Theta++) {
                    for (int i_gam1_gam2 = id_Triangle_gam1_gam2_min; i_gam1_gam2 < id_Triangle_gam1_gam2_max; i_gam1_gam2++) {
                        index = Compute_Full_Id(0, i_E, i_f, i_Phi, i_Theta, i_gam1_gam2);
                        
                        // index_temp = Compute_Full_Id_Single_proc(0, i_E, i_f, i_Phi, i_Theta, i_gam1_gam2);
                        
                        ratio_E = E_NON_GRAY[index];
                        N1_1 = N1_1_NON_GRAY[index];
                        N1_2 = N1_2_NON_GRAY[index];
                        N1_3 = N1_3_NON_GRAY[index];
                        gam1 = gam1_NON_GRAY[index];
                        gam2 = gam2_NON_GRAY[index];
                        
                        Jacobian_M2.Compute_W_array(I0_star_val, N1_1, N1_2, gam1, gam2);
    
                        switch (VAR_NUM) {
                            case N3_111_ENTRY:
                                d_q_RT_dU(Jacobian_M2.dqxxx_RT, Jacobian_M2.dN3_111_RT, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
                                if (var_index == INDEX_I0) {
                                    
                                } else if (var_index == INDEX_N1_1) {
                                    N3_Numerical = dN3_111_NON_GRAY_dN1_1[index];
                                    N3_Fit = Jacobian_M2.dqxxx_RT[1];
                                } else if (var_index == INDEX_N1_2) {
                                    N3_Numerical = dN3_111_NON_GRAY_dN1_2[index];
                                    N3_Fit = Jacobian_M2.dqxxx_RT[2];
                                } else if (var_index == INDEX_N2_11) {
                                    
                                } else if (var_index == INDEX_N2_12) {
                                    
                                } else if (var_index == INDEX_N2_22) {
                                
                                } else {
                                    cout << "Var index not specified !!!!!!!!!!!!!!!1" << endl;
                                }
                                break;
                            case N3_122_ENTRY:
                                d_q_RT_dU(Jacobian_M2.dqxyy_RT, Jacobian_M2.dN3_122_RT, ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
                                if (var_index == INDEX_I0) {
                                    
                                } else if (var_index == INDEX_N1_1) {
                                    N3_Numerical = dN3_122_NON_GRAY_dN1_1[index];
                                    N3_Fit = Jacobian_M2.dqxyy_RT[1];
                                } else if (var_index == INDEX_N1_2) {
                                    N3_Numerical = dN3_122_NON_GRAY_dN1_2[index];
                                    N3_Fit = Jacobian_M2.dqxyy_RT[2];
                                } else if (var_index == INDEX_N2_11) {
                                    
                                } else if (var_index == INDEX_N2_12) {
                                    
                                } else if (var_index == INDEX_N2_22) {
                                
                                } else {
                                    cout << "Var index not specified !!!!!!!!!!!!!!!1" << endl;
                                }
                                break;
                            case N3_123_ENTRY:
                                if (var_index == INDEX_I0) {
                                    
                                } else if (var_index == INDEX_N1_1) {
                                    
                                } else if (var_index == INDEX_N1_2) {
                                    
                                } else if (var_index == INDEX_N2_11) {
                                    
                                } else if (var_index == INDEX_N2_12) {
                                    
                                } else if (var_index == INDEX_N2_22) {
                                
                                } else {
                                    cout << "Var index not specified !!!!!!!!!!!!!!!1" << endl;
                                }
                                N3_Numerical = N3_123_NON_GRAY[index];
                                N3_Fit = Evaluate_N3_123(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                                break;
                            default:
                                cout << "********************* VAR_NUM not defined !!!!! ************************" << endl;
                                exit(0);
                        }
                        
                        L2_Norm_N3 += pow(N3_Fit - N3_Numerical, 2);
                        L_inf_Norm_N3 = max(L_inf_Norm_N3, fabs(N3_Fit - N3_Numerical));
                        
//                         if (id_proc == PRIMARY_ID) {
//                             cout << "index = " << index << "   " << "ratio_E = " << ratio_E << "   " << "N1_1 = " << N1_1 << "   " << "N1_2 = " << N1_2 << "   " << "N1_3 = " << N1_3 << "   " << "gam1 = " << gam1 << "   " << "gam2 = " << gam2 << "   " << "N3_Fit = " << N3_Fit << "    " << "N3_Numerical = " << N3_Numerical << "    " << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << endl;
//                         }
                        
                        if (L_inf_Norm_N3 == fabs(N3_Fit - N3_Numerical)) {
                            max_err_N3_N1_1 = N1_1;
                            max_err_N3_N1_2 = N1_2;
                            max_err_N3_N1_3 = N1_3;
                            max_err_N3_gam1 = gam1;
                            max_err_N3_gam2 = gam2;
                        }
                    }
                }
            }
        }
    }
    
    // cout << "id_proc = " << id_proc << "  " << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << endl; 
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Gather all maximum entropy solutions to primary processor 
    MPI_Gather(&L2_Norm_N3, 1, MPI_LONG_DOUBLE, L2_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Gather(&L_inf_Norm_N3, 1, MPI_LONG_DOUBLE, L_inf_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Compute L2 and L_inf norm among processors
    L2_Norm_N3 = 0.0;
    L_inf_Norm_N3 = 0.0;
    if (id_proc == PRIMARY_ID) {
        for (int id_prc = 0; id_prc < num_proc; id_prc++) {
            L2_Norm_N3 += L2_Norm_N3_array[id_prc];
            L_inf_Norm_N3 = max(L_inf_Norm_N3, fabs(L_inf_Norm_N3_array[id_prc]));
        }
        
        L2_Norm_N3 /= N_Points_E*N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
        L2_Norm_N3 = sqrt(L2_Norm_N3);
        
        cout << "Interpolation Stats ........." << "L2_Norm_N3 = " << L2_Norm_N3 << "     "  << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << "    " << "max_err_N3_N1_1 = " << max_err_N3_N1_1 << "    " << "max_err_N3_N1_2 = " << max_err_N3_N1_2 << "    " << "max_err_N3_N1_3 = " << max_err_N3_N1_3 << "    " << "max_err_N3_gam1 = " << max_err_N3_gam1 << "    " << "max_err_N3_gam2 = " << max_err_N3_gam2 << endl;
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        delete[] L2_Norm_N3_array; L2_Norm_N3_array = NULL;
        delete[] L_inf_Norm_N3_array; L_inf_Norm_N3_array = NULL;
    }
}

// ******************************************************************************************
// This routine computes the L1 and L2 errors of our interpolative-based approximations of the
// the closing fluxes for the gray or the non-gray M2 closure, for a given set of evaluation 
// points uniformly distributed in the realizable space for angular moments up to second-order
// ******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Compute_L_ONE_L_TWO_Errors_N3_ijk(N3_Non_Gray_M2_3D_RT_Cheby &N3_3D_RT_Unif, const int &flag_Write_Output, ofstream &out_L_inf_Norm, const int &VAR_NUM) {
    long double max_err_N3_N1_1, max_err_N3_N1_2, max_err_N3_N1_3, max_err_N3_gam1, max_err_N3_gam2;
    long double L2_Norm_N3, L_inf_Norm_N3;
    long double ratio_E, N1_1, N1_2, N1_3, gam1, gam2;
    long double N3_Fit, N3_Numerical;
    long double Mobius_Scale_Actual, Mobius_Scale_Fit;
    int index;
    int id_E_min, id_E_max;
    int id_f_min, id_f_max;
    int id_phi_min, id_phi_max;
    int id_theta_min, id_theta_max;
    int id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max;
    long double *L2_Norm_N3_array = NULL, *L_inf_Norm_N3_array = NULL;
    
    long double x_L_N3[N3_3D_RT_Unif.N_pts_Mob_Scale], weight_L_N3[N3_3D_RT_Unif.N_pts_Mob_Scale];
    
    if (id_proc == PRIMARY_ID) {
        L2_Norm_N3_array = new long double[num_proc];
        L_inf_Norm_N3_array = new long double[num_proc];
    }
    
    L2_Norm_N3 = 0.0;
    L_inf_Norm_N3 = 0.0;
    
    gen_laguerre_ek_compute ( N3_3D_RT_Unif.N_pts_Mob_Scale, 0.0, x_L_N3, weight_L_N3 );
    
    N3_3D_RT_Unif.Compute_MPI_Processes_Max_Min_Indexes(id_E_min, id_E_max, id_f_min, id_f_max, id_phi_min, id_phi_max, id_theta_min, id_theta_max, id_Triangle_gam1_gam2_min, id_Triangle_gam1_gam2_max);
    
    for (int id_Mobius = 0; id_Mobius < N3_3D_RT_Unif.N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_N3[id_Mobius];
        for (int i_E = id_E_min; i_E < id_E_max; i_E++) {
            for (int i_f = id_f_min; i_f < id_f_max; i_f++) {
                for (int i_Phi = id_phi_min; i_Phi < id_phi_max; i_Phi++) {
                    for (int i_Theta = id_theta_min; i_Theta < id_theta_max; i_Theta++) {
                        for (int i_gam1_gam2 = id_Triangle_gam1_gam2_min; i_gam1_gam2 < id_Triangle_gam1_gam2_max; i_gam1_gam2++) {
                            index = N3_3D_RT_Unif.Compute_Full_Id(id_Mobius, i_E, i_f, i_Phi, i_Theta, i_gam1_gam2);
                                
                            ratio_E = N3_3D_RT_Unif.E_NON_GRAY[index];
                            N1_1 = N3_3D_RT_Unif.N1_1_NON_GRAY[index];
                            N1_2 = N3_3D_RT_Unif.N1_2_NON_GRAY[index];
                            N1_3 = N3_3D_RT_Unif.N1_3_NON_GRAY[index];
                            gam1 = N3_3D_RT_Unif.gam1_NON_GRAY[index];
                            gam2 = N3_3D_RT_Unif.gam2_NON_GRAY[index];
                            
                            Mobius_Scale_Fit = Evaluate_Length_Scale_N3_ijk(N1_1, N1_2, N1_3, gam1, gam2);
                            
                            if (fabs(ratio_E) == 1.0) {
                                ratio_E = ratio_E;
                            } else {
                                //                     cout << "ratio_E Before = " << ratio_E << endl;
                                ratio_E = Inverse_Mobius_Transformation(ratio_E, Mobius_Scale_Actual);
                                //                     cout << "I0 Before = " << ratio_E << endl;
                                ratio_E = Mobius_Transformation(ratio_E, Mobius_Scale_Fit);
                                //                     cout << "ratio_E after = " << ratio_E << endl;
                            }
                            
                            if (VAR_NUM == N3_111_ENTRY) {
                                N3_Numerical = N3_3D_RT_Unif.f_N3_111_NON_GRAY[index];
                                N3_Fit = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_111_ENTRY);
                                // N3_Numerical = N3_3D_RT_Unif.N3_111_NON_GRAY[index];
                                // N3_Fit = Evaluate_N3_111(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            } else if (VAR_NUM == N3_112_ENTRY) {
                                N3_Numerical = N3_3D_RT_Unif.N3_112_NON_GRAY[index];
                                N3_Fit = Evaluate_N3_112(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            } else if (VAR_NUM == N3_122_ENTRY) {
                                N3_Numerical = N3_3D_RT_Unif.f_N3_122_NON_GRAY[index];
                                N3_Fit = Evaluate_g_N3_ijk(ratio_E, N1_1, N1_2, N1_3, gam1, gam2, N3_122_ENTRY);
                                // N3_Numerical = N3_3D_RT_Unif.N3_122_NON_GRAY[index];
                                // N3_Fit = Evaluate_N3_122(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            } else if (VAR_NUM == N3_123_ENTRY) {
                                N3_Numerical = N3_3D_RT_Unif.N3_123_NON_GRAY[index];
                                N3_Fit = Evaluate_N3_123(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            } else if (VAR_NUM == N3_222_ENTRY) {
                                N3_Numerical = N3_3D_RT_Unif.N3_222_NON_GRAY[index];
                                N3_Fit = Evaluate_N3_222(ratio_E, N1_1, N1_2, N1_3, gam1, gam2);
                            } else {
                                cout << "********************* VAR_NUM not defined !!!!! ************************" << endl;
                                exit(0);
                            }
                            
                            if (N3_Numerical !=  N3_Numerical) {
                                cout << "id_proc = " << id_proc << "  " << "index = " << index << "  " << "VAR_NUM = " << VAR_NUM << "  " << "ratio_E = " << ratio_E << "  " << "N1_1 = " << N1_1 << "  " << "N1_2 = " << N1_2 << "  " << "N1_3 = " << N1_3 << "  " << "gam1 = " << gam1 << "  " << "gam2 = " << gam2 << "  " << "N3_Fit = " << N3_Fit << "  " << "N3_Numerical = " << N3_Numerical << endl;
                                
                                cout << "id_proc_E = " << MPI_proc_parameters.id_proc_E << "  " << "id_proc_f = " << MPI_proc_parameters.id_proc_f << "  " << "id_proc_phi = " << MPI_proc_parameters.id_proc_phi << "  " << "id_proc_theta = " << MPI_proc_parameters.id_proc_theta << "  " << "id_proc_Triangle_gam1_gam2 = " << MPI_proc_parameters.id_proc_Triangle_gam1_gam2 << endl;
        
                                cout << "id_E_min = " << id_E_min << "  " << "id_E_max = " << id_E_max << "  " << "id_f_min = " << id_f_min << "  " << "id_f_max = " << id_f_max << "  " << "id_phi_min = " << id_phi_min << "  " << "id_phi_max = " << id_phi_max << "  " << "id_theta_min = " << id_theta_min << "  " << "id_theta_max = " << id_theta_max << "  " << "id_Triangle_gam1_gam2_min = " << id_Triangle_gam1_gam2_min << "  " << "id_Triangle_gam1_gam2_max = " << id_Triangle_gam1_gam2_max << endl;
                                
                                exit(0);
                            }
                            
                            L2_Norm_N3 += pow(N3_Fit - N3_Numerical, 2);
                            L_inf_Norm_N3 = max(L_inf_Norm_N3, fabs(N3_Fit - N3_Numerical));
                                
                            if (L_inf_Norm_N3 == fabs(N3_Fit - N3_Numerical)) {
                                max_err_N3_N1_1 = N1_1;
                                max_err_N3_N1_2 = N1_2;
                                max_err_N3_N1_3 = N1_3;
                                max_err_N3_gam1 = gam1;
                                max_err_N3_gam2 = gam2;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Gather all maximum entropy solutions to primary processor 
    MPI_Gather(&L2_Norm_N3, 1, MPI_LONG_DOUBLE, L2_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    MPI_Gather(&L_inf_Norm_N3, 1, MPI_LONG_DOUBLE, L_inf_Norm_N3_array, 1, MPI_LONG_DOUBLE, PRIMARY_ID, MPI_COMM_WORLD);
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Compute L2 and L_inf norm among processors
    L2_Norm_N3 = 0.0;
    L_inf_Norm_N3 = 0.0;
    if (id_proc == PRIMARY_ID) {
        for (int id_prc = 0; id_prc < N3_3D_RT_Unif.num_proc_used; id_prc++) {
            L2_Norm_N3 += L2_Norm_N3_array[id_prc];
            L_inf_Norm_N3 = max(L_inf_Norm_N3, fabs(L_inf_Norm_N3_array[id_prc]));
        }
        
        L2_Norm_N3 /= N3_3D_RT_Unif.N_pts_Mob_Scale*N3_3D_RT_Unif.N_Points_E*N3_3D_RT_Unif.N_Points_f*N3_3D_RT_Unif.N_Points_Phi*N3_3D_RT_Unif.N_Points_Theta*N3_3D_RT_Unif.N_Points_Triangle_gam1_gam2;
        L2_Norm_N3 = sqrt(L2_Norm_N3);
        
        N3_3D_RT_Unif.L_inf_Norm_N3 = L_inf_Norm_N3;
        
        N3_3D_RT_Unif.L2_Norm_N3 = L2_Norm_N3;
        
        if (flag_Write_Output) {
            out_L_inf_Norm << setw(16) << L_inf_Norm_N3 << setw(16) << L2_Norm_N3 << endl;    
        }
        
        cout << "Convergence Stats ........." << "L2_Norm_N3 = " << L2_Norm_N3 << "     "  << "L_inf_Norm_N3 = " << L_inf_Norm_N3 << "    " << "max_err_N3_N1_1 = " << max_err_N3_N1_1 << "    " << "max_err_N3_N1_2 = " << max_err_N3_N1_2 << "    " << "max_err_N3_N1_3 = " << max_err_N3_N1_3 << "    " << "max_err_N3_gam1 = " << max_err_N3_gam1 << "    " << "max_err_N3_gam2 = " << max_err_N3_gam2 << endl;
    }
    
    // Wait for all the processors running to get to this point before continuing
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        delete[] L2_Norm_N3_array; L2_Norm_N3_array = NULL;
        delete[] L_inf_Norm_N3_array; L_inf_Norm_N3_array = NULL;
    }
}

// ******************************************************************************************
// This routine computes maximum-entropy solutions for sets of angular moments up to second-
// order associated with the interpolation nodes for the purpose of our interpolative-based
// procedure, with values of the Length scale based on the solution of the non-linear
// least-squares solver
// ******************************************************************************************
void N3_Non_Gray_M2_3D_RT_Cheby :: Precompute_Final_Max_Ent_Solution(const int &Maximum_Entropy_Solution_Regime) {
    int Node_Distribution_E, Node_Distribution_f, Node_Distribution_gam1, Node_Distribution_gam2, Node_Distribution_gam3;
    MN_Var_Num_Points Num_points;
    int NVARS = N_Points_f*N_Points_Triangle_gam1_gam2;
    int Optim_Flag;
    int Npts_Temp;
    int size_rec = MPI_proc_parameters.size_rec;
    
    record_N3 *rec_N3_local, *rec_N3_global;
    
    M2_State_Param M2_State(9);
    
    rec_N3_local = new record_N3[size_rec];
    
    Npts_Temp = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
    rec_N3_global = new record_N3[Npts_Temp];
    
    fstream in1_out;
    
    Node_Distribution_E = CHEBYSHEV_SECOND_KIND_DISTRIBUTION; //CHEBYSHEV_FIRST_KIND_DISTRIBUTION; //
    Node_Distribution_f = CHEBYSHEV_FIRST_KIND_DISTRIBUTION; //CHEBYSHEV_SECOND_KIND_DISTRIBUTION;
    
    // cout << "Change back distribution !!!!" << endl;
    Node_Distribution_Theta_Phi = SPHERICAL_HARMONIC_DISTRIBUTION; //UNIFORM_DISTRIBUTION; //
    Node_Distribution_gam1 = CHEBYSHEV_FIRST_KIND_DISTRIBUTION; //CHEBYSHEV_SECOND_KIND_DISTRIBUTION;
    
    Num_points.E = N_Points_E;
    Num_points.f = N_Points_f;
    Num_points.phi = N_Points_Phi;
    Num_points.theta = N_Points_Theta;
    Num_points.gam1 = N_Points_gam1;
    Num_points.N_Points_Triangle_gam1_gam2 = N_Points_Triangle_gam1_gam2;
    Num_points.Order_SH = 0; //Order_SH; //
    Num_points.N_pts_SH = 1; //N_Coeffs_SH; //
    Num_points.Maximum_Entropy_Solution_Regime = Maximum_Entropy_Solution_Regime;
    
    Num_points.Length_Scale_Dist_Type = LENGTH_SCALE_DIST_FIT;
    
    M2_State.Dimension = 3;
    M2_State.Problem_Type = Problem_Type;
    M2_State.id_proc = id_proc;
    M2_State.num_proc = num_proc;
    M2_State.Node_Dist_E = Node_Distribution_E;
    M2_State.Node_Dist_f = Node_Distribution_f;
    M2_State.Node_Dist_Phi_Theta = Node_Distribution_Theta_Phi;
    M2_State.Node_Dist_gam1 = Node_Distribution_gam1;
    M2_State.Triangle_Domain = FULL_TRIANGLE;
    M2_State.display = false; //true; //
    M2_State.proc_display = PRIMARY_ID;
    
    M2_State.Rec_Moms_Test.Copy(Rec_Moms_Test);
    
    if (id_proc == PRIMARY_ID) {
        cout << "Number of Points ..........................." << "N_Points_E = " << Num_points.E << "     "  << "N_Points_f = " << Num_points.f << "     "  << "N_Points_Phi = " << Num_points.phi << "     "  << "N_Points_Theta = " << Num_points.theta << "     "  << "N_Points_gam1 = " << Num_points.gam1 << endl;
    }
    
    for (int i = 0; i < NVARS; i++) {
        // cout << "Coefficient_Matrix_Fit_Mob_Scale = " << Coefficient_Mobius_Scale_N3_ijk[i] << endl;
    }
    
    Optim_Flag = Optimization_Algorithm_Array(rec_N3_local, M2_State, &Num_points, &MPI_proc_parameters, 1, Coefficient_Mobius_Scale_N3_ijk);
    
    if (id_proc == PRIMARY_ID) {
        cout << "Done solving maximum entropy problem !!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Gather all maximum entropy solutions to primary processor 
    MPI_Gather(rec_N3_local, size_rec, rec_N3_type, rec_N3_global, size_rec, rec_N3_type, PRIMARY_ID, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        cout << "Done gathering maximum entropy solution !!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    
    if (id_proc == PRIMARY_ID) {
        Reorder_Maximum_Entropy_Data_MPI(rec_N3_global);
    }
    
    if (id_proc == PRIMARY_ID) {
        cout << "Done reordering maximum entropy solution !!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    
    int Npts_Total = Compute_Num_Pts_Full();
    // Broadcast rec_N3_global to other processors
    MPI_Bcast(rec_N3_global, Npts_Total, rec_N3_type, PRIMARY_ID, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id_proc == PRIMARY_ID) {
        cout << "Done broadcasting maximum entropy solution !!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    
    if (!Rec_Moms_Test.flag_Moments_Test && iteration == LAST_ITERATION) {
        if (id_proc == PRIMARY_ID) {
            Write_Maximum_Entropy_Data(rec_N3_global);
        }
        if (id_proc == PRIMARY_ID) {
            cout << "Writing maximum entropy solution completed !!!!!!!!!!!!!!!!!!!!!" << endl;
        }
    }
    
    Read_Maximum_Entropy_Data(rec_N3_global);
    
    if (id_proc == PRIMARY_ID) {
        cout << "Done rereading maximum entropy solution !!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    
    delete[] rec_N3_local; rec_N3_local = NULL;
    delete[] rec_N3_global; rec_N3_global = NULL;
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Reorder_Maximum_Entropy_Data_MPI(record_N3 *rec_N3_global) {
    int num_proc_per_var_E, num_proc_per_var_f, num_proc_per_var_phi, num_proc_per_var_theta, num_proc_per_var_Triangle_gam1_gam2;
    int N_Points_E_per_proc, N_Points_f_per_proc, N_Points_phi_per_proc, N_Points_theta_per_proc, N_Points_Triangle_gam1_gam2_Per_Proc;
    int index, index_MPI;
    int size_rec = MPI_proc_parameters.size_rec;
    int Npts_Temp;
    record_N3 *rec_N3_global_temp;
    Npts_Temp = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
    rec_N3_global_temp = new record_N3[Npts_Temp];
    
    for (int i = 0; i < Npts_Temp; i++) {
        rec_N3_global_temp[i] = rec_N3_global[i];
    }
    
    num_proc_per_var_E = MPI_proc_parameters.num_proc_per_var_E;
    num_proc_per_var_f = MPI_proc_parameters.num_proc_per_var_f;
    num_proc_per_var_phi = MPI_proc_parameters.num_proc_per_var_phi;
    num_proc_per_var_theta = MPI_proc_parameters.num_proc_per_var_theta;
    num_proc_per_var_Triangle_gam1_gam2 = MPI_proc_parameters.num_proc_per_var_Triangle_gam1_gam2;
    
    N_Points_E_per_proc = MPI_proc_parameters.N_Points_E_per_proc;
    N_Points_f_per_proc = MPI_proc_parameters.N_Points_f_per_proc;
    N_Points_phi_per_proc = MPI_proc_parameters.N_Points_Phi_per_proc;
    N_Points_theta_per_proc = MPI_proc_parameters.N_Points_Theta_per_proc;
    N_Points_Triangle_gam1_gam2_Per_Proc = MPI_proc_parameters.N_Points_Triangle_gam1_gam2_Per_Proc;
    
    int increment = 0;
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        for (int i_proc_E = 0; i_proc_E < num_proc_per_var_E; i_proc_E++) {
            for (int id_E = 0; id_E < N_Points_E_per_proc; id_E++) {
                for (int i_proc_f = 0; i_proc_f < num_proc_per_var_f; i_proc_f++) {
                    for (int id_f = 0; id_f < N_Points_f_per_proc; id_f++) {
                        for (int i_proc_phi = 0; i_proc_phi < num_proc_per_var_phi; i_proc_phi++) {
                            for (int id_phi = 0; id_phi < N_Points_phi_per_proc; id_phi++) {
                                for (int i_proc_theta = 0; i_proc_theta < num_proc_per_var_theta; i_proc_theta++) {
                                    for (int id_theta = 0; id_theta < N_Points_theta_per_proc; id_theta++) {
                                        for (int i_proc_Triangle_gam1_gam2 = 0; i_proc_Triangle_gam1_gam2 < num_proc_per_var_Triangle_gam1_gam2; i_proc_Triangle_gam1_gam2++) {
                                        for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2_Per_Proc; i_gam1_gam2++) {
                                            index_MPI = i_proc_E*num_proc_per_var_f + i_proc_f;
                                            index_MPI = index_MPI*num_proc_per_var_phi + i_proc_phi;
                                            index_MPI = index_MPI*num_proc_per_var_theta + i_proc_theta;
                                            index_MPI = index_MPI*num_proc_per_var_Triangle_gam1_gam2 + i_proc_Triangle_gam1_gam2;
                                            index_MPI *= size_rec;
                                            
                                            index = id_E*N_Points_f_per_proc + id_f;
                                            index = index*N_Points_phi_per_proc + id_phi;
                                            index = index*N_Points_theta_per_proc + id_theta;
                                            index = index*N_Points_Triangle_gam1_gam2_Per_Proc + i_gam1_gam2;
                                            
                                            index_MPI += index;
                                            
                                            rec_N3_global[increment] = rec_N3_global_temp[index_MPI];
                                            increment++;
                                            
//                                             cout << "id_Mobius = " << id_Mobius << "   " << "ratio_I0 = " << rec_N3_global[increment].ratio_I0 << "   " << "I0 = " << rec_N3_global[increment].I0 << "   " << "N1_1 = " << rec_N3_global[increment].N1_1 << "   " << "N1_2 = " << rec_N3_global[increment].N1_2 << "   " << "N1_3 = " << rec_N3_global[increment].N1_3 << "   " << "gam1 = " << rec_N3_global[increment].gam1 << "   " << "gam2 = " << rec_N3_global[increment].gam2 << endl;
                                            
                                        } // end for i_gam1_gam2
                                        } // end for i_proc_Triangle_gam1_gam2
                                    } // end for id_theta
                                } // end for i_proc_theta
                            } // end for id_phi
                        } // end for i_proc_phi
                    } // end for id_f
                } // end for i_proc_f
            } // end for id_E
        } // end for i_proc_E
    } // end for id_Mobius
    
    if (increment != Npts_Temp) {
        cout << "increment = " << increment << "  " << "Npts_Temp = " << Npts_Temp << endl;
        exit(0);
    }
    
    delete[] rec_N3_global_temp;
    rec_N3_global_temp = NULL;
}

void N3_Non_Gray_M2_3D_RT_Cheby :: Reorder_Maximum_Entropy_Data_f_N3_ijk_MPI(long double *f_N3_111_NON_GRAY, long double *f_N3_122_NON_GRAY, long double *f_N3_123_NON_GRAY) {
    int num_proc_per_var_E, num_proc_per_var_f, num_proc_per_var_phi, num_proc_per_var_theta, num_proc_per_var_Triangle_gam1_gam2;
    int N_Points_E_per_proc, N_Points_f_per_proc, N_Points_phi_per_proc, N_Points_theta_per_proc, N_Points_Triangle_gam1_gam2_Per_Proc;
    int index, index_MPI;
    int size_rec = MPI_proc_parameters.size_rec;
    int Npts_Temp;
    long double *f_N3_111_NON_GRAY_temp, *f_N3_122_NON_GRAY_temp, *f_N3_123_NON_GRAY_temp;
    Npts_Temp = N_pts_Mob_Scale*N_Points_E*N_Points_f*N_Points_Phi*N_Points_Theta*N_Points_Triangle_gam1_gam2;
    f_N3_111_NON_GRAY_temp = new long double[Npts_Temp];
    f_N3_122_NON_GRAY_temp = new long double[Npts_Temp];
    f_N3_123_NON_GRAY_temp = new long double[Npts_Temp];
    
    for (int i = 0; i < Npts_Temp; i++) {
        f_N3_111_NON_GRAY_temp[i] = f_N3_111_NON_GRAY[i];
        f_N3_122_NON_GRAY_temp[i] = f_N3_122_NON_GRAY[i];
        f_N3_123_NON_GRAY_temp[i] = f_N3_123_NON_GRAY[i];
    }
    
    num_proc_per_var_E = MPI_proc_parameters.num_proc_per_var_E;
    num_proc_per_var_f = MPI_proc_parameters.num_proc_per_var_f;
    num_proc_per_var_phi = MPI_proc_parameters.num_proc_per_var_phi;
    num_proc_per_var_theta = MPI_proc_parameters.num_proc_per_var_theta;
    num_proc_per_var_Triangle_gam1_gam2 = MPI_proc_parameters.num_proc_per_var_Triangle_gam1_gam2;
    
    N_Points_E_per_proc = MPI_proc_parameters.N_Points_E_per_proc;
    N_Points_f_per_proc = MPI_proc_parameters.N_Points_f_per_proc;
    N_Points_phi_per_proc = MPI_proc_parameters.N_Points_Phi_per_proc;
    N_Points_theta_per_proc = MPI_proc_parameters.N_Points_Theta_per_proc;
    N_Points_Triangle_gam1_gam2_Per_Proc = MPI_proc_parameters.N_Points_Triangle_gam1_gam2_Per_Proc;
    
    int increment = 0;
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        for (int i_proc_E = 0; i_proc_E < num_proc_per_var_E; i_proc_E++) {
            for (int id_E = 0; id_E < N_Points_E_per_proc; id_E++) {
                for (int i_proc_f = 0; i_proc_f < num_proc_per_var_f; i_proc_f++) {
                    for (int id_f = 0; id_f < N_Points_f_per_proc; id_f++) {
                        for (int i_proc_phi = 0; i_proc_phi < num_proc_per_var_phi; i_proc_phi++) {
                            for (int id_phi = 0; id_phi < N_Points_phi_per_proc; id_phi++) {
                                for (int i_proc_theta = 0; i_proc_theta < num_proc_per_var_theta; i_proc_theta++) {
                                    for (int id_theta = 0; id_theta < N_Points_theta_per_proc; id_theta++) {
                                        for (int i_proc_Triangle_gam1_gam2 = 0; i_proc_Triangle_gam1_gam2 < num_proc_per_var_Triangle_gam1_gam2; i_proc_Triangle_gam1_gam2++) {
                                        for (int i_gam1_gam2 = 0; i_gam1_gam2 < N_Points_Triangle_gam1_gam2_Per_Proc; i_gam1_gam2++) {
                                            index_MPI = i_proc_E*num_proc_per_var_f + i_proc_f;
                                            index_MPI = index_MPI*num_proc_per_var_phi + i_proc_phi;
                                            index_MPI = index_MPI*num_proc_per_var_theta + i_proc_theta;
                                            index_MPI = index_MPI*num_proc_per_var_Triangle_gam1_gam2 + i_proc_Triangle_gam1_gam2;
                                            index_MPI *= size_rec;
                                            
                                            index = id_E*N_Points_f_per_proc + id_f;
                                            index = index*N_Points_phi_per_proc + id_phi;
                                            index = index*N_Points_theta_per_proc + id_theta;
                                            index = index*N_Points_Triangle_gam1_gam2_Per_Proc + i_gam1_gam2;
                                            
                                            index_MPI += index;
                                            
                                            f_N3_111_NON_GRAY[increment] = f_N3_111_NON_GRAY_temp[index_MPI];
                                            f_N3_122_NON_GRAY[increment] = f_N3_122_NON_GRAY_temp[index_MPI];
                                            f_N3_123_NON_GRAY[increment] = f_N3_123_NON_GRAY_temp[index_MPI];
                                            increment++;
                                        } // end for i_gam1_gam2
                                        } // end for i_proc_Triangle_gam1_gam2
                                    } // end for id_theta
                                } // end for i_proc_theta
                            } // end for id_phi
                        } // end for i_proc_phi
                    } // end for id_f
                } // end for i_proc_f
            } // end for id_E
        } // end for i_proc_E
    } // end for id_Mobius
    
    if (increment != Npts_Temp) {
        cout << "increment = " << increment << "  " << "Npts_Temp = " << Npts_Temp << endl;
        exit(0);
    }
    
    delete[] f_N3_111_NON_GRAY_temp;
    f_N3_111_NON_GRAY_temp = NULL;
    
    delete[] f_N3_122_NON_GRAY_temp;
    f_N3_122_NON_GRAY_temp = NULL;
    
    delete[] f_N3_123_NON_GRAY_temp;
    f_N3_123_NON_GRAY_temp = NULL;
}

// //********************************************************************************************************
// // This routine sets up values of the weighting function, g_N3_ijk, of the polynomial interpolant, 
// // f_N3_ijk, for the purpose of our interpolative-based approximations of the closing fluxes for either
// // the gray or the non-gray M2 closure
// //********************************************************************************************************
// void N3_Non_Gray_M2_3D_RT_Cheby :: SetupInterpolant_Values_HL_LL() {
//     int index;
//     long double norm_f, norm_f_2, gam3;
//     
//     if (N_Points_E != 1) {
//         cout << "Problem with number of points in HL or LL, N_Points_E = " << N_Points_E << endl;
//         exit(0);
//     }
//     
//     for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
//         for (int i_Phi = 0 ; i_Phi < N_Points_Phi; i_Phi++) {
//             for (int i_Theta = 0 ; i_Theta < N_Points_Theta; i_Theta++) {
//                 for (int i_Cheby_gam1_gam2 = 0; i_Cheby_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_Cheby_gam1_gam2++) {
//                     index = (i_Cheby_f*N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
//                     index = index*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2;
//                     
//                     norm_f_2 = pow(N1_1_NON_GRAY[index],2) + pow(N1_2_NON_GRAY[index],2) + pow(N1_3_NON_GRAY[index],2);
//                     norm_f = sqrt(norm_f_2);
//                     
//                     // Compute f_N3_111
//                     if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(fabs(1.0 - gam1_NON_GRAY[index]) < 1.0e-8)) {
//                         // gam1 > 0 and gam1 < 1
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_111_NON_GRAY[index] = -d2_N3_111_NON_GRAY_dnorm_f_dN1_1[index]/2.0;
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dN1_1[index]/(1.0 - norm_f_2);
//                             }
//                         } else {
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_111_NON_GRAY[index] = (3.0/2.0)*pow(N1_1_NON_GRAY[index], 2) - dN3_111_NON_GRAY_dnorm_f[index]/(2.0*N1_1_NON_GRAY[index]);
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_111_NON_GRAY[index] = (N3_111_NON_GRAY[index] - pow(N1_1_NON_GRAY[index], 3));
//                                 f_N3_111_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                             }
//                         }
//                         f_N3_111_NON_GRAY[index] = f_N3_111_NON_GRAY[index] - gam1_NON_GRAY[index];
//                         f_N3_111_NON_GRAY[index] /= gam1_NON_GRAY[index]*(1.0 - gam1_NON_GRAY[index]);
//                     } else if (gam1_NON_GRAY[index] < 1.0e-8) {
//                         // gam1 = 0
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_111_NON_GRAY[index] = -d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_111_NON_GRAY[index] = d2_N3_111_NON_GRAY_dgam1_dN1_1[index]/(1.0 - norm_f_2);
//                             }
//                         } else {
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_111_NON_GRAY[index] = -d2_N3_111_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dgam1[index];
//                                 f_N3_111_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                             }
//                         }
//                         f_N3_111_NON_GRAY[index] = (f_N3_111_NON_GRAY[index] - 1.0)/(1.0 - gam1_NON_GRAY[index]);
//                     } else if (fabs(1.0 - gam1_NON_GRAY[index]) < 1.0e-8) {
//                         // gam1 = 1
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_111_NON_GRAY[index] = -d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;  
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_111_NON_GRAY[index] = d2_N3_111_NON_GRAY_dgam1_dN1_1[index]/(1.0 - norm_f_2);
//                             }
//                         } else {
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_111_NON_GRAY[index] = -d2_N3_111_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dgam1[index];
//                                 f_N3_111_NON_GRAY[index] = N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                             }
//                         }
//                         f_N3_111_NON_GRAY[index] = (1.0 - f_N3_111_NON_GRAY[index])/gam1_NON_GRAY[index];
//                     } else {
//                         cout << "Error in Setup Interpolant Values for N3 111" << endl;
//                         exit(0);
//                     }
//                     
//                     // Compute f_N3_122
//                     if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8)) {
//                         // gam1 > 0 and gam1 < 1 and gam2 > 0 and gam2 < 1 - gam1
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_122_NON_GRAY[index] = pow(N1_2_NON_GRAY[index], 2) - d2_N3_122_NON_GRAY_dnorm_f_dN1_1[index]/2.0;
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dN1_1[index] - pow(N1_2_NON_GRAY[index], 2);
//                                 f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
//                             }
//                         } else {
//                             // N1_1 != 0
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_122_NON_GRAY[index] = (3.0/2.0)*pow(N1_2_NON_GRAY[index], 2) - dN3_122_NON_GRAY_dnorm_f[index]/(2.0*N1_1_NON_GRAY[index]);
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_122_NON_GRAY[index] = (N3_122_NON_GRAY[index] - N1_1_NON_GRAY[index]*pow(N1_2_NON_GRAY[index], 2));
//                                 f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                             }
//                         }
//                         f_N3_122_NON_GRAY[index] = f_N3_122_NON_GRAY[index] - gam2_NON_GRAY[index];
//                         f_N3_122_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam2_NON_GRAY[index];
//                     } else if ((gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8)) {
//                         // gam1 = 0 and gam2 > 0 and gam2 < 1 - gam1
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam1_dN1_1[index];
//                                 f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
//                             }
//                         } else {
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_122_NON_GRAY[index] = -d2_N3_122_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dgam1[index];
//                                 f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                             }
//                         }
//                         f_N3_122_NON_GRAY[index] /= gam2_NON_GRAY[index];
//                     } else if (!(gam1_NON_GRAY[index] < 1.0e-8) && (gam2_NON_GRAY[index] < 1.0e-8)) {
//                         // gam1 > 0 and gam1 < 1 and gam2 = 0
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1[index]/2.0;
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam2_dN1_1[index];
//                                 f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
//                             }
//                         } else {
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_122_NON_GRAY[index] = -d2_N3_122_NON_GRAY_dnorm_f_dgam2[index]/(2.0*N1_1_NON_GRAY[index]);
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dgam2[index];
//                                 f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                             }
//                         }
//                         f_N3_122_NON_GRAY[index] = (f_N3_122_NON_GRAY[index] - 1.0)/gam1_NON_GRAY[index];
//                     } else if ((gam1_NON_GRAY[index] < 1.0e-8) && (gam2_NON_GRAY[index] < 1.0e-8)) {
//                         // gam1 = 0 && gam2 = 0
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_122_NON_GRAY[index] = -d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1[index]/2.0;
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_122_NON_GRAY[index] = d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1[index];
//                                 f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
//                             }
//                         } else {
//                             if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                 // norm_f = 1
//                                 f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2[index]/(2.0*N1_1_NON_GRAY[index]);
//                             } else {
//                                 // norm_f < 1
//                                 f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam1_dgam2[index];
//                                 f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                             }
//                         }
//                     } else {
//                         cout << "Error in Setup Interpolant Values for N3 122" << endl;
//                         exit(0);
//                     }
//                     
//                     // Compute f_N3_123
//                     gam3 = 1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index];
//                     if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8) && !(fabs(1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index]) < 1.0e-8)) {
//                         // gam1 > 0 and gam2 > 0 and gam3 > 0
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_2 = 0
//                                 if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_3 = 0
//                                     f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3[index] - 1.0;
//                                 } else {
//                                     // N1_3 >= 0
//                                     f_N3_123_NON_GRAY[index] = (d2_N3_123_NON_GRAY_dN1_1_dN1_2[index]/N1_3_NON_GRAY[index]) - 1.0;  
//                                 }
//                             } else {
//                                 // N1_2 >= 0
//                                 if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_3 = 0
//                                     f_N3_123_NON_GRAY[index] = (d2_N3_123_NON_GRAY_dN1_1_dN1_3[index]/N1_2_NON_GRAY[index]) - 1.0;
//                                 } else {
//                                     // N1_3 >= 0
//                                     f_N3_123_NON_GRAY[index] = (dN3_123_NON_GRAY_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index])) - 1.0;  
//                                 }
//                             }
//                         } else {
//                             // N1_1 >= 0
//                             if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_2 = 0
//                                 if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_3 = 0
//                                     f_N3_123_NON_GRAY[index] = (d2_N3_123_NON_GRAY_dN1_2_dN1_3[index]/N1_1_NON_GRAY[index]) - 1.0;
//                                 } else {
//                                     // N1_3 >= 0
//                                     f_N3_123_NON_GRAY[index] = (dN3_123_NON_GRAY_dN1_2[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index])) - 1.0;  
//                                 }
//                             } else {
//                                 // N1_2 >= 0
//                                 if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_3 = 0
//                                     f_N3_123_NON_GRAY[index] = (dN3_123_NON_GRAY_dN1_3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index])) - 1.0;
//                                 } else {
//                                     // N1_3 >= 0
//                                     f_N3_123_NON_GRAY[index] = (N3_123_NON_GRAY[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index])) - 1.0;  
//                                 }
//                             }
//                         }
//                         f_N3_123_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam2_NON_GRAY[index]*gam3;
//                     } else if (gam1_NON_GRAY[index] < 1.0e-8) {
//                         // gam1 = 0
//                         if (gam2_NON_GRAY[index] < 1.0e-8) {
//                             // gam2 = 0
//                             if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_1 = 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2[index]/N1_3_NON_GRAY[index];  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2[index]/N1_2_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             } else {
//                                 // N1_1 >= 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2[index]/N1_1_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam1_dgam2[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             }
//                             f_N3_123_NON_GRAY[index] /= gam3;
//                         } else if (fabs(1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index]) < 1.0e-8) {
//                             // gam3 = 0
//                             if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_1 = 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3[index]/N1_3_NON_GRAY[index];  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3[index]/N1_2_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             } else {
//                                 // N1_1 >= 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3[index]/N1_1_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam1_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             }
//                             f_N3_123_NON_GRAY[index] /= gam2_NON_GRAY[index];
//                         } else {
//                             // gam1 = 0 && gam2 > 0 && gam3 > 0
//                             if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_1 = 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1[index]/N1_3_NON_GRAY[index];  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1[index]/N1_2_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam1_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]); 
//                                     }
//                                 }
//                             } else {
//                                 // N1_1 >= 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1[index]/N1_1_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_2_dgam1[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_3_dgam1[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dgam1[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             }
//                             f_N3_123_NON_GRAY[index] /= gam2_NON_GRAY[index]*gam3;
//                         }
//                     } else if (gam2_NON_GRAY[index] < 1.0e-8) {
//                         // gam1 > 0 and gam2 = 0
//                         if (fabs(1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index]) < 1.0e-8) {
//                             // gam3 = 0
//                             if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_1 = 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3[index]/N1_3_NON_GRAY[index];  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3[index]/N1_2_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             } else {
//                                 // N1_1 >= 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3[index]/N1_1_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam2_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             }
//                             f_N3_123_NON_GRAY[index] /= gam1_NON_GRAY[index];
//                         } else {
//                             // gam1 > 0 && gam2 = 0 && gam3 > 0
//                             if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_1 = 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2[index]/N1_3_NON_GRAY[index];  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2[index]/N1_2_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam2_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             } else {
//                                 // N1_1 >= 0
//                                 if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_2 = 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2[index]/N1_1_NON_GRAY[index];
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_2_dgam2[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 } else {
//                                     // N1_2 >= 0
//                                     if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_3 = 0
//                                         f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_3_dgam2[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                     } else {
//                                         // N1_3 >= 0
//                                         f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dgam2[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                     }
//                                 }
//                             }
//                             f_N3_123_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam3;
//                         }
//                     } else if (fabs(1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index]) < 1.0e-8) {
//                         // gam1 > 0 and gam2 > 0 and gam3 = 0
//                         if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                             // N1_1 = 0
//                             if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_2 = 0
//                                 if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_3 = 0
//                                     f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3[index];
//                                 } else {
//                                     // N1_3 >= 0
//                                     f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3[index]/N1_3_NON_GRAY[index];  
//                                 }
//                             } else {
//                                 // N1_2 >= 0
//                                 if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_3 = 0
//                                     f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3[index]/N1_2_NON_GRAY[index];
//                                 } else {
//                                     // N1_3 >= 0
//                                     f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_1_dgam3[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                 }
//                             }
//                         } else {
//                             // N1_1 >= 0
//                             if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                 // N1_2 = 0
//                                 if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_3 = 0
//                                     f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3[index]/N1_1_NON_GRAY[index];
//                                 } else {
//                                     // N1_3 >= 0
//                                     f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_2_dgam3[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                 }
//                             } else {
//                                 // N1_2 >= 0
//                                 if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_3 = 0
//                                     f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_3_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                 } else {
//                                     // N1_3 >= 0
//                                     f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                 }
//                             }
//                         }
//                         f_N3_123_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam2_NON_GRAY[index];
//                     } else {
//                         cout << "Error in Setup Interpolant Values for N3 123" << endl;
//                         exit(0);
//                     }
//                     
//                     // cout /*<< "index = " << index << "   "*/ << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "   " << "N3_111 = " << N3_111_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY[index] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index] << "  " << "dN3_122 = " << d2_N3_122_NON_GRAY_dgam2_dN1_1[index] << endl;
//                     
//                     //                                                     if (fabs(f_N3_111_NON_GRAY[index]) > 1.0 || fabs(f_N3_122_NON_GRAY[index]) > 1.0) {
// //                     if (f_N3_111_NON_GRAY[index] != f_N3_111_NON_GRAY[index]) {
// //                         cout << "N1_1 = " << N1_1_NON_GRAY[index] << "  " << "N1_2 = " << N1_2_NON_GRAY[index] << "  " << "N1_3 = " << N1_3_NON_GRAY[index] << "  " << "gam1 = " << gam1_NON_GRAY[index] << "  " << "gam2 = " << gam2_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY[index] /*<< "  " << "N3_122 = " << N3_122_NON_GRAY[index] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index]*/ << endl;
// //                         exit(0);
// //                     }
//                 } // end for i_Cheby_gam1_gam2
//             } // end for i_Theta
//         } // end for i_Phi
//     } // end for i_Cheby_f
// }

// void N3_Non_Gray_M2_3D_RT_Cheby :: SetupInterpolant_Values_BE(N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_HL, N3_Non_Gray_M2_3D_RT_Cheby &M2_3D_Data_N3_LL) {
//     int index, index_HL_LL;
//     long double norm_f, norm_f_2, gam3;
//     long double g_N3_111_HL, g_N3_111_LL;
//     long double g_N3_122_HL, g_N3_122_LL;
//     long double g_N3_123_HL, g_N3_123_LL;
//     
//     for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
//         for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
//             for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
//                 for (int i_Phi = 0 ; i_Phi < N_Points_Phi; i_Phi++) {
//                     for (int i_Theta = 0 ; i_Theta < N_Points_Theta; i_Theta++) {
//                         for (int i_Cheby_gam1_gam2 = 0; i_Cheby_gam1_gam2 < N_Points_Triangle_gam1_gam2; i_Cheby_gam1_gam2++) {
//                             index = id_Mobius*N_Points_E + i_Cheby_E;
//                             index = index*N_Points_f + i_Cheby_f;
//                             index = (index*N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
//                             index = index*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2;
//                             
//                             index_HL_LL = (i_Cheby_f*N_Points_Phi + i_Phi)*N_Points_Theta + i_Theta;
//                             index_HL_LL = index_HL_LL*N_Points_Triangle_gam1_gam2 + i_Cheby_gam1_gam2;
//                             
//                             norm_f_2 = pow(N1_1_NON_GRAY[index],2) + pow(N1_2_NON_GRAY[index],2) + pow(N1_3_NON_GRAY[index],2);
//                             norm_f = sqrt(norm_f_2);
//                             
//                             g_N3_111_HL = M2_3D_Data_N3_HL.f_N3_111_NON_GRAY[index_HL_LL];
//                             g_N3_111_LL = M2_3D_Data_N3_LL.f_N3_111_NON_GRAY[index_HL_LL];
//                             
//                             g_N3_122_HL = M2_3D_Data_N3_HL.f_N3_122_NON_GRAY[index_HL_LL];
//                             g_N3_122_LL = M2_3D_Data_N3_LL.f_N3_122_NON_GRAY[index_HL_LL];
//                             
//                             g_N3_123_HL = M2_3D_Data_N3_HL.f_N3_123_NON_GRAY[index_HL_LL];
//                             g_N3_123_LL = M2_3D_Data_N3_LL.f_N3_123_NON_GRAY[index_HL_LL];
//                             
//                             // Compute f_N3_111
//                             if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(fabs(1.0 - gam1_NON_GRAY[index]) < 1.0e-8)) {
//                                 // gam1 > 0 and gam1 < 1
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY[index] = -d2_N3_111_NON_GRAY_dnorm_f_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dN1_1[index]/(1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY[index] = (3.0/2.0)*pow(N1_1_NON_GRAY[index], 2) - dN3_111_NON_GRAY_dnorm_f[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY[index] = (N3_111_NON_GRAY[index] - pow(N1_1_NON_GRAY[index], 3));
//                                         f_N3_111_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_111_NON_GRAY[index] = f_N3_111_NON_GRAY[index] - gam1_NON_GRAY[index];
//                                 f_N3_111_NON_GRAY[index] /= gam1_NON_GRAY[index]*(1.0 - gam1_NON_GRAY[index]);
//                             } else if ((gam1_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 = 0
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY[index] = -d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY[index] = d2_N3_111_NON_GRAY_dgam1_dN1_1[index]/(1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY[index] = -d2_N3_111_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dgam1[index];
//                                         f_N3_111_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_111_NON_GRAY[index] = (f_N3_111_NON_GRAY[index] - 1.0)/(1.0 - gam1_NON_GRAY[index]);
//                             } else if (fabs(1.0 - gam1_NON_GRAY[index]) < 1.0e-8) {
//                                 // gam1 = 1
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY[index] = -d3_N3_111_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;  
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY[index] = d2_N3_111_NON_GRAY_dgam1_dN1_1[index]/(1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_111_NON_GRAY[index] = -d2_N3_111_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_111_NON_GRAY[index] = dN3_111_NON_GRAY_dgam1[index];
//                                         f_N3_111_NON_GRAY[index] = N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_111_NON_GRAY[index] = (1.0 - f_N3_111_NON_GRAY[index])/gam1_NON_GRAY[index];
//                             } else {
//                                 cout << "Error in Setup Interpolant Values for N3 111" << endl;
//                                 exit(0);
//                             }
//                             
//                             // Compute f_N3_122
//                             if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 > 0 and gam1 < 1 and gam2 > 0 and gam2 < 1 - gam1
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY[index] = pow(N1_2_NON_GRAY[index], 2) - d2_N3_122_NON_GRAY_dnorm_f_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dN1_1[index] - pow(N1_2_NON_GRAY[index], 2);
//                                         f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     // N1_1 != 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY[index] = (3.0/2.0)*pow(N1_2_NON_GRAY[index], 2) - dN3_122_NON_GRAY_dnorm_f[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY[index] = (N3_122_NON_GRAY[index] - N1_1_NON_GRAY[index]*pow(N1_2_NON_GRAY[index], 2));
//                                         f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_122_NON_GRAY[index] = f_N3_122_NON_GRAY[index] - gam2_NON_GRAY[index];
//                                 f_N3_122_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam2_NON_GRAY[index];
//                             } else if ((gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 = 0 and gam2 > 0 and gam2 < 1 - gam1
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam1_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam1_dN1_1[index];
//                                         f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY[index] = -d2_N3_122_NON_GRAY_dnorm_f_dgam1[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dgam1[index];
//                                         f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_122_NON_GRAY[index] /= gam2_NON_GRAY[index];
//                             } else if (!(gam1_NON_GRAY[index] < 1.0e-8) && (gam2_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 > 0 and gam1 < 1 and gam2 = 0
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam2_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam2_dN1_1[index];
//                                         f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY[index] = -d2_N3_122_NON_GRAY_dnorm_f_dgam2[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY[index] = dN3_122_NON_GRAY_dgam2[index];
//                                         f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                                 f_N3_122_NON_GRAY[index] = (f_N3_122_NON_GRAY[index] - 1.0)/gam1_NON_GRAY[index];
//                             } else if ((gam1_NON_GRAY[index] < 1.0e-8) && (gam2_NON_GRAY[index] < 1.0e-8)) {
//                                 // gam1 = 0 && gam2 = 0
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY[index] = -d4_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2_dN1_1[index]/2.0;
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY[index] = d3_N3_122_NON_GRAY_dgam1_dgam2_dN1_1[index];
//                                         f_N3_122_NON_GRAY[index] /= (1.0 - norm_f_2);
//                                     }
//                                 } else {
//                                     if (fabs(1.0 - norm_f) < 1.0e-8) {
//                                         // norm_f = 1
//                                         f_N3_122_NON_GRAY[index] = -d3_N3_122_NON_GRAY_dnorm_f_dgam1_dgam2[index]/(2.0*N1_1_NON_GRAY[index]);
//                                     } else {
//                                         // norm_f < 1
//                                         f_N3_122_NON_GRAY[index] = d2_N3_122_NON_GRAY_dgam1_dgam2[index];
//                                         f_N3_122_NON_GRAY[index] /= N1_1_NON_GRAY[index]*(1.0 - norm_f_2);
//                                     }
//                                 }
//                             } else {
//                                 cout << "Error in Setup Interpolant Values for N3 122" << endl;
//                                 exit(0);
//                             }
//                             
//                             // Compute f_N3_123
//                             gam3 = 1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index];
//                             if (!(gam1_NON_GRAY[index] < 1.0e-8) && !(gam2_NON_GRAY[index] < 1.0e-8) && !(fabs(1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index]) < 1.0e-8)) {
//                                 // gam1 > 0 and gam2 > 0 and gam3 > 0
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_2 = 0
//                                         if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_3 = 0
//                                             f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3[index] - 1.0;
//                                         } else {
//                                             // N1_3 >= 0
//                                             f_N3_123_NON_GRAY[index] = (d2_N3_123_NON_GRAY_dN1_1_dN1_2[index]/N1_3_NON_GRAY[index]) - 1.0;  
//                                         }
//                                     } else {
//                                         // N1_2 >= 0
//                                         if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_3 = 0
//                                             f_N3_123_NON_GRAY[index] = (d2_N3_123_NON_GRAY_dN1_1_dN1_3[index]/N1_2_NON_GRAY[index]) - 1.0;
//                                         } else {
//                                             // N1_3 >= 0
//                                             f_N3_123_NON_GRAY[index] = (dN3_123_NON_GRAY_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index])) - 1.0;  
//                                         }
//                                     }
//                                 } else {
//                                     // N1_1 >= 0
//                                     if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_2 = 0
//                                         if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_3 = 0
//                                             f_N3_123_NON_GRAY[index] = (d2_N3_123_NON_GRAY_dN1_2_dN1_3[index]/N1_1_NON_GRAY[index]) - 1.0;
//                                         } else {
//                                             // N1_3 >= 0
//                                             f_N3_123_NON_GRAY[index] = (dN3_123_NON_GRAY_dN1_2[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index])) - 1.0;  
//                                         }
//                                     } else {
//                                         // N1_2 >= 0
//                                         if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_3 = 0
//                                             f_N3_123_NON_GRAY[index] = (dN3_123_NON_GRAY_dN1_3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index])) - 1.0;
//                                         } else {
//                                             // N1_3 >= 0
//                                             f_N3_123_NON_GRAY[index] = (N3_123_NON_GRAY[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index])) - 1.0;  
//                                         }
//                                     }
//                                 }
//                                 f_N3_123_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam2_NON_GRAY[index]*gam3;
//                             } else if (gam1_NON_GRAY[index] < 1.0e-8) {
//                                 // gam1 = 0
//                                 if (gam2_NON_GRAY[index] < 1.0e-8) {
//                                     // gam2 = 0
//                                     if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_1 = 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam2[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam2[index]/N1_3_NON_GRAY[index];  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam2[index]/N1_2_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dgam1_dgam2_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     } else {
//                                         // N1_1 >= 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam2[index]/N1_1_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam2[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam2[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam1_dgam2[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     }
//                                     f_N3_123_NON_GRAY[index] /= gam3;
//                                 } else if (fabs(1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index]) < 1.0e-8) {
//                                     // gam3 = 0
//                                     if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_1 = 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1_dgam3[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1_dgam3[index]/N1_3_NON_GRAY[index];  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1_dgam3[index]/N1_2_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dgam1_dN1_1_dgam3[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     } else {
//                                         // N1_1 >= 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1_dgam3[index]/N1_1_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dgam1_dgam3[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_3_dgam1_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam1_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     }
//                                     f_N3_123_NON_GRAY[index] /= gam2_NON_GRAY[index];
//                                 } else {
//                                     // gam1 = 0 && gam2 > 0 && gam3 > 0
//                                     if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_1 = 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam1[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam1[index]/N1_3_NON_GRAY[index];  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam1[index]/N1_2_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam1_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     } else {
//                                         // N1_1 >= 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam1[index]/N1_1_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_2_dgam1[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_3_dgam1[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dgam1[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     }
//                                     f_N3_123_NON_GRAY[index] /= gam2_NON_GRAY[index]*gam3;
//                                 }
//                             } else if (gam2_NON_GRAY[index] < 1.0e-8) {
//                                 // gam1 > 0 and gam2 = 0
//                                 if (fabs(1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index]) < 1.0e-8) {
//                                     // gam3 = 0
//                                     if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_1 = 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d5_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2_dgam3[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2_dgam3[index]/N1_3_NON_GRAY[index];  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2_dgam3[index]/N1_2_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dgam2_dN1_1_dgam3[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     } else {
//                                         // N1_1 >= 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2_dgam3[index]/N1_1_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dgam2_dgam3[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_3_dgam2_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam2_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     }
//                                     f_N3_123_NON_GRAY[index] /= gam1_NON_GRAY[index];
//                                 } else {
//                                     // gam1 > 0 && gam2 = 0 && gam3 > 0
//                                     if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_1 = 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam2[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam2[index]/N1_3_NON_GRAY[index];  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam2[index]/N1_2_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dgam2_dN1_1[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     } else {
//                                         // N1_1 >= 0
//                                         if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_2 = 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam2[index]/N1_1_NON_GRAY[index];
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_2_dgam2[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         } else {
//                                             // N1_2 >= 0
//                                             if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                                 // N1_3 = 0
//                                                 f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_3_dgam2[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                             } else {
//                                                 // N1_3 >= 0
//                                                 f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dgam2[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                             }
//                                         }
//                                     }
//                                     f_N3_123_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam3;
//                                 }
//                             } else if (fabs(1.0 - gam1_NON_GRAY[index] - gam2_NON_GRAY[index]) < 1.0e-8) {
//                                 // gam1 > 0 and gam2 > 0 and gam3 = 0
//                                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                                     // N1_1 = 0
//                                     if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_2 = 0
//                                         if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_3 = 0
//                                             f_N3_123_NON_GRAY[index] = d4_N3_123_NON_GRAY_dN1_1_dN1_2_dN1_3_dgam3[index];
//                                         } else {
//                                             // N1_3 >= 0
//                                             f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_2_dgam3[index]/N1_3_NON_GRAY[index];  
//                                         }
//                                     } else {
//                                         // N1_2 >= 0
//                                         if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_3 = 0
//                                             f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_1_dN1_3_dgam3[index]/N1_2_NON_GRAY[index];
//                                         } else {
//                                             // N1_3 >= 0
//                                             f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_1_dgam3[index]/(N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                         }
//                                     }
//                                 } else {
//                                     // N1_1 >= 0
//                                     if (fabs(N1_2_NON_GRAY[index]) < 1.0e-8) {
//                                         // N1_2 = 0
//                                         if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_3 = 0
//                                             f_N3_123_NON_GRAY[index] = d3_N3_123_NON_GRAY_dN1_2_dN1_3_dgam3[index]/N1_1_NON_GRAY[index];
//                                         } else {
//                                             // N1_3 >= 0
//                                             f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_2_dgam3[index]/(N1_1_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                         }
//                                     } else {
//                                         // N1_2 >= 0
//                                         if (fabs(N1_3_NON_GRAY[index]) < 1.0e-8) {
//                                             // N1_3 = 0
//                                             f_N3_123_NON_GRAY[index] = d2_N3_123_NON_GRAY_dN1_3_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]);
//                                         } else {
//                                             // N1_3 >= 0
//                                             f_N3_123_NON_GRAY[index] = dN3_123_NON_GRAY_dgam3[index]/(N1_1_NON_GRAY[index]*N1_2_NON_GRAY[index]*N1_3_NON_GRAY[index]);  
//                                         }
//                                     }
//                                 }
//                                 f_N3_123_NON_GRAY[index] /= gam1_NON_GRAY[index]*gam2_NON_GRAY[index];
//                             } else {
//                                 cout << "Error in Setup Interpolant Values for N3 123" << endl;
//                                 exit(0);
//                             }
//                             
//                             // cout /*<< "index = " << index << "   "*/ << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "   " << "N3_111 = " << N3_111_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY[index] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index] << "  " << "dN3_122 = " << d2_N3_122_NON_GRAY_dgam2_dN1_1[index] << endl;
//                             
//                             //                                                     if (fabs(f_N3_111_NON_GRAY[index]) > 1.0 || fabs(f_N3_122_NON_GRAY[index]) > 1.0) {
//                             if (f_N3_111_NON_GRAY[index] != f_N3_111_NON_GRAY[index]) {
//                                 cout /*<< "index = " << index << "   "*/ << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "N1_2 = " << N1_2_NON_GRAY[index] << "   " << "N1_3 = " << N1_3_NON_GRAY[index] << "   " << "gam1 = " << gam1_NON_GRAY[index] << "   " << "gam2 = " << gam2_NON_GRAY[index] << "  " << "f_N3_111 = " << f_N3_111_NON_GRAY[index] << "  " << "N3_122 = " << N3_122_NON_GRAY[index] << "  " << "f_N3_122 = " << f_N3_122_NON_GRAY[index] << "  " << "dN3_111_dN1_1 = " << dN3_111_NON_GRAY_dN1_1[index] << "  " << "dN3_122_dN1_1 = " << dN3_122_NON_GRAY_dN1_1[index] << endl;
//                                 exit(0);
//                             }
//                         } // end for i_Cheby_gam1_gam2
//                     } // end for i_Theta
//                 } // end for i_Phi
//             } // end for i_Cheby_f
//         } // end for i_Cheby_E
//     } // end for id_Mobius
// }
