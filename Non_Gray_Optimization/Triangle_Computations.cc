#ifndef _NG_MN_Model_3D_OPTIM_H_INCLUDED
#include "NG_MN_Model_3D_OPTIM.h"
#endif // _NG_MN_Model_3D_OPTIM_H_INCLUDED

//****************************************************************************************
// The following routines compute the parameters for the equation of the line passing 
// through the barycenter of the triangle (P1 P2 P3) and any point within the triangle
// which does not coincide with the barycenter
//****************************************************************************************
long double a_median(const M2_State_Param &M2_State) {
    long double a_value;
    
    a_value = (1.0 - 3.0*M2_State.gamma_2)/(1.0 - 3.0*M2_State.gamma_1);
    
    return a_value;
}

long double b_median(const M2_State_Param &M2_State) {
    long double b_value;
    
    b_value = (M2_State.gamma_2 - M2_State.gamma_1)/(1.0 - 3.0*M2_State.gamma_1);
    
    return b_value;
}

//****************************************************************************************
// This routines computes the coordinates of the points corresponding to the intersection
// between the line passing through the barycenter of the triangle and any point along one 
// of the edge of the triangle with the opposite edge of the triangle
//****************************************************************************************
void Triangle_Median_Edge_Intersection_Points(M2_State_Param &M2_State) {
    // FACE_A ==> (P2 P3)
    // FACE_B ==> (P1 P3)
    // FACE_C ==> (P1 P2)
    
    long double a_median, b_median;
    a_median = M2_State.a_median;
    b_median = M2_State.b_median;
    // Equation of the line is : gam2 = a * gam1 + b
    // We can plug gam1 = 0 in the latter, yielding, gam2 = b
    // This means that if 0 <= b <= 1, then the latter is feasible and there is
    // therefore interesection with the edge (P2 P3) ==> FACE_A
    
    if (b_median >= 0.0 && b_median <= 1.0) {
        // Then we have an intersection with face A : gam1 = 0
        M2_State.gam1_bound_min = 0.0;
        M2_State.gam2_bound_min = b_median;
        
        if (b_median <= 0.5) {
            // Then we also have an intersection with face C : gam2 = 1 - gam1
            M2_State.gam1_bound_max = (1.0 - b_median)/(1.0 + a_median);
            M2_State.gam2_bound_max = (a_median + b_median)/(1.0 + a_median);
        } else {
            // Then we also have an intersection with face B : gam2 = 0
            M2_State.gam1_bound_max = -b_median/a_median;
            M2_State.gam2_bound_max = 0.0;
        }
    } else {
        // Then we have an intersection with face B : gam2 = 0
        M2_State.gam1_bound_min = -b_median/a_median;
        M2_State.gam2_bound_min = 0.0;
        
        if (M2_State.gam1_bound_min <= 0.5) {
            // Then we also have an intersection with face C : gam2 = 1 - gam1
            M2_State.gam1_bound_max = (1.0 - b_median)/(1.0 + a_median);
            M2_State.gam2_bound_max = (a_median + b_median)/(1.0 + a_median);
        } else {
            // Then we also have an intersection with face A : gam1 = 0
            M2_State.gam1_bound_max = 0.0;
            M2_State.gam2_bound_max = b_median;
        }
    }
}

//****************************************************************************************
// This routines computes the components of the unit vector in the direction of the
// line passing through the barycenter of the triangle and any point along one 
// of the edge of the triangle
//****************************************************************************************
void Median_Line_Unit_Vector(long double &u1, long double &u2, M2_State_Param &M2_State) {
    // The equation of the line is: gam2 = a_median * gam1 + b_median
    long double hypothenuse;
    long double gam1_G, gam2_G, gam3_G;
    
    gam1_G = 1.0/3.0; 
    gam2_G = 1.0/3.0;
    gam3_G = 1.0 - gam1_G - gam2_G;
    
    hypothenuse = pow(gam1_G - M2_State.gam1_knot, 2) + pow(gam2_G - M2_State.gam2_knot, 2);
    hypothenuse = sqrt(hypothenuse);
    
    u1 = (gam1_G - M2_State.gam1_knot)/hypothenuse;
    u2 = (gam2_G - M2_State.gam2_knot)/hypothenuse;
}

void Median_Line_Unit_Vector_Plane_Gam1_Gam3(long double &u1, long double &u2, M2_State_Param &M2_State) {
    // The equation of the line is: gam2 = a_median * gam1 + b_median
    long double hypothenuse;
    long double gam1_G, gam2_G, gam3_G;
    long double gam3_knot;
    
    gam1_G = 1.0/3.0; 
    gam2_G = 1.0/3.0;
    gam3_G = 1.0 - gam1_G - gam2_G;
    
    gam3_knot = 1.0 - M2_State.gam1_knot - M2_State.gam2_knot;
    
    hypothenuse = pow(gam1_G - M2_State.gam1_knot, 2) + pow(gam3_G - gam3_knot, 2);
    hypothenuse = sqrt(hypothenuse);
    
    u1 = (gam1_G - M2_State.gam1_knot)/hypothenuse;
    u2 = (gam3_G - gam3_knot)/hypothenuse;
}
