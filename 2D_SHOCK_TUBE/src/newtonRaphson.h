#pragma once

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include "math.h"
#include "node.h"
#include "inverse.h"
#include "matmul.h"
#include "grid.h"
#include "lbmkernerls.h"
//#define CHUNK_SIZE 1

double func1(node& p, double cv){
    //Lagrangian Multipliers
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2];                                              
    //constraint1: Total energy
    double f1  = p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) +\
                 p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) +\
                 p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) +\
                 p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) +\
                 p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) +\
                 p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) +\
                 p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) +\
                 p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) +\
                 p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) ) +\
                (-2.0 * p.rho  * (p.T * cv + (p.ux*p.ux + p.uy*p.uy) * 0.5));

    return f1;
}

double func2(node& p, double cv){
    //Lagrangian Multipliers
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2];                                                

    double f2  =  p.cx[0] * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) +\
                  p.cx[1] * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) +\
                  p.cx[2] * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) +\
                  p.cx[3] * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) +\
                  p.cx[4] * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) +\
                  p.cx[5] * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) +\
                  p.cx[6] * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) +\
                  p.cx[7] * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) +\
                  p.cx[8] * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) ) +\
                  (-2.0  * p.rho * p.ux   * (p.T * cv + (p.ux*p.ux + p.uy*p.uy) * 0.5 + p.T));

    return f2;
}

double func3(const node& p, const double cv){
    //Lagrangian Multipliers
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2];                                             
    //constraint1: Trace of third order Maxell Boltzmann MOment
    double f3  =  p.cy[0] * p.rho * p.w[0] * exp(1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0])) + \
                  p.cy[1] * p.rho * p.w[1] * exp(1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1])) + \
                  p.cy[2] * p.rho * p.w[2] * exp(1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2])) + \
                  p.cy[3] * p.rho * p.w[3] * exp(1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3])) + \
                  p.cy[4] * p.rho * p.w[4] * exp(1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4])) + \
                  p.cy[5] * p.rho * p.w[5] * exp(1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5])) + \
                  p.cy[6] * p.rho * p.w[6] * exp(1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6])) + \
                  p.cy[7] * p.rho * p.w[7] * exp(1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7])) + \
                  p.cy[8] * p.rho * p.w[8] * exp(1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8])) + \
                  (-2.0  * p.rho * p.uy   * (p.T * cv + (p.ux*p.ux + p.uy*p.uy) * 0.5 + p.T));

    return f3;
}

double j11(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 

    double j1  = 1.0 * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) +\
                 1.0 * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) +\
                 1.0 * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) +\
                 1.0 * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) +\
                 1.0 * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) +\
                 1.0 * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) +\
                 1.0 * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) +\
                 1.0 * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) +\
                 1.0 * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) );
   return j1;
}

double j12(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    // double j2 = 2*x2;
    double j2  = 1.0 * p.cx[0] * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) +\
                 1.0 * p.cx[1] * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) +\
                 1.0 * p.cx[2] * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) +\
                 1.0 * p.cx[2] * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) +\
                 1.0 * p.cx[4] * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) +\
                 1.0 * p.cx[5] * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) +\
                 1.0 * p.cx[6] * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) +\
                 1.0 * p.cx[7] * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) +\
                 1.0 * p.cx[8] * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) );
   return j2;
}

double j13(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 

    double j3  = 1.0 * p.cy[0] * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) +\
                 1.0 * p.cy[1] * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) +\
                 1.0 * p.cy[2] * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) +\
                 1.0 * p.cy[3] * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) +\
                 1.0 * p.cy[4] * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) +\
                 1.0 * p.cy[5] * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) +\
                 1.0 * p.cy[6] * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) +\
                 1.0 * p.cy[7] * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) +\
                 1.0 * p.cy[8] * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) );
    return j3;
}

double j21(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    // double j3 = 2*x1;
    double j3  = 1.0 * p.cx[0] * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) +\
                 1.0 * p.cx[1] * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) +\
                 1.0 * p.cx[2] * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) +\
                 1.0 * p.cx[3] * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) +\
                 1.0 * p.cx[4] * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) +\
                 1.0 * p.cx[5] * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) +\
                 1.0 * p.cx[6] * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) +\
                 1.0 * p.cx[7] * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) +\
                 1.0 * p.cx[8] * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) );
    return j3;
}

double j22(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    // double j4 = -1.0;
    double j4  = 1.0 * p.cx[0] * p.cx[0] * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) +\
                 1.0 * p.cx[1] * p.cx[1] * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) +\
                 1.0 * p.cx[2] * p.cx[2] * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) +\
                 1.0 * p.cx[3] * p.cx[3] * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) +\
                 1.0 * p.cx[4] * p.cx[4] * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) +\
                 1.0 * p.cx[5] * p.cx[5] * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) +\
                 1.0 * p.cx[6] * p.cx[6] * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) +\
                 1.0 * p.cx[7] * p.cx[7] * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) +\
                 1.0 * p.cx[8] * p.cx[8] * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) );
    return j4;
}

double j23(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 

    double j6  = 1.0 * p.cx[0] * p.cy[0] * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) + \
                 1.0 * p.cx[1] * p.cy[1] * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) + \
                 1.0 * p.cx[2] * p.cy[2] * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) + \
                 1.0 * p.cx[3] * p.cy[3] * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) + \
                 1.0 * p.cx[4] * p.cy[4] * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) + \
                 1.0 * p.cx[5] * p.cy[5] * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) + \
                 1.0 * p.cx[6] * p.cy[6] * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) + \
                 1.0 * p.cx[7] * p.cy[7] * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) + \
                 1.0 * p.cx[8] * p.cy[8] * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) );
    return j6;
}

double j31(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 

    double j7  = 1.0 * p.cy[0] * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) + \
                 1.0 * p.cy[1] * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) + \
                 1.0 * p.cy[2] * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) + \
                 1.0 * p.cy[3] * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) + \
                 1.0 * p.cy[4] * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) + \
                 1.0 * p.cy[5] * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) + \
                 1.0 * p.cy[6] * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) + \
                 1.0 * p.cy[7] * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) + \
                 1.0 * p.cy[8] * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) );
    return j7;
}

double j32(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 

    double j8  = 1.0 * p.cx[0] * p.cy[0] * p.rho * p.w[0] * exp(1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0])) + \
                 1.0 * p.cx[1] * p.cy[1] * p.rho * p.w[1] * exp(1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1])) + \
                 1.0 * p.cx[2] * p.cy[2] * p.rho * p.w[2] * exp(1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2])) + \
                 1.0 * p.cx[3] * p.cy[3] * p.rho * p.w[3] * exp(1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3])) + \
                 1.0 * p.cx[4] * p.cy[4] * p.rho * p.w[4] * exp(1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4])) + \
                 1.0 * p.cx[5] * p.cy[5] * p.rho * p.w[5] * exp(1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5])) + \
                 1.0 * p.cx[6] * p.cy[6] * p.rho * p.w[6] * exp(1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6])) + \
                 1.0 * p.cx[7] * p.cy[7] * p.rho * p.w[7] * exp(1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7])) + \
                 1.0 * p.cx[8] * p.cy[8] * p.rho * p.w[8] * exp(1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]));
    return j8;
}

double j33(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 

    double j9  = 1.0 * p.cy[0] * p.cy[0] * p.rho * p.w[0] * exp( 1.0*(0.0 + x1 + x2*p.cx[0] + x3*p.cy[0]) ) + \
                 1.0 * p.cy[1] * p.cy[1] * p.rho * p.w[1] * exp( 1.0*(0.0 + x1 + x2*p.cx[1] + x3*p.cy[1]) ) + \
                 1.0 * p.cy[2] * p.cy[2] * p.rho * p.w[2] * exp( 1.0*(0.0 + x1 + x2*p.cx[2] + x3*p.cy[2]) ) + \
                 1.0 * p.cy[3] * p.cy[3] * p.rho * p.w[3] * exp( 1.0*(0.0 + x1 + x2*p.cx[3] + x3*p.cy[3]) ) + \
                 1.0 * p.cy[4] * p.cy[4] * p.rho * p.w[4] * exp( 1.0*(0.0 + x1 + x2*p.cx[4] + x3*p.cy[4]) ) + \
                 1.0 * p.cy[5] * p.cy[5] * p.rho * p.w[5] * exp( 1.0*(0.0 + x1 + x2*p.cx[5] + x3*p.cy[5]) ) + \
                 1.0 * p.cy[6] * p.cy[6] * p.rho * p.w[6] * exp( 1.0*(0.0 + x1 + x2*p.cx[6] + x3*p.cy[6]) ) + \
                 1.0 * p.cy[7] * p.cy[7] * p.rho * p.w[7] * exp( 1.0*(0.0 + x1 + x2*p.cx[7] + x3*p.cy[7]) ) + \
                 1.0 * p.cy[8] * p.cy[8] * p.rho * p.w[8] * exp( 1.0*(0.0 + x1 + x2*p.cx[8] + x3*p.cy[8]) );
    return j9;
}

void allfuncs(const node& p, const double cv, double* tmp){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2];                                              
    
    double j1  = 0.0, j2  = 0.0, j3  = 0.0, j4  = 0.0, j5  = 0.0, j6  = 0.0, j7  = 0.0, j8  = 0.0, j9  = 0.0;                                   
    
    double f1  = (-2.0 * p.rho  * (p.T * cv + (p.ux*p.ux + p.uy*p.uy) * 0.5));
    double f2  =  (-2.0  * p.rho * p.ux   * (p.T * cv + (p.ux*p.ux + p.uy*p.uy) * 0.5 + p.T));
    double f3  = (-2.0  * p.rho * p.uy   * (p.T * cv + (p.ux*p.ux + p.uy*p.uy) * 0.5 + p.T));

    for(int i = 0; i < 9; ++i){ 
        double tmp = p.rho * p.w[i] * exp( 1.0*(0.0 + x1 + x2*p.cx[i] + x3*p.cy[i]) );
        f1 += tmp;
        f2 += p.cx[i] * tmp;
        f3 += p.cy[i] * tmp;
        j1 += 1.0 * tmp;
        j2 += 1.0 * p.cx[i] * tmp;
        j3 += 1.0 * p.cy[i] * tmp;
        j4 += 1.0 * p.cx[i] * tmp;
        j5 += 1.0 * p.cx[i] * p.cx[i] * tmp;
        j6 += 1.0 * p.cx[i] * p.cy[i] * tmp;
        j7 += 1.0 * p.cy[i] * tmp;
        j8 += 1.0 * p.cx[i] * p.cy[i] * tmp;
        j9 += 1.0 * p.cy[i] * p.cy[i] * tmp;
    }

    tmp[0] = f1;
    tmp[1] = f2;
    tmp[2] = f3;
    tmp[3] = j1;
    tmp[4] = j2;
    tmp[5] = j3;
    tmp[6] = j4;
    tmp[7] = j5;
    tmp[8] = j6;
    tmp[9] = j7;
    tmp[10] = j8;
    tmp[11] = j9;
}

void norm(double& err, double* a1, double* a2){

    err  = (a1[0] - a2[0]) * (a1[0] - a2[0]) +\
           (a1[1] - a2[1]) * (a1[1] - a2[1]) +\
           (a1[2] - a2[2]) * (a1[2] - a2[2]) ;
    err = pow(err,0.5);
}

//Newton Raphson Algorithm
// void newtonRaphson(grid& g,const double& cv,int flagIndex,double tol=1e-14)
void newtonRaphson(grid& g,state d,double tol=1e-14){
    #pragma omp parallel
    {
        #pragma omp for schedule(static) nowait
        for (int i = 0; i < g.NY * g.NX; ++i){
            if(g.domain[i].flag == d.fluid_flag){

                double err = 0.1;
                int iter   = 0;
                int MAX_ITER = 5000;
                double F[3]{0.0, 0.0, 0.0};
                double J[9]{0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0};
                Matmul mul(3);
                Inverse D(9); //create Inverse 3*3 = 9 array
                while ( (err > tol)  && (iter<MAX_ITER) ){//

                    double tmp_res[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
                    allfuncs(g.domain[i],d.CV, tmp_res);
                    F[0] = tmp_res[0];
                    F[1] = tmp_res[1];
                    F[2] = tmp_res[2];
                    J[0] = tmp_res[3];
                    J[1] = tmp_res[4];
                    J[2] = tmp_res[5];
                    J[3] = tmp_res[6];
                    J[4] = tmp_res[7];
                    J[5] = tmp_res[8];
                    J[6] = tmp_res[9];
                    J[7] = tmp_res[10];
                    J[8] = tmp_res[11];

                    // F[0] = func1(g.domain[i],d.CV);
                    // F[1] = func2(g.domain[i],d.CV);
                    // F[2] = func3(g.domain[i],d.CV); 
                    // J[0] = j11(g.domain[i]);
                    // J[1] = j12(g.domain[i]);
                    // J[2] = j13(g.domain[i]);
                    // J[3] = j21(g.domain[i]);
                    // J[4] = j22(g.domain[i]);
                    // J[5] = j23(g.domain[i]);
                    // J[6] = j31(g.domain[i]);
                    // J[7] = j32(g.domain[i]);
                    // J[8] = j33(g.domain[i]);

                    D.assign(J); //assign jacobian to inverse object
                
                    D.Det();     //get the determinant
                
                    if(isnan(D.DET)){
                        std::cout << ">>>>>>>>>>>>>>Error: SINGULAR MATRIC" << std::endl;
                    }

                    D.Inv();

                    mul.Mul(D.INV,F);
                    F[0] = g.domain[i].LM[0] - mul.res[0];
                    F[1] = g.domain[i].LM[1] - mul.res[1];
                    F[2] = g.domain[i].LM[2] - mul.res[2];
                    norm(err,F,g.domain[i].LM);      

                    g.domain[i].LM[0] = F[0];
                    g.domain[i].LM[1] = F[1];
                    g.domain[i].LM[2] = F[2];

                    iter+=1;

                }
                g.domain[i].nr = iter;
            }
        }
    }

}
