#pragma once

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include "math.h"
#include "node.h"
// #include "inverse.h"
#include "matmul.h"
#include "grid.h"
#include "lbmkernerls.h"
#include <gsl/gsl_linalg.h>


/*Residual function definition*/
/*constraint1 : density*/
double func0(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            
    /*constraint1: density*/
    double fEqSum = 0.0;
    for( int q = 0; q <21 ; ++q){
        fEqSum = fEqSum + p.rho * exp( -1.0 * ( 1.0 + x0  
                                                    + x1 * p.cx[q]
                                                    + x2 * p.cy[q]
                                                    + x3 * p.cx[q] * p.cx[q]
                                                    + x4 * p.cx[q] * p.cy[q]
                                                    + x5 * p.cy[q] * p.cy[q]
                                                    + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                    + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f0 = fEqSum - p.rho;
    // std::cout << f0 << std::endl;
    return f0;
}
/*constraint2 : x momentum*/
double func1(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum - p.rho*p.ux;
    // std::cout << f1 << std::endl; 
    return f1;
}
/*constraint3 : y momentum*/
double func2(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum - p.rho*p.uy;
    // std::cout << f2 << std::endl;
    return f2;
}
/*constraint4 : second order moment*/
double func3(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + p.cx[q] * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum - p.rho * (p.ux*p.ux + p.T);
    return f3;
}
/*constraint5 : second order moment*/
double func4(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + p.cx[q] * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum - p.rho * (p.ux*p.uy);
    return f4;
}
/*constraint6 : second order moment*/
double func5(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + p.cy[q] * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum - p.rho * (p.uy*p.uy + p.T);
    return f5;
}
/*constraint7 : second order moment*/
double func6(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum - 2.0 * p.rho * p.ux * ( (0.5 * (p.ux*p.ux + p.uy*p.uy + 2.0*p.T)) + p.T );
    return f6;
}
/*constraint8 : second order moment*/
double func7(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum - 2.0 * p.rho * p.uy * ( (0.5 * (p.ux*p.ux + p.uy*p.uy + 2.0*p.T)) + p.T );
    return f7;
}

/*first roq jacobian*/
double j00(node& p){
    double x0 = p.LM[0]; double x1 = p.LM[1]; double x2 = p.LM[2]; double x3 = p.LM[3]; double x4 = p.LM[4]; double x5 = p.LM[5]; double x6 = p.LM[6]; double x7 = p.LM[7];  
    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.rho * exp( -1.0 * ( 1.0 + x0  
                                                                + x1 * p.cx[q]
                                                                + x2 * p.cy[q]
                                                                + x3 * p.cx[q] * p.cx[q]
                                                                + x4 * p.cx[q] * p.cy[q]
                                                                + x5 * p.cy[q] * p.cy[q]
                                                                + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double j0 = fEqSum;
    return j0;
}

double j01(node& p){
    double x0 = p.LM[0]; double x1 = p.LM[1]; double x2 = p.LM[2]; double x3 = p.LM[3]; double x4 = p.LM[4]; double x5 = p.LM[5]; double x6 = p.LM[6]; double x7 = p.LM[7];  
    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0  
                                                                            + x1 * p.cx[q]
                                                                            + x2 * p.cy[q]
                                                                            + x3 * p.cx[q] * p.cx[q]
                                                                            + x4 * p.cx[q] * p.cy[q]
                                                                            + x5 * p.cy[q] * p.cy[q]
                                                                            + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                            + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double j1 = fEqSum;
    return j1;
}

double j02(node& p){
    double x0 = p.LM[0]; double x1 = p.LM[1]; double x2 = p.LM[2]; double x3 = p.LM[3]; double x4 = p.LM[4]; double x5 = p.LM[5]; double x6 = p.LM[6]; double x7 = p.LM[7];  
    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0  
                                                                            + x1 * p.cx[q]
                                                                            + x2 * p.cy[q]
                                                                            + x3 * p.cx[q] * p.cx[q]
                                                                            + x4 * p.cx[q] * p.cy[q]
                                                                            + x5 * p.cy[q] * p.cy[q]
                                                                            + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                            + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double j2 = fEqSum;
    return j2;
}

double j03(node& p){
    double x0 = p.LM[0]; double x1 = p.LM[1]; double x2 = p.LM[2]; double x3 = p.LM[3]; double x4 = p.LM[4]; double x5 = p.LM[5]; double x6 = p.LM[6]; double x7 = p.LM[7];  
    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * (p.cx[q] * p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                            + x1 * p.cx[q]
                                                                                            + x2 * p.cy[q]
                                                                                            + x3 * p.cx[q] * p.cx[q]
                                                                                            + x4 * p.cx[q] * p.cy[q]
                                                                                            + x5 * p.cy[q] * p.cy[q]
                                                                                            + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                            + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double j3 = fEqSum;
    return j3;
}

double j04(node& p){
    double x0 = p.LM[0]; double x1 = p.LM[1]; double x2 = p.LM[2]; double x3 = p.LM[3]; double x4 = p.LM[4]; double x5 = p.LM[5]; double x6 = p.LM[6]; double x7 = p.LM[7];  
    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                          + x1 * p.cx[q]
                                                                          + x2 * p.cy[q]
                                                                          + x3 * p.cx[q] * p.cx[q]
                                                                          + x4 * p.cx[q] * p.cy[q]
                                                                          + x5 * p.cy[q] * p.cy[q]
                                                                          + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                          + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double j4 = fEqSum;
    return j4;
}

double j05(node& p){
    double x0 = p.LM[0]; double x1 = p.LM[1]; double x2 = p.LM[2]; double x3 = p.LM[3]; double x4 = p.LM[4]; double x5 = p.LM[5]; double x6 = p.LM[6]; double x7 = p.LM[7];  
    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * (p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                          + x1 * p.cx[q]
                                                                          + x2 * p.cy[q]
                                                                          + x3 * p.cx[q] * p.cx[q]
                                                                          + x4 * p.cx[q] * p.cy[q]
                                                                          + x5 * p.cy[q] * p.cy[q]
                                                                          + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                          + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double j5 = fEqSum;
    return j5;
}

double j06(node& p){
    double x0 = p.LM[0]; double x1 = p.LM[1]; double x2 = p.LM[2]; double x3 = p.LM[3]; double x4 = p.LM[4]; double x5 = p.LM[5]; double x6 = p.LM[6]; double x7 = p.LM[7];  
    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * (p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                         + x1 * p.cx[q]
                                                                                                         + x2 * p.cy[q]
                                                                                                         + x3 * p.cx[q] * p.cx[q]
                                                                                                         + x4 * p.cx[q] * p.cy[q]
                                                                                                         + x5 * p.cy[q] * p.cy[q]
                                                                                                         + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                         + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double j6 = fEqSum;
    return j6;
}

double j07(node& p){
    double x0 = p.LM[0]; double x1 = p.LM[1]; double x2 = p.LM[2]; double x3 = p.LM[3]; double x4 = p.LM[4]; double x5 = p.LM[5]; double x6 = p.LM[6]; double x7 = p.LM[7];  
    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * (p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                         + x1 * p.cx[q]
                                                                                                         + x2 * p.cy[q]
                                                                                                         + x3 * p.cx[q] * p.cx[q]
                                                                                                         + x4 * p.cx[q] * p.cy[q]
                                                                                                         + x5 * p.cy[q] * p.cy[q]
                                                                                                         + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                         + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double j7 = fEqSum;
    return j7;
}
/*seond row jacobian*/
double j10(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum ;
    return f1;
}

double j11(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum;
    return f1;
}

double j12(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum ;
    return f1;
}

double j13(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum ;
    return f1;
}

double j14(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum ;
    return f1;
}

double j15(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum ;
    return f1;
}

double j16(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                                        + x1 * p.cx[q]
                                                                                                                        + x2 * p.cy[q]
                                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum ;
    return f1;
}

double j17(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                                        + x1 * p.cx[q]
                                                                                                                        + x2 * p.cy[q]
                                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f1 = fEqSum ;
    return f1;
}
/*third row jacobain*/
double j20(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum;
    return f2;
}

double j21(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                            + x1 * p.cx[q]
                                                                            + x2 * p.cy[q]
                                                                            + x3 * p.cx[q] * p.cx[q]
                                                                            + x4 * p.cx[q] * p.cy[q]
                                                                            + x5 * p.cy[q] * p.cy[q]
                                                                            + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                            + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum;
    return f2;
}

double j22(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum;
    return f2;
}

double j23(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum;
    return f2;
}

double j24(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum;
    return f2;
}

double j25(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q]* (p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum;
    return f2;
}

double j26(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q]* (p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum;
    return f2;
}

double j27(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                              + x1 * p.cx[q]
                                                              + x2 * p.cy[q]
                                                              + x3 * p.cx[q] * p.cx[q]
                                                              + x4 * p.cx[q] * p.cy[q]
                                                              + x5 * p.cy[q] * p.cy[q]
                                                              + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                              + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f2 = fEqSum;
    return f2;
}
/*fourth row*/
double j30(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum ;
    return f3;
}

double j31(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * (p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum ;
    return f3;
}

double j32(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum ;
    return f3;
}

double j33(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * (p.cx[q] * p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum ;
    return f3;
}

double j34(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum ;
    return f3;
}

double j35(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * (p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum ;
    return f3;
}

double j36(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * (p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum ;
    return f3;
}

double j37(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cx[q] * (p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f3 = fEqSum ;
    return f3;
}
/*fifth row*/
double j40(node& p ){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum ;
    return f4;
}

double j41(node& p ){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cy[q] * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum ;
    return f4;
}

double j42(node& p ){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cy[q] * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum ;
    return f4;
}

double j43(node& p ){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cy[q] * (p.cx[q] * p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum ;
    return f4;
}

double j44(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cy[q] * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum ;
    return f4;
}

double j45(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cy[q] * (p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum ;
    return f4;
}

double j46(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cy[q] * (p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum ;
    return f4;
}

double j47(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * p.cy[q] * (p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f4 = fEqSum ;
    return f4;
}
/*sixt row*/
double j50(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum ;
    return f5;
}

double j51(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum ;
    return f5;
}

double j52(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum ;
    return f5;
}

double j53(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * (p.cx[q] * p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum ;
    return f5;
}

double j54(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum ;
    return f5;
}

double j55(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * (p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum ;
    return f5;
}

double j56(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * (p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum ;
    return f5;
}

double j57(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * p.cy[q] * (p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                        + x1 * p.cx[q]
                                                                        + x2 * p.cy[q]
                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f5 = fEqSum ;
    return f5;
}
/*seventh row*/
double j60(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum ;
    return f6;
}

double j61(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum ;
    return f6;
}

double j62(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum ;
    return f6;
}

double j63(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cx[q] * p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum ;
    return f6;
}

double j64(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum ;
    return f6;
}

double j65(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum ;
    return f6;
}

double j66(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum ;
    return f6;
}

double j67(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cx[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f6 = fEqSum ;
    return f6;
}
/*eighth row*/
double j70(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum ;
    return f7;
}

double j71(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.cx[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum ;
    return f7;
}

double j72(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.cy[q] * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum ;
    return f7;
}

double j73(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cx[q] * p.cx[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum ;
    return f7;
}

double j74(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cx[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum ;
    return f7;
}

double j75(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum ;
    return f7;
}

double j76(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * (p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum ;
    return f7;
}

double j77(node& p){
    /*Lagrangian Multipliers*/
    double x0 = p.LM[0];
    double x1 = p.LM[1];
    double x2 = p.LM[2];
    double x3 = p.LM[3];
    double x4 = p.LM[4];
    double x5 = p.LM[5];
    double x6 = p.LM[6];
    double x7 = p.LM[7];                                            

    double fEqSum = 0;
    for( int q = 0; q < 21; ++q){
        fEqSum = fEqSum + (-1.0) * p.cy[q] * (p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) * p.rho * exp( -1.0 * ( 1.0 + x0 + 
                                                                                                        + x1 * p.cx[q]
                                                                                                        + x2 * p.cy[q]
                                                                                                        + x3 * p.cx[q] * p.cx[q]
                                                                                                        + x4 * p.cx[q] * p.cy[q]
                                                                                                        + x5 * p.cy[q] * p.cy[q]
                                                                                                        + x6 * p.cx[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q])
                                                                                                        + x7 * p.cy[q] *(p.cx[q] * p.cx[q] + p.cy[q] * p.cy[q]) ) );
    }
    double f7 = fEqSum ;
    return f7;
}

void norm(double& err, double* a1, double* a2){
    /*3d*/
    err  = (a1[0] - a2[0]) * (a1[0] - a2[0]) +\
           (a1[1] - a2[1]) * (a1[1] - a2[1]) +\
           (a1[2] - a2[2]) * (a1[2] - a2[2]) +\
           (a1[3] - a2[3]) * (a1[3] - a2[3]) +\
           (a1[4] - a2[4]) * (a1[4] - a2[4]) +\
           (a1[5] - a2[5]) * (a1[5] - a2[5]) +\
           (a1[6] - a2[6]) * (a1[6] - a2[6]) +\
           (a1[7] - a2[7]) * (a1[7] - a2[7]) ;

    err = pow(err,0.5);
}

/*USING GSL*/

void inv8(int n,double* INV8,double* j)
{
    gsl_matrix_view m  = gsl_matrix_view_array(j, n, n);
    double temp[n*n];
    // double tempx[n];
    int s;
    // gsl_matrix_view  a  =  gsl_matrix_view_array(tempx,n,1);
    gsl_matrix_view inv = gsl_matrix_view_array(temp,n,n);
    gsl_permutation *p = gsl_permutation_alloc(n);
    gsl_linalg_LU_decomp (&m.matrix, p, &s);    
    // double det = gsl_linalg_LU_det(&m.matrix,s);
    // std::cout << "determinant : "<< 1.0/det << "\n";
    // double det = gsl_matrix_get 
    // det = linalg_LU_det (LU, s)
    gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);
    // int k =0;
    for (int i = 0; i < n; ++i){
        for (int j= 0; j < n; ++j){
            INV8[n*i+j]  = gsl_matrix_get(&inv.matrix,i,j);
            // printf(j==2?"%6.3f\n":"%6.3f ",gsl_matrix_get(&inv.matrix,i,j));
            // k=k+1;
        }
    }
    gsl_permutation_free(p);
}
/*Newton Raphson Algorithm*/
// void newtonRaphson(node * grid, const double cv,const double NX,const double NY,const double tol){double* a,double tol=1e-12
void newtonRaphson(grid& g,const double& cv,double tol=1e-12){
    // std::cout <<"newton raphson running.......\n";
    // initializeLagrangians(g);
    // double* J    = (double*)malloc(sizeof(double)* 64);
    // double* INV8 = (double*)malloc(sizeof(double)* 64);

    for (size_t i = 0; i < g.NY * g.NX; ++i){
        // int i = 42;
        // if(g.domain[i].flag == 0){
            // initializeLagrangians(g);
            // int i = 0;
            // std::cout << g.domain[5].LM[i]<< "<>"<< g.domain[i].LM[1] <<"<>"<<  g.domain[i].LM[2] <<std::endl;
            double err = 0.1;//np.array([0.1,0.1,0.1])
            int iter   = 0;
            int MAX_ITER = 3000;
            double F[8]{0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0};

            double J[64]{0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0};

            double INV8[64]{0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0};
            Matmul mul(8);
            while ( (err > tol) && (iter<MAX_ITER)  ){//  iter< 1  && (iter<MAX_ITER) 

                // std::cout <<"               "<<" Newtons ITER:  " << iter<< std::endl;//# print("ITER",iter)
                F[0] = func0(g.domain[i],cv);
                F[1] = func1(g.domain[i],cv);
                F[2] = func2(g.domain[i],cv);
                F[3] = func3(g.domain[i],cv);
                F[4] = func4(g.domain[i],cv);
                F[5] = func5(g.domain[i],cv);
                F[6] = func6(g.domain[i],cv);
                F[7] = func7(g.domain[i],cv); //(g.domain[i],cv);
                // std::cout << g.domain[i].rho<< std::endl; 

                // std::cout << F[0]<< "<>"<< F[1] << "<>"<< F[2]<<"<>"<<  F[3]<< "<>"<< F[4] << "<>"<< F[5] <<"<>"<< F[6] <<"<>"<< F[7]<< std::endl;
                // std::cout << g.domain[i].LM[0]<< "<>"<< g.domain[i].LM[1] <<"<>"<<  g.domain[i].LM[2] <<std::endl;
                //0                        1                        2                         3                           4                       5                         6                         7
                J[0] = j00(g.domain[i]);   J[1] = j01(g.domain[i]); J[2]  = j02(g.domain[i]); J[3]  = j03(g.domain[i]); J[4]  = j04(g.domain[i]); J[5]  = j05(g.domain[i]); J[6]  = j06(g.domain[i]); J[7]  = j07(g.domain[i]);
                
                J[8] = j10(g.domain[i]);   J[9] = j11(g.domain[i]); J[10] = j12(g.domain[i]); J[11] = j13(g.domain[i]); J[12] = j14(g.domain[i]); J[13] = j15(g.domain[i]); J[14] = j16(g.domain[i]); J[15] = j17(g.domain[i]);
                
                J[16] = j20(g.domain[i]); J[17] = j21(g.domain[i]); J[18] = j22(g.domain[i]); J[19] = j23(g.domain[i]); J[20] = j24(g.domain[i]); J[21] = j25(g.domain[i]); J[22] = j26(g.domain[i]); J[23] = j27(g.domain[i]);
                
                J[24] = j30(g.domain[i]); J[25] = j31(g.domain[i]); J[26] = j32(g.domain[i]); J[27] = j33(g.domain[i]); J[28] = j34(g.domain[i]); J[29] = j35(g.domain[i]); J[30] = j36(g.domain[i]); J[31] = j37(g.domain[i]);
                
                J[32] = j40(g.domain[i]); J[33] = j41(g.domain[i]); J[34] = j42(g.domain[i]); J[35] = j43(g.domain[i]); J[36] = j44(g.domain[i]); J[37] = j45(g.domain[i]); J[38] = j46(g.domain[i]); J[39] = j47(g.domain[i]);
                
                J[40] = j50(g.domain[i]); J[41] = j51(g.domain[i]); J[42] = j52(g.domain[i]); J[43] = j53(g.domain[i]); J[44] = j54(g.domain[i]); J[45] = j55(g.domain[i]); J[46] = j56(g.domain[i]); J[47] = j57(g.domain[i]);
                
                J[48] = j60(g.domain[i]); J[49] = j61(g.domain[i]); J[50] = j62(g.domain[i]); J[51] = j63(g.domain[i]); J[52] = j64(g.domain[i]); J[53] = j65(g.domain[i]); J[54] = j66(g.domain[i]); J[55] = j67(g.domain[i]);
                
                J[56] = j70(g.domain[i]); J[57] = j71(g.domain[i]); J[58] = j72(g.domain[i]); J[59] = j73(g.domain[i]); J[60] = j74(g.domain[i]); J[61] = j75(g.domain[i]); J[62] = j76(g.domain[i]); J[63] = j77(g.domain[i]);

                inv8(8,INV8,J);

                mul.Mul(INV8,F);

                F[0] = g.domain[i].LM[0] - mul.res[0];
                F[1] = g.domain[i].LM[1] - mul.res[1];
                F[2] = g.domain[i].LM[2] - mul.res[2];
                F[3] = g.domain[i].LM[3] - mul.res[3];
                F[4] = g.domain[i].LM[4] - mul.res[4];
                F[5] = g.domain[i].LM[5] - mul.res[5];
                F[6] = g.domain[i].LM[6] - mul.res[6];
                F[7] = g.domain[i].LM[7] - mul.res[7];

                norm(err,F,g.domain[i].LM);      // err = abs(new_estimate - init_guess);
                // std::cout << err << std::endl;
                // norm(err,F,a);
                // a[0] = F[0];
                // a[1] = F[1];
                // a[2] = F[2];
                g.domain[i].LM[0] = F[0];
                g.domain[i].LM[1] = F[1];
                g.domain[i].LM[2] = F[2];
                g.domain[i].LM[3] = F[3];
                g.domain[i].LM[4] = F[4];
                g.domain[i].LM[5] = F[5];
                g.domain[i].LM[6] = F[6];
                g.domain[i].LM[7] = F[7];
                // std::cout << g.domain[i].LM[0] << " " <<g.domain[i].LM[1]<< std::endl; 
                iter+=1;
                // if(iter == MAX_ITER){std::cout << ">MAX iter reached<"<<std::endl;}
                // }
            }
        // g.domain[i].LM[0] = F[0];
        // g.domain[i].LM[1] = F[1];
        // g.domain[i].LM[2] = F[2];
        // g.domain[i].LM[3] = F[3];
        // g.domain[i].LM[4] = F[4];
        // g.domain[i].LM[5] = F[5];
        // g.domain[i].LM[6] = F[6];
        // g.domain[i].LM[7] = F[7];
                // std::cout << g.domain[i].LM[0] << " <>" <<g.domain[i].LM[1] <<"  <>" <<g.domain[i].LM[2] << std::endl; 
        // }
    }

    // free(J);
    // free(INV8);

}


// /*newton rapson boundary we can use tis one to calculate equilibrium populations at the boundary*/
// void newtonRaphsonOneNode(grid& g,const double& cv, int i, double tol=1e-14){
//     // std::cout <<"newton raphson running.......\n";
//     // initializeLagrangians(g);
    
//     // for (size_t i = 0; i < g.NY * g.NX; ++i){
//         if(g.domain[i].flag == 0){
//             // initializeLagrangians(g);
//             //g.domain[i].rho
//             // int i = 0;
//             // std::cout << g.domain[5].LM[i]<< "<>"<< g.domain[i].LM[1] <<"<>"<<  g.domain[i].LM[2] <<std::endl;
//             double err = 0.1;//np.array([0.1,0.1,0.1])
//             int iter   = 0;
//             int MAX_ITER = 5000;
//             double F[3]{0.0, 0.0, 0.0};
//             double J[9]{0.0, 0.0, 0.0,
//                         0.0, 0.0, 0.0,
//                         0.0, 0.0, 0.0};
//             Matmul mul(3);
//             Inverse D(9); //create Inverse 3*3 = 9 array
//             while ( (err > tol)  && (iter<MAX_ITER) ){//

//                 // std::cout <<"               "<<" Newtons ITER:  " << iter<< std::endl;//# print("ITER",iter)
//                 F[0] = func1(g.domain[i],cv);
//                 F[1] = func2(g.domain[i],cv);
//                 F[2] = func3(g.domain[i],cv); //(g.domain[i],cv); 
//                 // std::cout << F[0]<< "<>"<< F[1] << "<>"<< F[2]<<std::endl;
//                 // std::cout << g.domain[i].LM[0]<< "<>"<< g.domain[i].LM[1] <<"<>"<<  g.domain[i].LM[2] <<std::endl;
//                 J[0] = j11(g.domain[i]);
//                 J[1] = j12(g.domain[i]);
//                 J[2] = j13(g.domain[i]);
//                 J[3] = j21(g.domain[i]);
//                 J[4] = j22(g.domain[i]);
//                 J[5] = j23(g.domain[i]);
//                 J[6] = j31(g.domain[i]);
//                 J[7] = j32(g.domain[i]);
//                 J[8] = j33(g.domain[i]);

//                 D.assign(J); //assign jacobian to inverse object
//                 // std::cout<< J[5] << " "<< J[6] << " "<< J[7] << " "<< J[8] << " "<< std::endl;
//                 D.Det();     //get the determinant
//                 // std::cout << "det: " << D.DET << std::endl;
//                 if(isnan(D.DET)){
//                     std::cout << ">>>>>>>>>>>>>>Error: SINGULAR MATRIC" << std::endl;
//                 }
//                 // if(abs(D.DET) >= 1e-10){
//                 // std::cout <<"Lagrangians initialized with previous values" <<std::endl;
//                 D.Inv();
//                 // std::cout << D << std::endl;
//                 // std::cout<< D.INV[0] << " "<< D.INV[1] << " "<< D.INV[2] << " "<< D.INV[3] << " "<< std::endl;
//                 mul.Mul(D.INV,F);
//                 F[0] = g.domain[i].LM[0] - mul.res[0];
//                 F[1] = g.domain[i].LM[1] - mul.res[1];
//                 F[2] = g.domain[i].LM[2] - mul.res[2];
//                 norm(err,F,g.domain[i].LM);      // err = abs(new_estimate - init_guess);
//                 // std::cout << err << std::endl;
//                 // norm(err,F,a);
//                 // a[0] = F[0];
//                 // a[1] = F[1];
//                 // a[2] = F[2];
//                 g.domain[i].LM[0] = F[0];
//                 g.domain[i].LM[1] = F[1];
//                 g.domain[i].LM[2] = F[2];
//                 // std::cout << g.domain[i].LM[0] << " " <<g.domain[i].LM[1]<< std::endl; 
//                 iter+=1;
//                 // if(iter == MAX_ITER){std::cout << ">MAX iter reached<"<<std::endl;}
//                 // }
//             }
//                 // std::cout << g.domain[i].LM[0] << " <>" <<g.domain[i].LM[1] <<"  <>" <<g.domain[i].LM[2] << std::endl; 
//         }

//     // }
    

// }

/*Newton Raphson 2DAlgorithm*/
// void newtonRaphson(double* array0, double tol = 1e-12){
        
//         double err = 0.1;//np.array([0.1,0.1,0.1])
//         size_t iter = 0;
//         double F[2]{0.0,0.0};
//         double J[4]{0.0,0.0,
//                     0.0,0.0};
//         Matmul mul(2);
//         Inverse  D(4); //create Inverse object
//         while ( (err > tol) ){
//             std::cout << "iter = "<<iter<<std::endl;
//             F[0] = func1(array0[0],array0[1]);
//             F[1] = func2(array0[0],array0[1]);

//             J[0] = j11(array0[0],array0[1]);
//             J[1] = j12(array0[0],array0[1]);
//             J[2] = j21(array0[0],array0[1]);
//             J[3] = j22(array0[0],array0[1]);

//             D.assign(J); //assign jacobian to inverse object
//             D.Det();     //get the determinant
//             if(abs(D.DET) >= 1e-10){
//                 D.Inv();
//                 mul.Mul(D.INV,F);
//                 F[0] = array0[0] - mul.res[0];
//                 F[1] = array0[1] - mul.res[1];
//                 norm(err,F,array0);
//                 array0[0] = F[0];
//                 array0[1] = F[1];
//                 // 
//                 iter+=1;}
//         }

//         // return 0;
// }



// double func1(double x,double y, double z){

//     double f1  = x*x + y*y + exp(z*z) -3.0;
//     return f1;
// }

// double func2(double x,double y, double z){
//     double f2  = x*x + y*y - z - 1.0;
//     return f2;
// }

// double func3(double x,double y, double z){
//     double f3  = exp(x) + y - z - 3.0; 
//     return f3;
// }

// double j11(double x,double y, double z){
//     double j1 = 2.0*x;
//     return j1;

// }

// double j12(double x,double y, double z){
//     double j2  = 2.0*y;
//     return j2;
// }

// double j13(double x,double y, double z){
//     double j3  = 2.0*z*exp(z*z);
//     return j3;
// }


// double j21(double x,double y, double z){
//     double j4  = 2.0*x;
//     return j4;
// }

// double j22(double x,double y, double z){
//     double j5  = 2.0*y;
//     return j5;
// }

// double j23(double x,double y, double z){
//     double j6  = -1.0;
//     return j6;
// }

// double j31(double x,double y, double z){
//     double j7  = exp(x);
//     return j7;
// }

// double j32(double x,double y, double z){
//     double j8  = 1.0;
//     return j8;
// }

// double j33(double x,double y, double z){
//     double j9  = -1.0;
//     return j9;
// }
// // double* constraints(double* F,const node& p, double cv){
// //     F[0] = func1(p,cv);
// //     F[1] = func2(p,cv);
// //     F[2] = func3(p,cv);
// //     return F;
// // }

// // double* Jacobian(double* J,const node& p){
// //     J[0]  =  j11(p);
// //     J[1]  =  j12(p);
// //     J[2]  =  j13(p);
// //     J[3]  =  j21(p);
// //     J[4]  =  j22(p);
// //     J[5]  =  j23(p);
// //     J[6]  =  j31(p);
// //     J[7]  =  j32(p);
// //     J[8]  =  j33(p);
// //     return J;
// // }




// double randZeroToOne()
// {
//     return (double)rand() / (double)RAND_MAX;//rand() / (RAND_MAX + 1.);
// }