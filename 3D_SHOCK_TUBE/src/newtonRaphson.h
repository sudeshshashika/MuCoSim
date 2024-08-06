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
#include "state.h"
#include <gsl/gsl_linalg.h>

/*Residual function definition*/
/*constraint1: Total energy*/
double func1(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                              
    // x*x + y*y + exp(z*z) -3.0;
    // std::cout << p.cx[1]<< std::endl;
    double f1Sum = 0.0;
    for(int q=0; q<27; ++q){
        f1Sum = f1Sum + p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                        + x2*p.cx[q] 
                                                        + x3*p.cy[q]
                                                        + x4*p.cz[q]) );
    }
    double f1  = f1Sum - 2.0 * p.rho  * (p.T * cv + (p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) * 0.5);
    return f1;
}
/*constraint2: ux: trace of Maxwell boltzmann moment*/
double func2(node& p, double cv){
    /*Lagrangian Multipliers*/
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2];
    double x4    = p.LM[3];                                                
    // std::cout << p.cx[5] << std::endl;
    // double f2 = x1 * x1 - x2 +1.0;
    // std::cout << p.cx[1]<< std::endl;
    double f2Sum = 0.0;
    for(int q=0; q<27; ++q){
        f2Sum = f2Sum + p.cx[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
    double f2  = f2Sum - 2.0  * p.rho * p.ux   * ( p.T * cv + (p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) * 0.5 + p.T );
    return f2;
}
/*constraint3: uy: trace of Maxwell boltzmann moment*/
double func3(const node& p, const double cv){
    /*Lagrangian Multipliers*/
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f3Sum = 0.0;
    for(int q=0; q<27; ++q){
        f3Sum = f3Sum + p.cy[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
    double f3  =  f3Sum - 2.0  * p.rho * p.uy   * (p.T * cv + (p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) * 0.5 + p.T);
    return f3;
}
/*constraint4: uz: trace of Maxwell boltzmann moment*/
double func4(const node& p, const double cv){
    /*Lagrangian Multipliers*/
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f4Sum = 0.0;
    for(int q=0; q<27; ++q){
        f4Sum = f4Sum + p.cz[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
    double f4  =  f4Sum - 2.0  * p.rho * p.uz   * (p.T * cv + (p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) * 0.5 + p.T);
    return f4;
}
/*Jacobian Definitions*/
//first row
double j11(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2];
    double x4    = p.LM[3]; 
    double f1Sum = 0.0;
    for(int q=0; q<27; ++q){
        f1Sum = f1Sum + p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                        + x2*p.cx[q] 
                                                        + x3*p.cy[q]
                                                        + x4*p.cz[q]) );
    }
   return f1Sum;
}
double j12(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3]; 
    double f1Sum = 0.0;
    for(int q=0; q<27; ++q){
        f1Sum = f1Sum + p.cx[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
   return f1Sum;
}
double j13(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3]; 
    double f1Sum = 0.0;
    for(int q=0; q<27; ++q){
        f1Sum = f1Sum + p.cy[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
   return f1Sum;
}
double j14(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3]; 
    double f1Sum = 0.0;
    for(int q=0; q<27; ++q){
        f1Sum = f1Sum + p.cz[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
   return f1Sum;
}
//second row
double j21(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                                
    double f2Sum = 0.0;
    for(int q=0; q<27; ++q){
        f2Sum = f2Sum + p.cx[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
    return f2Sum;
}
double j22(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                                
    double f2Sum = 0.0;
    for(int q=0; q<27; ++q){
        f2Sum = f2Sum + p.cx[q] * p.cx[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f2Sum;
}
double j23(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                                
    double f2Sum = 0.0;
    for(int q=0; q<27; ++q){
        f2Sum = f2Sum + p.cy[q] * p.cx[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f2Sum;
}
double j24(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                                
    double f2Sum = 0.0;
    for(int q=0; q<27; ++q){
        f2Sum = f2Sum + p.cz[q] * p.cx[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f2Sum;
}
// third row
double j31(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f3Sum = 0.0;
    for(int q=0; q<27; ++q){
        f3Sum = f3Sum + p.cy[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
    return f3Sum;
}
double j32(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f3Sum = 0.0;
    for(int q=0; q<27; ++q){
        f3Sum = f3Sum + p.cx[q] * p.cy[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f3Sum;
}
double j33(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f3Sum = 0.0;
    for(int q=0; q<27; ++q){
        f3Sum = f3Sum + p.cy[q] * p.cy[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f3Sum;
}
double j34(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f3Sum = 0.0;
    for(int q=0; q<27; ++q){
        f3Sum = f3Sum + p.cz[q] * p.cy[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f3Sum;
}
// fourth row
double j41(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f4Sum = 0.0;
    for(int q=0; q<27; ++q){
        f4Sum = f4Sum + p.cz[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                  + x2*p.cx[q] 
                                                                  + x3*p.cy[q]
                                                                  + x4*p.cz[q]) );
    }
    return f4Sum;
}
double j42(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f4Sum = 0.0;
    for(int q=0; q<27; ++q){
        f4Sum = f4Sum + p.cx[q] * p.cz[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f4Sum;
}
double j43(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f4Sum = 0.0;
    for(int q=0; q<27; ++q){
        f4Sum = f4Sum + p.cy[q] * p.cz[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f4Sum;
}
double j44(node& p){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2]; 
    double x4    = p.LM[3];                                             
    double f4Sum = 0.0;
    for(int q=0; q<27; ++q){
        f4Sum = f4Sum + p.cz[q] * p.cz[q] * p.rho * p.w[q] * exp( 1.0* (0.0 + x1 
                                                                            + x2*p.cx[q] 
                                                                            + x3*p.cy[q]
                                                                            + x4*p.cz[q]) );
    }
    return f4Sum;
}
// 

// 
void norm(double& err, double* a1, double* a2){
    /*3d*/
    err  = (a1[0] - a2[0]) * (a1[0] - a2[0]) +\
           (a1[1] - a2[1]) * (a1[1] - a2[1]) +\
           (a1[2] - a2[2]) * (a1[2] - a2[2]) +\
           (a1[3] - a2[3]) * (a1[3] - a2[3]) ;

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
// 
void allfuncs(const node& p, const double cv, double* tmp){
    double x1    = p.LM[0];
    double x2    = p.LM[1];
    double x3    = p.LM[2];  
    double x4    = p.LM[3];                                            
    
    double j0  = 0.0, j1  = 0.0, j2  = 0.0, j3  = 0.0, j4  = 0.0, j5  = 0.0, j6  = 0.0, j7  = 0.0, j8  = 0.0;
    double j9  = 0.0, j10  = 0.0, j11  = 0.0, j12  = 0.0, j13  = 0.0, j14  = 0.0, j15  = 0.0;  

    double f1  = (-2.0  * p.rho         * (p.T * cv + (p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) * 0.5));
    double f2  = (-2.0  * p.rho * p.ux  * (p.T * cv + (p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) * 0.5 + p.T));
    double f3  = (-2.0  * p.rho * p.uy  * (p.T * cv + (p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) * 0.5 + p.T));
    double f4  = (-2.0  * p.rho * p.uz  * (p.T * cv + (p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) * 0.5 + p.T));


    for(int i = 0; i < 27; ++i){ 
        double tmp = p.rho * p.w[i] * exp(1.0*(0.0  + x1 
                                                    + x2*p.cx[i] 
                                                    + x3*p.cy[i]
                                                    + x4*p.cz[i]));
        f1 += tmp;
        f2 += p.cx[i] * tmp;
        f3 += p.cy[i] * tmp;
        f4 += p.cz[i] * tmp;

        j0 += 1.0 * tmp;
        j1 += 1.0 * p.cx[i] * tmp;
        j2 += 1.0 * p.cy[i] * tmp;
        j3 += 1.0 * p.cz[i] * tmp;

        j4 += 1.0 * p.cx[i] * tmp;
        j5 += 1.0 * p.cx[i] * p.cx[i] * tmp;
        j6 += 1.0 * p.cx[i] * p.cy[i] * tmp;
        j7 += 1.0 * p.cx[i] * p.cz[i] * tmp;

        j8  += 1.0 * p.cy[i] * tmp;
        j9  += 1.0 * p.cy[i] * p.cx[i] * tmp;
        j10 += 1.0 * p.cy[i] * p.cy[i] * tmp;
        j11 += 1.0 * p.cy[i] * p.cz[i] * tmp;

        j12  += 1.0 * p.cz[i] * tmp;
        j13  += 1.0 * p.cz[i] * p.cx[i] * tmp;
        j14  += 1.0 * p.cz[i] * p.cy[i] * tmp;
        j15  += 1.0 * p.cz[i] * p.cz[i] * tmp;
    }

    tmp[0] = f1;
    tmp[1] = f2;
    tmp[2] = f3;
    tmp[3] = f4;

    tmp[4] = j0;
    tmp[5] = j1;
    tmp[6] = j2;
    tmp[7] = j3;

    tmp[8]  = j4;
    tmp[9]  = j5;
    tmp[10] = j6;
    tmp[11] = j7;

    tmp[12] = j8;
    tmp[13] = j9;
    tmp[14] = j10;
    tmp[15] = j11;

    tmp[16] = j12;
    tmp[17] = j13;
    tmp[18] = j14;
    tmp[19] = j15;
}
// 
void newtonRaphson(grid& g,const double& cv,int flagIndex,double tol=1e-12){
    #pragma omp parallel for
    for (int i = 0; i < g.NZ * g.NY * g.NX; ++i){
        if(g.domain[i].flag == flagIndex){
            double err   = 0.1;//np.array([0.1,0.1,0.1])
            int iter     = 0;
            int MAX_ITER = 5000;
            double F[4]    { 0.0, 0.0, 0.0, 0.0};
            double J[16]   { 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0};
            double INV8[16]{ 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0};
            Matmul mul(4);
            // Inverse D(16); //create Inverse 3*3 = 9 array
            while ( (err > tol)  && (iter<MAX_ITER) ){//

               
                double tmp_res[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
                allfuncs(g.domain[i],cv, tmp_res);
                F[0] = tmp_res[0];
                F[1] = tmp_res[1];
                F[2] = tmp_res[2];
                F[3] = tmp_res[3];

                J[0]  = tmp_res[4];
                J[1]  = tmp_res[5];
                J[2]  = tmp_res[6];
                J[3]  = tmp_res[7];
                J[4]  = tmp_res[8];
                J[5]  = tmp_res[9];
                J[6]  = tmp_res[10];
                J[7]  = tmp_res[11];
                J[8]  = tmp_res[12]; 
                J[9]  = tmp_res[13]; 
                J[10] = tmp_res[14]; 
                J[11] = tmp_res[15]; 
                J[12] = tmp_res[16]; 
                J[13] = tmp_res[17]; 
                J[14] = tmp_res[18];
                J[15] = tmp_res[19]; 
 
                inv8(4,INV8,J);

                mul.Mul(INV8,F);
                F[0] = g.domain[i].LM[0] - mul.res[0];
                F[1] = g.domain[i].LM[1] - mul.res[1];
                F[2] = g.domain[i].LM[2] - mul.res[2];
                F[3] = g.domain[i].LM[3] - mul.res[3];

                norm(err,F,g.domain[i].LM);      // err = abs(new_estimate - init_guess);

                g.domain[i].LM[0] = F[0];
                g.domain[i].LM[1] = F[1];
                g.domain[i].LM[2] = F[2];
                g.domain[i].LM[3] = F[3];

                iter+=1;

            }
            g.domain[i].nr = iter;
        }
        
    }

}
