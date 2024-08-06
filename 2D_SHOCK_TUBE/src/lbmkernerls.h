#pragma once
#include <cmath> 
#include <omp.h>
#include "math.h"
#include "node.h"
#include "grid.h"
#include "auxliray.h"
#include "state.h"
//#define CHUNK_SIZE 1

void initializeLagrangians(grid& g){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int i = 0; i < g.NX * g.NY; ++i){
            g.domain[i].LM[0] = 0.5;
            g.domain[i].LM[1] = 0.5;
            g.domain[i].LM[2] = 0.5;
    }
}

void initializationMacro(grid& g,state d,int flagIndex){
    
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY; ++j){
        for(int i =0; i < g.NX; ++i){
            if(g(j,i).flag == flagIndex){
                g(j,i).ux  = d.u_x; 
                g(j,i).uy  = d.u_y;
                if(i <=  g.NX/2){
                    g(j,i).rho = d.rho_l;
                    g(j,i).T   = d.T_l;
                }
                if(i > g.NX/2){
                    g(j,i).rho = d.rho_r;
                    g(j,i).T   = d.T_r;
                }
            }
        }
    }
}

//Initialize f populations
void initilizeDDF(grid& g){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY * g.NX; ++j){
        for(int i = 0; i < 9; ++i){
            // if(g.domain[j].flag ==0){
                g.domain[j].fold[i]  = g.domain[j].feq[i];
                g.domain[j].fnew[i]  = g.domain[j].feq[i];
                g.domain[j].gold[i]  = g.domain[j].geq[i];            
                g.domain[j].gnew[i]  = g.domain[j].geq[i];
            // }
        }
    }
}
// Set temperature dependent weights
void setWeights(grid& g,int flagIndex=0){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY * g.NX; ++j ){
        // if(g.domain[j].flag == flagIndex){ //flag zeros is the fluid cells
            g.domain[j].w[0]= (1.0 - g.domain[j].T) * (1.0 - g.domain[j].T);  //# 0 *  0

            g.domain[j].w[1]= (0.5 * g.domain[j].T) * (1.0 - g.domain[j].T);  //# 1 *  0 
            g.domain[j].w[2]= (0.5 * g.domain[j].T) * (1.0 - g.domain[j].T);  //# 0 *  1
            g.domain[j].w[3]= (0.5 * g.domain[j].T) * (1.0 - g.domain[j].T);  //#-1 *  0 
            g.domain[j].w[4]= (0.5 * g.domain[j].T) * (1.0 - g.domain[j].T);  //# 0 * -1

            g.domain[j].w[5]= (0.5 * g.domain[j].T) * (0.5 * g.domain[j].T);  //# 1 *  1
            g.domain[j].w[6]= (0.5 * g.domain[j].T) * (0.5 * g.domain[j].T);  //#-1 *  1
            g.domain[j].w[7]= (0.5 * g.domain[j].T) * (0.5 * g.domain[j].T);  //#-1 * -1
            g.domain[j].w[8]= (0.5 * g.domain[j].T) * (0.5 * g.domain[j].T);  //# 1 * -1
        // }
    }
}
// Calculate the macroscopic parameters
void calc_parameters_from_PDF(grid& g,const double& cv, int flagIndex){
    #pragma omp parallel
    {
        #pragma omp for //schedule(static,CHUNK_SIZE)
        for (int j = 0; j < g.NY * g.NX; ++j){
            if(g.domain[j].flag == flagIndex){
                //Density
                g.domain[j].rho = g.domain[j].fold[0] +\
                                g.domain[j].fold[1] +\
                                g.domain[j].fold[2] +\
                                g.domain[j].fold[3] +\
                                g.domain[j].fold[4] +\
                                g.domain[j].fold[5] +\
                                g.domain[j].fold[6] +\
                                g.domain[j].fold[7] +\
                                g.domain[j].fold[8] ;
                if(isnan(g.domain[j].rho)){std::cout << "-nan- encountrerd in rho" << std::endl;}
                //Velocity C
                g.domain[j].ux =  ( g.domain[j].fold[0] * g.domain[j].cx[0]  +\
                                    g.domain[j].fold[1] * g.domain[j].cx[1]  +\
                                    g.domain[j].fold[2] * g.domain[j].cx[2]  +\
                                    g.domain[j].fold[3] * g.domain[j].cx[3]  +\
                                    g.domain[j].fold[4] * g.domain[j].cx[4]  +\
                                    g.domain[j].fold[5] * g.domain[j].cx[5]  +\
                                    g.domain[j].fold[6] * g.domain[j].cx[6]  +\
                                    g.domain[j].fold[7] * g.domain[j].cx[7]  +\
                                    g.domain[j].fold[8] * g.domain[j].cx[8]  )/g.domain[j].rho;
                if(isnan(g.domain[j].ux)){std::cout << "-nan- encountrerd in ux" << std::endl;}
                //Velocity Y
                g.domain[j].uy =  ( g.domain[j].fold[0] * g.domain[j].cy[0]  +\
                                    g.domain[j].fold[1] * g.domain[j].cy[1]  +\
                                    g.domain[j].fold[2] * g.domain[j].cy[2]  +\
                                    g.domain[j].fold[3] * g.domain[j].cy[3]  +\
                                    g.domain[j].fold[4] * g.domain[j].cy[4]  +\
                                    g.domain[j].fold[5] * g.domain[j].cy[5]  +\
                                    g.domain[j].fold[6] * g.domain[j].cy[6]  +\
                                    g.domain[j].fold[7] * g.domain[j].cy[7]  +\
                                    g.domain[j].fold[8] * g.domain[j].cy[8]  )/g.domain[j].rho;
                if(isnan(g.domain[j].ux)){std::cout << "-nan- encountrerd in uy" << std::endl;}
                //Energy
                double gg     =   ( g.domain[j].gold[0] +\
                                    g.domain[j].gold[1] +\
                                    g.domain[j].gold[2] +\
                                    g.domain[j].gold[3] +\
                                    g.domain[j].gold[4] +\
                                    g.domain[j].gold[5] +\
                                    g.domain[j].gold[6] +\
                                    g.domain[j].gold[7] +\
                                    g.domain[j].gold[8] )/ (2.0*g.domain[j].rho);
                // Velocity magnitude
                double velMag  = g.domain[j].ux * g.domain[j].ux + g.domain[j].uy * g.domain[j].uy;
                //Temperature
                g.domain[j].T   = ( gg - (velMag * 0.5) ) / cv; 
                if(isnan(g.domain[j].T)){std::cout << "-nan- encountrerd in T" << std::endl;}
            }
        }
    }
}
// Calculate F equilibrium population
void calcFeq(grid& g,int flagIndex){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    
    for(int j = 0; j < g.NY * g.NX; ++j){
        if(g.domain[j].flag == flagIndex){
            g.domain[j].feq[0] = 1.0 * g.domain[j].rho * ( 1.0 - (g.domain[j].ux - g.domain[j].ushift) * (g.domain[j].ux - g.domain[j].ushift) - g.domain[j].T)                                           * ( 1.0 - g.domain[j].uy * g.domain[j].uy - g.domain[j].T);                     // 0 * 0

            g.domain[j].feq[1] = 0.5 * g.domain[j].rho * ( 1.0 * (g.domain[j].ux - g.domain[j].ushift) + (g.domain[j].ux - g.domain[j].ushift) * (g.domain[j].ux - g.domain[j].ushift) + g.domain[j].T)   * ( 1.0 - g.domain[j].uy * g.domain[j].uy - g.domain[j].T);                     // 1 *  0
            g.domain[j].feq[2] = 0.5 * g.domain[j].rho * ( 1.0 - (g.domain[j].ux - g.domain[j].ushift) * (g.domain[j].ux - g.domain[j].ushift) - g.domain[j].T)                                           * ( 1.0 * g.domain[j].uy + g.domain[j].uy * g.domain[j].uy + g.domain[j].T);    // 0 *  1       
            g.domain[j].feq[3] = 0.5 * g.domain[j].rho * (-1.0 * (g.domain[j].ux - g.domain[j].ushift) + (g.domain[j].ux - g.domain[j].ushift) * (g.domain[j].ux - g.domain[j].ushift) + g.domain[j].T)   * ( 1.0 - g.domain[j].uy * g.domain[j].uy - g.domain[j].T);                     //-1 *  0
            g.domain[j].feq[4] = 0.5 * g.domain[j].rho * ( 1.0 - (g.domain[j].ux - g.domain[j].ushift) * (g.domain[j].ux - g.domain[j].ushift) - g.domain[j].T)                                           * (-1.0 * g.domain[j].uy + g.domain[j].uy * g.domain[j].uy + g.domain[j].T);    // 0 * -1
            
            g.domain[j].feq[5] = 0.25* g.domain[j].rho * ( 1.0 * (g.domain[j].ux - g.domain[j].ushift) + (g.domain[j].ux - g.domain[j].ushift) *   (g.domain[j].ux - g.domain[j].ushift) + g.domain[j].T) * ( 1.0 * g.domain[j].uy + g.domain[j].uy * g.domain[j].uy + g.domain[j].T);  // 1 *  1
            g.domain[j].feq[6] = 0.25* g.domain[j].rho * (-1.0 * (g.domain[j].ux - g.domain[j].ushift) + (g.domain[j].ux - g.domain[j].ushift) *   (g.domain[j].ux - g.domain[j].ushift) + g.domain[j].T) * ( 1.0 * g.domain[j].uy + g.domain[j].uy * g.domain[j].uy + g.domain[j].T);  //-1 *  1
            g.domain[j].feq[7] = 0.25* g.domain[j].rho * (-1.0 * (g.domain[j].ux - g.domain[j].ushift) + (g.domain[j].ux - g.domain[j].ushift) *   (g.domain[j].ux - g.domain[j].ushift) + g.domain[j].T) * (-1.0 * g.domain[j].uy + g.domain[j].uy * g.domain[j].uy + g.domain[j].T);  //-1 * -1
            g.domain[j].feq[8] = 0.25* g.domain[j].rho * ( 1.0 * (g.domain[j].ux - g.domain[j].ushift) + (g.domain[j].ux - g.domain[j].ushift) *   (g.domain[j].ux - g.domain[j].ushift) + g.domain[j].T) * (-1.0 * g.domain[j].uy + g.domain[j].uy * g.domain[j].uy + g.domain[j].T);  // 1 * -1
        }
    }
}
/*calculate G equilibrium population*/
void calcGeq(grid& g,int flagIndex){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY * g.NX; ++j){
        if(g.domain[j].flag == flagIndex){
        g.domain[j].geq[0] = g.domain[j].rho * g.domain[j].w[0] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[0] + g.domain[j].LM[2] * g.domain[j].cy[0])); 
        g.domain[j].geq[1] = g.domain[j].rho * g.domain[j].w[1] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[1] + g.domain[j].LM[2] * g.domain[j].cy[1])); 
        g.domain[j].geq[2] = g.domain[j].rho * g.domain[j].w[2] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[2] + g.domain[j].LM[2] * g.domain[j].cy[2]));
        g.domain[j].geq[3] = g.domain[j].rho * g.domain[j].w[3] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[3] + g.domain[j].LM[2] * g.domain[j].cy[3])); 
        g.domain[j].geq[4] = g.domain[j].rho * g.domain[j].w[4] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[4] + g.domain[j].LM[2] * g.domain[j].cy[4])); 
        g.domain[j].geq[5] = g.domain[j].rho * g.domain[j].w[5] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[5] + g.domain[j].LM[2] * g.domain[j].cy[5]));
        g.domain[j].geq[6] = g.domain[j].rho * g.domain[j].w[6] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[6] + g.domain[j].LM[2] * g.domain[j].cy[6])); 
        g.domain[j].geq[7] = g.domain[j].rho * g.domain[j].w[7] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[7] + g.domain[j].LM[2] * g.domain[j].cy[7])); 
        g.domain[j].geq[8] = g.domain[j].rho * g.domain[j].w[8] * exp(1.0*(g.domain[j].LM[0] + g.domain[j].LM[1] * g.domain[j].cx[8] + g.domain[j].LM[2] * g.domain[j].cy[8]));
        }
    }
}
/*calculate Pressure tensor*/
void pTensor(grid& g){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY * g.NX; ++j){
        if(g.domain[j].flag == 0){
        //xx
            g.domain[j].P[0] =  g.domain[j].cx[0] * g.domain[j].cx[0] * g.domain[j].fold[0] +\
                                g.domain[j].cx[1] * g.domain[j].cx[1] * g.domain[j].fold[1] +\
                                g.domain[j].cx[2] * g.domain[j].cx[2] * g.domain[j].fold[2] +\
                                g.domain[j].cx[3] * g.domain[j].cx[3] * g.domain[j].fold[3] +\
                                g.domain[j].cx[4] * g.domain[j].cx[4] * g.domain[j].fold[4] +\
                                g.domain[j].cx[5] * g.domain[j].cx[5] * g.domain[j].fold[5] +\
                                g.domain[j].cx[6] * g.domain[j].cx[6] * g.domain[j].fold[6] +\
                                g.domain[j].cx[7] * g.domain[j].cx[7] * g.domain[j].fold[7] +\
                                g.domain[j].cx[8] * g.domain[j].cx[8] * g.domain[j].fold[8] ;
        //xy                    
            g.domain[j].P[1] =  g.domain[j].cx[0] * g.domain[j].cy[0] * g.domain[j].fold[0] +\
                                g.domain[j].cx[1] * g.domain[j].cy[1] * g.domain[j].fold[1] +\
                                g.domain[j].cx[2] * g.domain[j].cy[2] * g.domain[j].fold[2] +\
                                g.domain[j].cx[3] * g.domain[j].cy[3] * g.domain[j].fold[3] +\
                                g.domain[j].cx[4] * g.domain[j].cy[4] * g.domain[j].fold[4] +\
                                g.domain[j].cx[5] * g.domain[j].cy[5] * g.domain[j].fold[5] +\
                                g.domain[j].cx[6] * g.domain[j].cy[6] * g.domain[j].fold[6] +\
                                g.domain[j].cx[7] * g.domain[j].cy[7] * g.domain[j].fold[7] +\
                                g.domain[j].cx[8] * g.domain[j].cy[8] * g.domain[j].fold[8] ;
        //yx
            g.domain[j].P[2] =  g.domain[j].cx[0] * g.domain[j].cy[0] * g.domain[j].fold[0] +\
                                g.domain[j].cx[1] * g.domain[j].cy[1] * g.domain[j].fold[1] +\
                                g.domain[j].cx[2] * g.domain[j].cy[2] * g.domain[j].fold[2] +\
                                g.domain[j].cx[3] * g.domain[j].cy[3] * g.domain[j].fold[3] +\
                                g.domain[j].cx[4] * g.domain[j].cy[4] * g.domain[j].fold[4] +\
                                g.domain[j].cx[5] * g.domain[j].cy[5] * g.domain[j].fold[5] +\
                                g.domain[j].cx[6] * g.domain[j].cy[6] * g.domain[j].fold[6] +\
                                g.domain[j].cx[7] * g.domain[j].cy[7] * g.domain[j].fold[7] +\
                                g.domain[j].cx[8] * g.domain[j].cy[8] * g.domain[j].fold[8] ;
        //yy
            g.domain[j].P[3] =  g.domain[j].cy[0] * g.domain[j].cy[0] * g.domain[j].fold[0] +\
                                g.domain[j].cy[1] * g.domain[j].cy[1] * g.domain[j].fold[1] +\
                                g.domain[j].cy[2] * g.domain[j].cy[2] * g.domain[j].fold[2] +\
                                g.domain[j].cy[3] * g.domain[j].cy[3] * g.domain[j].fold[3] +\
                                g.domain[j].cy[4] * g.domain[j].cy[4] * g.domain[j].fold[4] +\
                                g.domain[j].cy[5] * g.domain[j].cy[5] * g.domain[j].fold[5] +\
                                g.domain[j].cy[6] * g.domain[j].cy[6] * g.domain[j].fold[6] +\
                                g.domain[j].cy[7] * g.domain[j].cy[7] * g.domain[j].fold[7] +\
                                g.domain[j].cy[8] * g.domain[j].cy[8] * g.domain[j].fold[8] ;
        }
    }
}
/*calculate Equilibrium pressure tensor*/
void pEqTensor(grid& g){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY * g.NX; ++j){
        if(g.domain[j].flag == 0){
        //xx
            g.domain[j].Peq[0] =  g.domain[j].cx[0] * g.domain[j].cx[0] * g.domain[j].feq[0] +\
                                g.domain[j].cx[1] * g.domain[j].cx[1] * g.domain[j].feq[1] +\
                                g.domain[j].cx[2] * g.domain[j].cx[2] * g.domain[j].feq[2] +\
                                g.domain[j].cx[3] * g.domain[j].cx[3] * g.domain[j].feq[3] +\
                                g.domain[j].cx[4] * g.domain[j].cx[4] * g.domain[j].feq[4] +\
                                g.domain[j].cx[5] * g.domain[j].cx[5] * g.domain[j].feq[5] +\
                                g.domain[j].cx[6] * g.domain[j].cx[6] * g.domain[j].feq[6] +\
                                g.domain[j].cx[7] * g.domain[j].cx[7] * g.domain[j].feq[7] +\
                                g.domain[j].cx[8] * g.domain[j].cx[8] * g.domain[j].feq[8] ;
        //xy                    
            g.domain[j].Peq[1] =  g.domain[j].cx[0] * g.domain[j].cy[0] * g.domain[j].feq[0] +\
                                g.domain[j].cx[1] * g.domain[j].cy[1] * g.domain[j].feq[1] +\
                                g.domain[j].cx[2] * g.domain[j].cy[2] * g.domain[j].feq[2] +\
                                g.domain[j].cx[3] * g.domain[j].cy[3] * g.domain[j].feq[3] +\
                                g.domain[j].cx[4] * g.domain[j].cy[4] * g.domain[j].feq[4] +\
                                g.domain[j].cx[5] * g.domain[j].cy[5] * g.domain[j].feq[5] +\
                                g.domain[j].cx[6] * g.domain[j].cy[6] * g.domain[j].feq[6] +\
                                g.domain[j].cx[7] * g.domain[j].cy[7] * g.domain[j].feq[7] +\
                                g.domain[j].cx[8] * g.domain[j].cy[8] * g.domain[j].feq[8] ;
        //yx
            g.domain[j].Peq[2] =  g.domain[j].cx[0] * g.domain[j].cy[0] * g.domain[j].feq[0] +\
                                g.domain[j].cx[1] * g.domain[j].cy[1] * g.domain[j].feq[1] +\
                                g.domain[j].cx[2] * g.domain[j].cy[2] * g.domain[j].feq[2] +\
                                g.domain[j].cx[3] * g.domain[j].cy[3] * g.domain[j].feq[3] +\
                                g.domain[j].cx[4] * g.domain[j].cy[4] * g.domain[j].feq[4] +\
                                g.domain[j].cx[5] * g.domain[j].cy[5] * g.domain[j].feq[5] +\
                                g.domain[j].cx[6] * g.domain[j].cy[6] * g.domain[j].feq[6] +\
                                g.domain[j].cx[7] * g.domain[j].cy[7] * g.domain[j].feq[7] +\
                                g.domain[j].cx[8] * g.domain[j].cy[8] * g.domain[j].feq[8] ;
        //yy
            g.domain[j].Peq[3] =  g.domain[j].cy[0] * g.domain[j].cy[0] * g.domain[j].feq[0] +\
                                g.domain[j].cy[1] * g.domain[j].cy[1] * g.domain[j].feq[1] +\
                                g.domain[j].cy[2] * g.domain[j].cy[2] * g.domain[j].feq[2] +\
                                g.domain[j].cy[3] * g.domain[j].cy[3] * g.domain[j].feq[3] +\
                                g.domain[j].cy[4] * g.domain[j].cy[4] * g.domain[j].feq[4] +\
                                g.domain[j].cy[5] * g.domain[j].cy[5] * g.domain[j].feq[5] +\
                                g.domain[j].cy[6] * g.domain[j].cy[6] * g.domain[j].feq[6] +\
                                g.domain[j].cy[7] * g.domain[j].cy[7] * g.domain[j].feq[7] +\
                                g.domain[j].cy[8] * g.domain[j].cy[8] * g.domain[j].feq[8] ;
        }
    }  
}
/*Quasi Euilibrium G*/
void calcQuasiEqG(grid& g){
    //P = P-Peq
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY * g.NX; ++j){
        if(g.domain[j].flag == 0){
            double pxx = g.domain[j].P[0] - g.domain[j].Peq[0];
            double pyx = g.domain[j].P[1] - g.domain[j].Peq[1];
            double pxy = g.domain[j].P[2] - g.domain[j].Peq[2];
            double pyy = g.domain[j].P[3] - g.domain[j].Peq[3];

            g.domain[j].qstr[0] = g.domain[j].geq[0] + ( g.domain[j].w[0] * (2.0 * g.domain[j].ux * (g.domain[j].cx[0]*pxx + g.domain[j].cy[0]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[0]*pxy + g.domain[j].cy[0]*pyy) ) ) / g.domain[j].T;

            g.domain[j].qstr[1] = g.domain[j].geq[1] + ( g.domain[j].w[1] * (2.0 * g.domain[j].ux * (g.domain[j].cx[1]*pxx + g.domain[j].cy[1]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[1]*pxy + g.domain[j].cy[1]*pyy) ) ) / g.domain[j].T;
            g.domain[j].qstr[2] = g.domain[j].geq[2] + ( g.domain[j].w[2] * (2.0 * g.domain[j].ux * (g.domain[j].cx[2]*pxx + g.domain[j].cy[2]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[2]*pxy + g.domain[j].cy[2]*pyy) ) ) / g.domain[j].T;
            g.domain[j].qstr[3] = g.domain[j].geq[3] + ( g.domain[j].w[3] * (2.0 * g.domain[j].ux * (g.domain[j].cx[3]*pxx + g.domain[j].cy[3]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[3]*pxy + g.domain[j].cy[3]*pyy) ) ) / g.domain[j].T;
            g.domain[j].qstr[4] = g.domain[j].geq[4] + ( g.domain[j].w[4] * (2.0 * g.domain[j].ux * (g.domain[j].cx[4]*pxx + g.domain[j].cy[4]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[4]*pxy + g.domain[j].cy[4]*pyy) ) ) / g.domain[j].T;

            g.domain[j].qstr[5] = g.domain[j].geq[5] + ( g.domain[j].w[5] * (2.0 * g.domain[j].ux * (g.domain[j].cx[5]*pxx + g.domain[j].cy[5]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[5]*pxy + g.domain[j].cy[5]*pyy) ) ) / g.domain[j].T;
            g.domain[j].qstr[6] = g.domain[j].geq[6] + ( g.domain[j].w[6] * (2.0 * g.domain[j].ux * (g.domain[j].cx[6]*pxx + g.domain[j].cy[6]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[6]*pxy + g.domain[j].cy[6]*pyy) ) ) / g.domain[j].T;
            g.domain[j].qstr[7] = g.domain[j].geq[7] + ( g.domain[j].w[7] * (2.0 * g.domain[j].ux * (g.domain[j].cx[7]*pxx + g.domain[j].cy[7]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[7]*pxy + g.domain[j].cy[7]*pyy) ) ) / g.domain[j].T;
            g.domain[j].qstr[8] = g.domain[j].geq[8] + ( g.domain[j].w[8] * (2.0 * g.domain[j].ux * (g.domain[j].cx[8]*pxx + g.domain[j].cy[8]*pyx) + 2.0 * g.domain[j].uy * (g.domain[j].cx[8]*pxy + g.domain[j].cy[8]*pyy) ) ) / g.domain[j].T; 
        }
    }
}
/*relaxation factors*/
void calcRelaxationFactors(grid& g, const double& sphr, const double& mu, const double& cv, const double& Pr){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j =0; j < g.NY * g.NX; ++j){
        if(g.domain[j].flag == 0){
            g.domain[j].omegaf  = 1.0 / ( (mu /(g.domain[j].rho * g.domain[j].T)) + 0.5 );
            g.domain[j].omegag  = 1.0 / ( (mu /(g.domain[j].rho * g.domain[j].T)) + 0.5 );
            g.domain[j].omegag2 = 1.0 / ( (mu/Pr) * (1.0 /(g.domain[j].rho * g.domain[j].T)) + 0.5 );//k/CP
            g.domain[j].tauf  = 1.0/ g.domain[j].omegaf;
            g.domain[j].taug  = 1.0/ g.domain[j].omegag;

        }
    }
}
/*Collision*/
void collision(grid& g){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY * g.NX; ++j){
        if(g.domain[j].flag == 0){

            g.domain[j].fold[0] = g.domain[j].fold[0] + g.domain[j].omegaf * (g.domain[j].feq[0] - g.domain[j].fold[0]);
            g.domain[j].fold[1] = g.domain[j].fold[1] + g.domain[j].omegaf * (g.domain[j].feq[1] - g.domain[j].fold[1]);
            g.domain[j].fold[2] = g.domain[j].fold[2] + g.domain[j].omegaf * (g.domain[j].feq[2] - g.domain[j].fold[2]);
            g.domain[j].fold[3] = g.domain[j].fold[3] + g.domain[j].omegaf * (g.domain[j].feq[3] - g.domain[j].fold[3]);
            g.domain[j].fold[4] = g.domain[j].fold[4] + g.domain[j].omegaf * (g.domain[j].feq[4] - g.domain[j].fold[4]);
            g.domain[j].fold[5] = g.domain[j].fold[5] + g.domain[j].omegaf * (g.domain[j].feq[5] - g.domain[j].fold[5]);
            g.domain[j].fold[6] = g.domain[j].fold[6] + g.domain[j].omegaf * (g.domain[j].feq[6] - g.domain[j].fold[6]);
            g.domain[j].fold[7] = g.domain[j].fold[7] + g.domain[j].omegaf * (g.domain[j].feq[7] - g.domain[j].fold[7]);
            g.domain[j].fold[8] = g.domain[j].fold[8] + g.domain[j].omegaf * (g.domain[j].feq[8] - g.domain[j].fold[8]);


            g.domain[j].gold[0] = g.domain[j].gold[0] + g.domain[j].omegag * (g.domain[j].geq[0] - g.domain[j].gold[0]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[0] - g.domain[j].geq[0]);
            g.domain[j].gold[1] = g.domain[j].gold[1] + g.domain[j].omegag * (g.domain[j].geq[1] - g.domain[j].gold[1]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[1] - g.domain[j].geq[1]);
            g.domain[j].gold[2] = g.domain[j].gold[2] + g.domain[j].omegag * (g.domain[j].geq[2] - g.domain[j].gold[2]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[2] - g.domain[j].geq[2]);
            g.domain[j].gold[3] = g.domain[j].gold[3] + g.domain[j].omegag * (g.domain[j].geq[3] - g.domain[j].gold[3]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[3] - g.domain[j].geq[3]);
            g.domain[j].gold[4] = g.domain[j].gold[4] + g.domain[j].omegag * (g.domain[j].geq[4] - g.domain[j].gold[4]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[4] - g.domain[j].geq[4]);
            g.domain[j].gold[5] = g.domain[j].gold[5] + g.domain[j].omegag * (g.domain[j].geq[5] - g.domain[j].gold[5]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[5] - g.domain[j].geq[5]);
            g.domain[j].gold[6] = g.domain[j].gold[6] + g.domain[j].omegag * (g.domain[j].geq[6] - g.domain[j].gold[6]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[6] - g.domain[j].geq[6]);
            g.domain[j].gold[7] = g.domain[j].gold[7] + g.domain[j].omegag * (g.domain[j].geq[7] - g.domain[j].gold[7]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[7] - g.domain[j].geq[7]);
            g.domain[j].gold[8] = g.domain[j].gold[8] + g.domain[j].omegag * (g.domain[j].geq[8] - g.domain[j].gold[8]) + (g.domain[j].omegag - g.domain[j].omegag2) * (g.domain[j].qstr[8] - g.domain[j].geq[8]);
        }
    }
}
/*shift assign*/
void shiftAsign(grid& g, double u_shift){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j =0;j<g.NY*g.NX;++j){
        g.domain[j].ushift = u_shift;
        for(int q=0; q<9;++q){
            g.domain[j].cx[q] = g.domain[j].cx[q] + u_shift;
            // g.domain[j].cy[q] = g.domain[j].cy[q] + 0.0;
        }
    }
}
/*reset fshift*/
void resetDDFshitS(grid& g){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j<g.NY*g.NX;++j){
        for(int q = 0; q<9;++q){
            /*interpolated values are assigned to this array and it will transfer back to fold and gold arrays*/
            g.domain[j].fshift[q] = g.domain[j].fold[q];
            g.domain[j].gshift[q] = g.domain[j].gold[q];
        }   
    }   
}

/*shift interpolation*/
void IDW(grid& g,double u_shift){
    /*need implementaition*/
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 1; j<g.NY-1; ++j){
        // if((u_shift > 0.0)){
            for(int i = 2; i < g.NX-1; ++i){ // from 2 
                if((u_shift > 0.0)){// && (g(j,i).flag == 0)
                    for(int q = 0; q < 9; ++q){
                        /*interpolated values are assigned to this array and it will transfer back to fold and gold arrays
                        here the shit is  in +x direction therefore extrapolation has to be done for the left fluid node near to the bounndary
                        therfore the the looping can be done from 2 to nx considering interpolation*/
                        g(j,i).fshift[q] = g(j,i-1).fold[q]*u_shift + g(j,i).fold[q]*(1.0 - u_shift);
                        g(j,i).gshift[q] = g(j,i-1).gold[q]*u_shift + g(j,i).gold[q]*(1.0 - u_shift);
                        /*for the cell near to left boundar extrapolation
                        we need to assume now all the old population are shifted by the shift velocity*/
                        // g(j,1).fshift[q] = (g(j,1).fold[q] * (1.0 + u_shift) + g(j,2).fold[q] * u_shift) / (1.0 + 2.0 * u_shift) ;
                        // g(j,1).gshift[q] = (g(j,1).gold[q] * (1.0 + u_shift) + g(j,2).gold[q] * u_shift) / (1.0 + 2.0 * u_shift) ;
                    }
                }
            }
        // }
    }
}
/*extrpolation*/
void IDWE(grid& g,double u_shift){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j<g.NY;++j){
        if((u_shift > 0.0)&&(g(j,1).flag == 0)){
                // if(g(j,1).flag == 0){
                    for(int q = 0; q < 9;++q){
                        g(j,1).fshift[q] = (g(j,1).fold[q] * (1.0 + u_shift)*1.0 + g(j,2).fold[q] * u_shift*1.0) / (1.0 + 2.0 * u_shift + 0.0) ;
                        g(j,1).gshift[q] = (g(j,1).gold[q] * (1.0 + u_shift)*1.0 + g(j,2).gold[q] * u_shift*1.0) / (1.0 + 2.0 * u_shift + 0.0) ; 
                    }
                }
        }
}
/*coorect the population values due to the shift*/
void interpolatedDDF(grid& g,double u_shift){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0;j < g.NY*g.NX;++j){
        if((u_shift > 0.0) ){//&& (g.domain[j].flag ==0)
            for(int q = 0; q < 9; ++q){
                g.domain[j].fold[q] = g.domain[j].fshift[q];
                g.domain[j].gold[q] = g.domain[j].gshift[q];
            }
        }
    }
}


/*Stream pull*/
void stream(grid& g){
    #pragma omp parallel
    { 
        #pragma omp for //schedule(static,CHUNK_SIZE)
        for(int j = 0; j < g.NY; ++j){
            for(int i = 0 ; i < g.NX; ++i){  //0 to nx-1

                if(g(j,i).flag == 0){
                    /*density population*/
                    g(j,i).fnew[0]   = g(j,i).fold[0];
                    g(j,i).fnew[1]   = g(j,i-1).fold[1];
                    g(j,i).fnew[2]   = g(j-1,i).fold[2];
                    g(j,i).fnew[3]   = g(j,i+1).fold[3];
                    g(j,i).fnew[4]   = g(j+1,i).fold[4];
                    g(j,i).fnew[5]   = g(j-1,i-1).fold[5];
                    g(j,i).fnew[6]   = g(j-1,i+1).fold[6];
                    g(j,i).fnew[7]   = g(j+1,i+1).fold[7];
                    g(j,i).fnew[8]   = g(j+1,i-1).fold[8];
                    /*energy population*/
                    g(j,i).gnew[0]   = g(j,i).gold[0];
                    g(j,i).gnew[1]   = g(j,i-1).gold[1];
                    g(j,i).gnew[2]   = g(j-1,i).gold[2];
                    g(j,i).gnew[3]   = g(j,i+1).gold[3];
                    g(j,i).gnew[4]   = g(j+1,i).gold[4];
                    g(j,i).gnew[5]   = g(j-1,i-1).gold[5];
                    g(j,i).gnew[6]   = g(j-1,i+1).gold[6];
                    g(j,i).gnew[7]   = g(j+1,i+1).gold[7];
                    g(j,i).gnew[8]   = g(j+1,i-1).gold[8];
    
                }
            }
        }
    }
}



void boundaryConditionInletOutlet(grid& g,grid& g_bc, double uin=0.0,double rhoout=1.0){
    // int endpX = (g.NX - 1);
    /*bounce back*/
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY ; ++j ){
        if(g(j,0).flag == 20){
            for(int q=0; q < 9; ++q){

                g(j,0).fold[q]  =  g_bc(j,1).feq[q];
                g(j,0).gold[q]  =  g_bc(j,1).geq[q];
 
            }
        }
        if(g(j,0).flag == 30){
            for(int q=0; q < 9; ++q){
                g(j,(g.NX - 1)).fold[q] = g(j,(g.NX - 1)-1).fold[q];
                g(j,(g.NX - 1)).gold[q] = g(j,(g.NX - 1)-1).gold[q];
                   
            }
        } 
    }
}


void boundaryConditiontWallLR(grid& g,grid& bg){
    //int endpX = g.NX - 1;
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 1; j < g.NY-1; ++j){//it was g.NY-1
        // if(g.domain[j].flag == 0){
        /*periodic bc*/
        // /*left*/
        //     g(j,0).  fold[1] = g(j,endpX-1).fold[1];
        //     g(j,0).  fold[5] = g(j,endpX-1).fold[5];
        //     g(j,0).  fold[8] = g(j,endpX-1).fold[8];

        //     g(j,0).gold[1]  = g(j,endpX-1).gold[1];
        //     g(j,0).gold[5]  = g(j,endpX-1).gold[5];
        //     g(j,0).gold[8]  = g(j,endpX-1).gold[8];
        // /*right*/
        //     g(j,endpX).  fold[3] = g(j,1).fold[3];
        //     g(j,endpX).  fold[6] = g(j,1).fold[6];
        //     g(j,endpX).  fold[7] = g(j,1).fold[7];
            
        //     g(j,endpX).  gold[3] = g(j,1).gold[3];
        //     g(j,endpX).  gold[6] = g(j,1).gold[6];
        //     g(j,endpX).  gold[7] = g(j,1).gold[7]; 
        /*bounce back*/
            // /*left*/
            // g(j,0).  fold[1] = g(j,1)  .fold[3];
            // g(j,0).  fold[5] = g(j+1,1).fold[7];
            // g(j,0).  fold[8] = g(j-1,1).fold[6];

            // g(j,0).  gold[1] = g(j,1)  .gold[3];
            // g(j,0).  gold[5] = g(j+1,1).gold[7];
            // g(j,0).  gold[8] = g(j-1,1).gold[6];
            // /*right*/
            // g(j,endpX).  fold[3] = g(j,endpX-1)  .fold[1];
            // g(j,endpX).  fold[6] = g(j+1,endpX-1).fold[8];
            // g(j,endpX).  fold[7] = g(j-1,endpX-1).fold[5];

            // g(j,endpX).  gold[3] = g(j,endpX-1)  .gold[1];
            // g(j,endpX).  gold[6] = g(j+1,endpX-1).gold[8];
            // g(j,endpX).  gold[7] = g(j-1,endpX-1).gold[5];

        /*dirichlet bc*/    
            for(int q = 0; q < 9; ++q){
                g(j,0).fold[q]     = g(j,1).fold[q];
                g(j,g.NX - 1).fold[q] = g(j,g.NX - 1-1).fold[q];
                g(j,0).gold[q]     = g(j,1).gold[q];
                g(j,g.NX - 1).gold[q] = g(j,g.NX - 1-1).gold[q]; 

                // g(j,1).fold[q]       = bg(j,1).feq[q];
                // g(j,endpX-1).fold[q] = bg(j,endpX-1).feq[q];
                // g(j,1).gold[q]       = bg(j,1).geq[q];
                // g(j,endpX-1).gold[q] = bg(j,endpX-1).geq[q];
            }
        }

        
        
    // }
}


void boundaryConditiontWallTP(grid& g){

    //int endpY = (g.NY - 1);
    // periodic
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int i = 0; i < g.NX; ++i){
        /*top*/
        g((g.NY - 1),i).  fold[4] = g(1,i)  .fold[4];
        g((g.NY - 1),i).  fold[7] = g(1,i)  .fold[7];
        g((g.NY - 1),i).  fold[8] = g(1,i)  .fold[8];

        g((g.NY - 1),i).  gold[4] = g(1,i)  .gold[4];
        g((g.NY - 1),i).  gold[7] = g(1,i)  .gold[7];
        g((g.NY - 1),i).  gold[8] = g(1,i)  .gold[8];
        /*bottom*/
        g(0,i).    fold[2] = g((g.NY - 1)-1,i)  .fold[2];
        g(0,i).    fold[5] = g((g.NY - 1)-1,i)  .fold[5];
        g(0,i).    fold[6] = g((g.NY - 1)-1,i)  .fold[6];

        g(0,i).    gold[2] = g((g.NY - 1)-1,i)  .gold[2];
        g(0,i).    gold[5] = g((g.NY - 1)-1,i)  .gold[5];
        g(0,i).    gold[6] = g((g.NY - 1)-1,i)  .gold[6];
    }
    /*****no slip wall****/
    // for(int i = 1; i < g.NX-1; ++i){
    //     /*top*/
    //     g(endpY,i).  fold[4] = g(endpY-1,i)  .fold[2];
    //     g(endpY,i).  fold[7] = g(endpY-1,i-1).fold[5];
    //     g(endpY,i).  fold[8] = g(endpY-1,i+1).fold[6];
    //     // /*bottom*/
    //     g(0,i).      fold[2] = g(1,i)  .fold[4];
    //     g(0,i).      fold[5] = g(1,i+1).fold[7];
    //     g(0,i).      fold[6] = g(1,i-1).fold[8];
    //     /*top*/
    //     g(endpY,i).  gold[4] = g(endpY-1,i)  .gold[2];
    //     g(endpY,i).  gold[7] = g(endpY-1,i-1).gold[5];
    //     g(endpY,i).  gold[8] = g(endpY-1,i+1).gold[6];
    //     // /*bottom*/
    //     g(0,i).      gold[2] = g(1,i)  .gold[4];
    //     g(0,i).      gold[5] = g(1,i+1).gold[7];
    //     g(0,i).      gold[6] = g(1,i-1).gold[8];
    // }
}
/*Boundary conditions for corners in the domain*/
void bcCorners(grid& g){
    int endpX = g.NX - 1;
    int endpY = g.NY - 1;
    g(0,0)        .fold[5] = g(1,1)      .fold[7];
    g(endpY,0)    .fold[8] = g(endpY-1,1).fold[6];
    g(0,endpX)    .fold[6] = g(1,endpX-1).fold[8];
    g(endpY,endpX).fold[7] = g(endpY-1,endpX-1).fold[5];

    g(0,0)        .gold[5] = g(1,1)      .gold[7];
    g(endpY,0)    .gold[8] = g(endpY-1,1).gold[6];
    g(0,endpX)    .gold[6] = g(1,endpX-1).gold[8];
    g(endpY,endpX).gold[7] = g(endpY-1,endpX-1).gold[5];

}
/*Knudsen number dependent stabilization*/
void KnudsenNumberStability(grid& g){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0 ; j <  g.NY * g.NX; ++j){
        double diffSumf = 0.0;
        double diffSumg = 0.0;

        if(g.domain[j].flag == 0){
            diffSumf =   abs(g.domain[j].fold[0] - g.domain[j].feq[0])/ g.domain[j].feq[0] +\
                                abs(g.domain[j].fold[1] - g.domain[j].feq[1])/ g.domain[j].feq[1] +\
                                abs(g.domain[j].fold[2] - g.domain[j].feq[2])/ g.domain[j].feq[2] +\
                                abs(g.domain[j].fold[3] - g.domain[j].feq[3])/ g.domain[j].feq[3] +\
                                abs(g.domain[j].fold[4] - g.domain[j].feq[4])/ g.domain[j].feq[4] +\
                                abs(g.domain[j].fold[5] - g.domain[j].feq[5])/ g.domain[j].feq[5] +\
                                abs(g.domain[j].fold[6] - g.domain[j].feq[6])/ g.domain[j].feq[6] +\
                                abs(g.domain[j].fold[7] - g.domain[j].feq[7])/ g.domain[j].feq[7] +\
                                abs(g.domain[j].fold[8] - g.domain[j].feq[8])/ g.domain[j].feq[8] ;

            diffSumg =   abs(g.domain[j].gold[0] - g.domain[j].geq[0])/ g.domain[j].geq[0] +\
                                abs(g.domain[j].gold[1] - g.domain[j].geq[1])/ g.domain[j].geq[1] +\
                                abs(g.domain[j].gold[2] - g.domain[j].geq[2])/ g.domain[j].geq[2] +\
                                abs(g.domain[j].gold[3] - g.domain[j].geq[3])/ g.domain[j].geq[3] +\
                                abs(g.domain[j].gold[4] - g.domain[j].geq[4])/ g.domain[j].geq[4] +\
                                abs(g.domain[j].gold[5] - g.domain[j].geq[5])/ g.domain[j].geq[5] +\
                                abs(g.domain[j].gold[6] - g.domain[j].geq[6])/ g.domain[j].geq[6] +\
                                abs(g.domain[j].gold[7] - g.domain[j].geq[7])/ g.domain[j].geq[7] +\
                                abs(g.domain[j].gold[8] - g.domain[j].geq[8])/ g.domain[j].geq[8] ;

            g.domain[j].epsilonf = (1.0/9.0) * diffSumf;  
            g.domain[j].epsilong = (1.0/9.0) * diffSumg;

            if(g.domain[j].epsilonf < 0.01){
                g.domain[j].alphaf = 1.0;
                g.domain[j].alphag = 1.0;}
            else if((g.domain[j].epsilonf >= 0.01) && (g.domain[j].epsilonf < 0.1)){
                g.domain[j].alphaf = 1.05;
                g.domain[j].alphag = 1.05;}
            else if((g.domain[j].epsilonf >= 0.1) && (g.domain[j].epsilonf < 1.0)){
                g.domain[j].alphaf = 1.35;
                g.domain[j].alphag = 1.35;}
            else if(g.domain[j].epsilonf >= 1.0){
                g.domain[j].alphaf = 1.0/g.domain[j].tauf;
                g.domain[j].alphag = 1.0/g.domain[j].taug;}

                

            g.domain[j].tauf = g.domain[j].tauf * g.domain[j].alphaf;
            g.domain[j].taug = g.domain[j].taug * g.domain[j].alphag;
            g.domain[j].omegaf = 1.0/g.domain[j].tauf;
            g.domain[j].omegag = 1.0/g.domain[j].taug;        
        }
  
    }
}
/*swapping old and new populations*/
void swap(grid& g){
    #pragma omp parallel for //schedule(static,CHUNK_SIZE)
    for(int j = 0; j < g.NY * g.NX ; ++j){
        for(int i = 0; i < 9; ++i){
            g.domain[j].fold[i] = g.domain[j].fnew[i];
            g.domain[j].gold[i] = g.domain[j].gnew[i];
        }
    } 
}

