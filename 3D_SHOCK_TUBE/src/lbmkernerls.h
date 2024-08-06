#pragma once  
#include <cmath>
#include <omp.h>
#include "math.h"
#include "node.h"
#include "grid.h"
#include "auxliray.h"
#include "state.h"

/*initialize lagrangians*/
void initializeLagrangians(grid& g){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        /*initialize the lgrangian multipliers if there's more multipliers initialize it too 
        if you want to initialize only the fluid cell just do it here*/
            g.domain[numberofNodes].LM[0] = 0.5;
            g.domain[numberofNodes].LM[1] = 0.5;
            g.domain[numberofNodes].LM[2] = 0.5;
            g.domain[numberofNodes].LM[3] = 0.5;
    }
}
/*initialize macroscopic parameters*/
void initializationMacro(grid& g,state d,int flagIndex){
    /*Initialize the domain accordingly all the cells in the domain is being initilized including the helper cells*/
    #pragma omp parallel for
    for(int k = 0; k < g.NZ; ++k){
        for(int j = 0; j < g.NY; ++j){
            for(int i =0; i < g.NX; ++i){
                if(g(k,j,i).flag == flagIndex){
                    g(k,j,i).ux  = d.u_x; 
                    g(k,j,i).uy  = d.u_y;
                    g(k,j,i).uz  = d.u_z;
                    // g(k,j,i).rho = rho0;
                    // g(k,j,i).T   = T0;
                    if(i <=  g.NX/2){
                        g(k,j,i).rho = d.rho_l;//0.5;
                        g(k,j,i).T   = d.T_l;//0.15;//0.2;//0.2;
                    }
                    if(i > g.NX/2){
                        g(k,j,i).rho = d.rho_r;//0.125;//2.0;
                        g(k,j,i).T   = d.T_r;//0.16;//0.025;
                    }
                }
            }
        }
    }
}
/*Initialize f populations*/
void initilizeDDF(grid& g){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        for(int q = 0; q < 27; ++q){ // now its 3d and we have 27 velocity sets
            // if(g.domain[j].flag ==0){
                g.domain[numberofNodes].fold[q]  = g.domain[numberofNodes].feq[q];
                g.domain[numberofNodes].fnew[q]  = g.domain[numberofNodes].feq[q];
                g.domain[numberofNodes].gold[q]  = g.domain[numberofNodes].geq[q];            
                g.domain[numberofNodes].gnew[q]  = g.domain[numberofNodes].geq[q];
            // }
        }
    }
}
/*set temperature dependent weights*/
void setWeights(grid& g,int flagIndex=0){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
            double zero     = 1.0 - g.domain[numberofNodes].T; 
            double plusOne  = g.domain[numberofNodes].T / 2.0;
            double minusOne = g.domain[numberofNodes].T / 2.0;
        // if(g.domain[j].flag == flagIndex){ //flag zeros is the fluid cells
            g.domain[numberofNodes].w[0] = zero     * zero     * zero;      //# 0 *  0 *  0 
            g.domain[numberofNodes].w[1] = plusOne  * zero     * zero;      //# 1 *  0 *  0 
            g.domain[numberofNodes].w[2] = minusOne * zero     * zero;      //#-1 *  0 *  0
            g.domain[numberofNodes].w[3] = zero     * plusOne  * zero;      //# 0 *  1 *  0 
            g.domain[numberofNodes].w[4] = zero     * minusOne * zero;      //# 0 * -1 *  0
            g.domain[numberofNodes].w[5] = zero     * zero     * plusOne;   //# 0 *  0 *  1
            g.domain[numberofNodes].w[6] = zero     * zero     * minusOne;  //# 0 *  0 * -1
            g.domain[numberofNodes].w[7] = plusOne  * plusOne  * zero;      //# 1 *  1 *  0
            g.domain[numberofNodes].w[8] = minusOne * minusOne * zero;      //#-1 * -1 *  0
        //
            g.domain[numberofNodes].w[9] = plusOne  * zero     * plusOne;   //# 1 *  0 *  1 
            g.domain[numberofNodes].w[10]= minusOne * zero     * minusOne;  //#-1 *  0 * -1 
            g.domain[numberofNodes].w[11]= zero     * plusOne  * plusOne;   //# 0 *  1 *  1
            g.domain[numberofNodes].w[12]= zero     * minusOne * minusOne;  //# 0 * -1 * -1 
            g.domain[numberofNodes].w[13]= plusOne  * minusOne * zero;      //# 1 * -1 *  0
            g.domain[numberofNodes].w[14]= minusOne * plusOne  * zero;      //#-1 *  1 *  0
            g.domain[numberofNodes].w[15]= plusOne  * zero     * minusOne;  //# 1 *  0 * -1
            g.domain[numberofNodes].w[16]= minusOne * zero     * plusOne;   //#-1 *  0 *  1
            g.domain[numberofNodes].w[17]= zero     * plusOne  * minusOne;  //# 0 *  1 * -1
        // 
            g.domain[numberofNodes].w[18]= zero     * minusOne  * plusOne;   //# 0 * -1 *  1 
            g.domain[numberofNodes].w[19]= plusOne  * plusOne   * plusOne;   //# 1 *  1 *  1 
            g.domain[numberofNodes].w[20]= minusOne * minusOne  * minusOne;  //#-1 * -1 * -1
            g.domain[numberofNodes].w[21]= plusOne  * plusOne   * minusOne;  //# 1 *  1 * -1 
            g.domain[numberofNodes].w[22]= minusOne * minusOne  * plusOne;   //#-1 * -1 *  1
            g.domain[numberofNodes].w[23]= plusOne  * minusOne  * plusOne;   //# 1 * -1 *  1
            g.domain[numberofNodes].w[24]= minusOne * plusOne   * minusOne;  //#-1 *  1 * -1
            g.domain[numberofNodes].w[25]= minusOne * plusOne   * plusOne;   //#-1 *  1 *  1
            g.domain[numberofNodes].w[26]= plusOne  * minusOne  * minusOne;  //# 1 * -1 * -1
        // }
    }
}
/*calculate the macroscopic parameters*/
/*Density*/
void density(grid& g,int flagIndex){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == flagIndex){
            double dSum = 0.0;
            for(int q = 0; q < 27; ++q){
                dSum = dSum + g.domain[numberofNodes].fold[q];
            }
            g.domain[numberofNodes].rho = dSum;
            if(isnan(g.domain[numberofNodes].rho)){std::cout << "-nan- encountrerd in rho" << std::endl;}
        }
    }
}
/*velocity*/
void velocity(grid& g,int flagIndex){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ*g.NY*g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == flagIndex){
            double uxSum = 0.0;
            double uySum = 0.0;
            double uzSum = 0.0;
            for(int q = 0; q < 27; ++q){
                uxSum = uxSum + g.domain[numberofNodes].fold[q] * g.domain[numberofNodes].cx[q];
                uySum = uySum + g.domain[numberofNodes].fold[q] * g.domain[numberofNodes].cy[q];
                uzSum = uzSum + g.domain[numberofNodes].fold[q] * g.domain[numberofNodes].cz[q];
            }
            g.domain[numberofNodes].ux = uxSum / g.domain[numberofNodes].rho;
            if(isnan(g.domain[numberofNodes].ux)){std::cout << "-nan- encountrerd in ux" << std::endl;}

            g.domain[numberofNodes].uy = uySum / g.domain[numberofNodes].rho;
            if(isnan(g.domain[numberofNodes].uy)){std::cout << "-nan- encountrerd in uy" << std::endl;}

            g.domain[numberofNodes].uz = uzSum / g.domain[numberofNodes].rho;
            if(isnan(g.domain[numberofNodes].uz)){std::cout << "-nan- encountrerd in uz" << std::endl;}
        }
    }
}
/*Density*/
void temperature(grid& g,const double& cv,int flagIndex){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ*g.NY*g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == flagIndex){
            double ESum = 0.0;
            for(int q = 0; q < 27; ++q){
                ESum = ESum + g.domain[numberofNodes].gold[q];
            }
            double E      = ESum / (2.0 * g.domain[numberofNodes].rho);
            double velMag =   g.domain[numberofNodes].ux * g.domain[numberofNodes].ux 
                            + g.domain[numberofNodes].uy * g.domain[numberofNodes].uy 
                            + g.domain[numberofNodes].uz * g.domain[numberofNodes].uz ;
            g.domain[numberofNodes].T = ( E - (velMag * 0.5) ) / cv;
            if(isnan(g.domain[numberofNodes].T)){std::cout << "-nan- encountrerd in T" << std::endl;}
        }
    }
}
/*helpers to calculate F population*/
double zeroF(double v, double vshift, double t){
    return 1.0 - (v - vshift)*(v - vshift) - t;
}
double plusOneF(double v, double vshift, double t){
    return (  1.0*(v-vshift) + (v-vshift)*(v-vshift) + t ) / 2.0;
}
double minusOneF(double v, double vshift, double t){
    return ( -1.0*(v-vshift) + (v-vshift)*(v-vshift) + t ) / 2.0;
}
/*calculate F equilibrium population*/
void calcFeq(grid& g,int flagIndex){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == flagIndex){
            double ax = zeroF    (g.domain[numberofNodes].ux,g.domain[numberofNodes].ushift,g.domain[numberofNodes].T); //zero
            double bx = plusOneF (g.domain[numberofNodes].ux,g.domain[numberofNodes].ushift,g.domain[numberofNodes].T); //plusOne
            double cx = minusOneF(g.domain[numberofNodes].ux,g.domain[numberofNodes].ushift,g.domain[numberofNodes].T); //minusOne

            double ay = zeroF    (g.domain[numberofNodes].uy,0.0,g.domain[numberofNodes].T);    //zero
            double by = plusOneF (g.domain[numberofNodes].uy,0.0,g.domain[numberofNodes].T);    //plusOne
            double cy = minusOneF(g.domain[numberofNodes].uy,0.0,g.domain[numberofNodes].T);    //minusOne

            double az = zeroF    (g.domain[numberofNodes].uz,0.0,g.domain[numberofNodes].T);    //zero
            double bz = plusOneF (g.domain[numberofNodes].uz,0.0,g.domain[numberofNodes].T);    //plusOne
            double cz = minusOneF(g.domain[numberofNodes].uz,0.0,g.domain[numberofNodes].T);    //minusOne

            g.domain[numberofNodes].feq[0] =  g.domain[numberofNodes].rho * ax * ay * az;     //# 0 *  0 *  0 
            g.domain[numberofNodes].feq[1] =  g.domain[numberofNodes].rho * bx * ay * az;     //# 1 *  0 *  0 
            g.domain[numberofNodes].feq[2] =  g.domain[numberofNodes].rho * cx * ay * az;     //#-1 *  0 *  0
            g.domain[numberofNodes].feq[3] =  g.domain[numberofNodes].rho * ax * by * az;     //# 0 *  1 *  0 
            g.domain[numberofNodes].feq[4] =  g.domain[numberofNodes].rho * ax * cy * az;     //# 0 * -1 *  0
            g.domain[numberofNodes].feq[5] =  g.domain[numberofNodes].rho * ax * ay * bz;     //# 0 *  0 *  1
            g.domain[numberofNodes].feq[6] =  g.domain[numberofNodes].rho * ax * ay * cz;     //# 0 *  0 * -1
            g.domain[numberofNodes].feq[7] =  g.domain[numberofNodes].rho * bx * by * az;     //# 1 *  1 *  0
            g.domain[numberofNodes].feq[8] =  g.domain[numberofNodes].rho * cx * cy * az;     //#-1 * -1 *  0
        //
            g.domain[numberofNodes].feq[9] =  g.domain[numberofNodes].rho * bx * ay * bz;     //# 1 *  0 *  1 
            g.domain[numberofNodes].feq[10]=  g.domain[numberofNodes].rho * cx * ay * cz;     //#-1 *  0 * -1 
            g.domain[numberofNodes].feq[11]=  g.domain[numberofNodes].rho * ax * by * bz;     //# 0 *  1 *  1
            g.domain[numberofNodes].feq[12]=  g.domain[numberofNodes].rho * ax * cy * cz;     //# 0 * -1 * -1 
            g.domain[numberofNodes].feq[13]=  g.domain[numberofNodes].rho * bx * cy * az;     //# 1 * -1 *  0
            g.domain[numberofNodes].feq[14]=  g.domain[numberofNodes].rho * cx * by * az;     //#-1 *  1 *  0
            g.domain[numberofNodes].feq[15]=  g.domain[numberofNodes].rho * bx * ay * cz;     //# 1 *  0 * -1
            g.domain[numberofNodes].feq[16]=  g.domain[numberofNodes].rho * cx * ay * bz;     //#-1 *  0 *  1
            g.domain[numberofNodes].feq[17]=  g.domain[numberofNodes].rho * ax * by * cz;     //# 0 *  1 * -1
        // 
            g.domain[numberofNodes].feq[18]=  g.domain[numberofNodes].rho * ax * cy * bz;   //# 0 * -1 *  1 
            g.domain[numberofNodes].feq[19]=  g.domain[numberofNodes].rho * bx * by * bz;   //# 1 *  1 *  1 
            g.domain[numberofNodes].feq[20]=  g.domain[numberofNodes].rho * cx * cy * cz;   //#-1 * -1 * -1
            g.domain[numberofNodes].feq[21]=  g.domain[numberofNodes].rho * bx * by * cz;   //# 1 *  1 * -1 
            g.domain[numberofNodes].feq[22]=  g.domain[numberofNodes].rho * cx * cy * bz;   //#-1 * -1 *  1
            g.domain[numberofNodes].feq[23]=  g.domain[numberofNodes].rho * bx * cy * bz;   //# 1 * -1 *  1
            g.domain[numberofNodes].feq[24]=  g.domain[numberofNodes].rho * cx * by * cz;   //#-1 *  1 * -1
            g.domain[numberofNodes].feq[25]=  g.domain[numberofNodes].rho * cx * by * bz;   //#-1 *  1 *  1
            g.domain[numberofNodes].feq[26]=  g.domain[numberofNodes].rho * bx * cy * cz;   //# 1 * -1 * -1      
        }
    }
}
/*calculate G equilibrium population*/
void calcGeq(grid& g,int flagIndex){
    #pragma omp parallel for
    for(int  numberofNodes= 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == flagIndex){
            for(int q = 0; q < 27; ++q){
                g.domain[numberofNodes].geq[q] = g.domain[numberofNodes].rho * g.domain[numberofNodes].w[q] 
                                                                             * exp( 1.0*(g.domain[numberofNodes].LM[0] 
                                                                                       + g.domain[numberofNodes].LM[1] * g.domain[numberofNodes].cx[q] 
                                                                                       + g.domain[numberofNodes].LM[2] * g.domain[numberofNodes].cy[q] 
                                                                                       + g.domain[numberofNodes].LM[3] * g.domain[numberofNodes].cz[q]) ); 
            }
        }
    }
}
/*calculate Pressure tensor*/
void pTensor(grid& g){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == 0){
            //pressure tensor
            double p0Sum   = 0.0;double p1Sum   = 0.0;double p2Sum   = 0.0;double p3Sum   = 0.0;double p4Sum   = 0.0;
            double p5Sum   = 0.0;double p6Sum   = 0.0;double p7Sum   = 0.0;double p8Sum   = 0.0;
            //equilibrium pressure tensor
            double p0eqSum = 0.0;double p1eqSum = 0.0;double p2eqSum = 0.0;double p3eqSum = 0.0;double p4eqSum = 0.0;
            double p5eqSum = 0.0;double p6eqSum = 0.0;double p7eqSum = 0.0;double p8eqSum = 0.0;
            for(int q=0; q < 27; ++q){
                p0Sum   = p0Sum + g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].fold[q];    //xx
                p1Sum   = p1Sum + g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].fold[q];    //xy
                p2Sum   = p2Sum + g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].fold[q];    //xz
                p3Sum   = p3Sum + g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].fold[q];    //yx
                p4Sum   = p4Sum + g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].fold[q];    //yy
                p5Sum   = p5Sum + g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].fold[q];    //yz
                p6Sum   = p6Sum + g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].fold[q];    //zx
                p7Sum   = p7Sum + g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].fold[q];    //zy
                p8Sum   = p8Sum + g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].fold[q];    //zz
            // 
                p0eqSum = p0eqSum + g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].fold[q];    //xx
                p1eqSum = p1eqSum + g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].fold[q];    //xy
                p2eqSum = p2eqSum + g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].fold[q];    //xz
                p3eqSum = p3eqSum + g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].fold[q];    //yx
                p4eqSum = p4eqSum + g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].fold[q];    //yy
                p5eqSum = p5eqSum + g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].fold[q];    //yz
                p6eqSum = p6eqSum + g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].cx[q] * g.domain[numberofNodes].fold[q];    //zx
                p7eqSum = p7eqSum + g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].cy[q] * g.domain[numberofNodes].fold[q];    //zy
                p8eqSum = p8eqSum + g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].cz[q] * g.domain[numberofNodes].fold[q];    //zz

            }
            g.domain[numberofNodes].P[0]   = p0Sum;
            g.domain[numberofNodes].P[1]   = p1Sum;
            g.domain[numberofNodes].P[2]   = p2Sum;
            g.domain[numberofNodes].P[3]   = p3Sum;
            g.domain[numberofNodes].P[4]   = p4Sum;
            g.domain[numberofNodes].P[5]   = p5Sum;
            g.domain[numberofNodes].P[6]   = p6Sum;
            g.domain[numberofNodes].P[7]   = p7Sum;
            g.domain[numberofNodes].P[8]   = p8Sum;
            // 
            g.domain[numberofNodes].Peq[0] = p0eqSum;
            g.domain[numberofNodes].Peq[1] = p1eqSum;
            g.domain[numberofNodes].Peq[2] = p2eqSum;
            g.domain[numberofNodes].Peq[3] = p3eqSum;
            g.domain[numberofNodes].Peq[4] = p4eqSum;
            g.domain[numberofNodes].Peq[5] = p5eqSum;
            g.domain[numberofNodes].Peq[6] = p6eqSum;
            g.domain[numberofNodes].Peq[7] = p7eqSum;
            g.domain[numberofNodes].Peq[8] = p8eqSum;
        }
    }
}
/*Quasi Euilibrium G*/
void calcQuasiEqG(grid& g){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == 0){
            //P = P-Peq
            double pxx = g.domain[numberofNodes].P[0] - g.domain[numberofNodes].Peq[0];
            double pxy = g.domain[numberofNodes].P[1] - g.domain[numberofNodes].Peq[1];
            double pxz = g.domain[numberofNodes].P[2] - g.domain[numberofNodes].Peq[2];
            double pyx = g.domain[numberofNodes].P[3] - g.domain[numberofNodes].Peq[3];
            double pyy = g.domain[numberofNodes].P[4] - g.domain[numberofNodes].Peq[4];
            double pyz = g.domain[numberofNodes].P[5] - g.domain[numberofNodes].Peq[5];
            double pzx = g.domain[numberofNodes].P[6] - g.domain[numberofNodes].Peq[6];
            double pzy = g.domain[numberofNodes].P[7] - g.domain[numberofNodes].Peq[7];
            double pzz = g.domain[numberofNodes].P[8] - g.domain[numberofNodes].Peq[8];
            // 
            for(int q=0; q < 27; ++q){
                g.domain[numberofNodes].qstr[q] = g.domain[numberofNodes].geq[q] + ( 2.0 * g.domain[numberofNodes].w[q] / g.domain[numberofNodes].T ) * (  g.domain[numberofNodes].ux * pxx * g.domain[numberofNodes].cx[q] + g.domain[numberofNodes].uy * pxy * g.domain[numberofNodes].cx[q] + g.domain[numberofNodes].uz * pxz * g.domain[numberofNodes].cx[q]
                                                                                                                                                          +g.domain[numberofNodes].ux * pyx * g.domain[numberofNodes].cy[q] + g.domain[numberofNodes].uy * pyy * g.domain[numberofNodes].cy[q] + g.domain[numberofNodes].uz * pyz * g.domain[numberofNodes].cy[q]
                                                                                                                                                          +g.domain[numberofNodes].ux * pzx * g.domain[numberofNodes].cz[q] + g.domain[numberofNodes].uy * pzy * g.domain[numberofNodes].cz[q] + g.domain[numberofNodes].uz * pzz * g.domain[numberofNodes].cz[q]);
            }
        }
    }
}
/*relaxation factors*/
void calcRelaxationFactors(grid& g, const double& sphr, const double& mu, const double& cv, const double& Pr){
    // cv = 1.0 /(SPHR - 1.0)
    // double CP = cv + 1.0;
    // double k  = CP * mu / Pr;
    #pragma omp parallel for
    for(int numberofNodes =0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == 0){
            g.domain[numberofNodes].omegaf  = 1.0 / ( (mu /(g.domain[numberofNodes].rho * g.domain[numberofNodes].T)) + 0.5 );
            g.domain[numberofNodes].omegag  = 1.0 / ( (mu /(g.domain[numberofNodes].rho * g.domain[numberofNodes].T)) + 0.5 );
            g.domain[numberofNodes].omegag2 = 1.0 / ( (mu/Pr) * (1.0 /(g.domain[numberofNodes].rho * g.domain[numberofNodes].T)) + 0.5 );
            g.domain[numberofNodes].tauf    = 1.0 / g.domain[numberofNodes].omegaf;
            g.domain[numberofNodes].taug    = 1.0 / g.domain[numberofNodes].omegag;

        }
    }
}
/*Collision*/
void collision(grid& g){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        if(g.domain[numberofNodes].flag == 0){
            for(int q=0; q < 27; ++q){
                g.domain[numberofNodes].fold[q] = g.domain[numberofNodes].fold[q] + g.domain[numberofNodes].omegaf * (g.domain[numberofNodes].feq[q] - g.domain[numberofNodes].fold[q]);
                g.domain[numberofNodes].gold[q] = g.domain[numberofNodes].gold[q] + g.domain[numberofNodes].omegag * (g.domain[numberofNodes].geq[q] - g.domain[numberofNodes].gold[q]) + (g.domain[numberofNodes].omegag - g.domain[numberofNodes].omegag2) * (g.domain[numberofNodes].qstr[q] - g.domain[numberofNodes].geq[q]);
            }
        }
    }
}
// /*shift assign*/
void shiftAsign(grid& g, double u_shift){
    #pragma omp parallel for
    for(int numberofNodes =0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        g.domain[numberofNodes].ushift = u_shift;
        for(int q = 0; q < 27; ++q){
            g.domain[numberofNodes].cx[q] = g.domain[numberofNodes].cx[q] + u_shift;
            // g.domain[j].cy[q] = g.domain[j].cy[q] + 0.0;
        }
    }
}
// /*reset fshift*/
void resetDDFshitS(grid& g){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        for(int q = 0; q < 27; ++q){
            /*interpolated values are assigned to this array and it will transfer back to fold and gold arrays*/
            g.domain[numberofNodes].fshift[q] = g.domain[numberofNodes].fold[q];
            g.domain[numberofNodes].gshift[q] = g.domain[numberofNodes].gold[q];
        }   
    }   
}

/*shift interpolation*/
void IDW(grid& g,double u_shift){
    /*need implementaition*/
    #pragma omp parallel for
    for(int k = 1; k < g.NZ-1; ++k){
        for(int j = 1; j <g.NY-1; ++j){
            // if((u_shift > 0.0)){
                for(int i = 2; i < g.NX-1; ++i){ // from 2 to last fluid node because the shift is positive until now
                    if((u_shift > 0.0)){// && (g(j,i).flag == 0)
                        for(int q = 0; q < 27; ++q){
                            /*interpolated values are assigned to this array and it will transfer back to fold and gold arrays
                            here the shit is  in +x direction therefore extrapolation has to be done for the left fluid node near to the bounndary
                            therfore the the looping can be done from 2 to nx considering interpolation*/
                            g(k,j,i).fshift[q] = g(k,j,i-1).fold[q]*u_shift + g(k,j,i).fold[q]*(1.0 - u_shift);
                            g(k,j,i).gshift[q] = g(k,j,i-1).gold[q]*u_shift + g(k,j,i).gold[q]*(1.0 - u_shift);

                        }
                    }
                }
            // }
        }
    }
}
/*extrpolation*/
void IDWE(grid& g,double u_shift){
    /*need implementaition*/
    #pragma omp parallel for
    for(int k = 0; k < g.NZ; ++k){
        for(int j = 0; j < g.NY; ++j){
            if((u_shift > 0.0)&&(g(k,j,1).flag == 0)){
                    // if(g(j,1).flag == 0){
                        for(int q = 0; q < 27; ++q){
                            g(k,j,1).fshift[q] = (g(k,j,1).fold[q] * (1.0 + u_shift)*1.0 + g(k,j,2).fold[q] * u_shift*1.0) / (1.0 + 2.0 * u_shift + 0.0) ;
                            g(k,j,1).gshift[q] = (g(k,j,1).gold[q] * (1.0 + u_shift)*1.0 + g(k,j,2).gold[q] * u_shift*1.0) / (1.0 + 2.0 * u_shift + 0.0) ;

                        }
                    }
            }
    }
}
/*coorect the population values due to the shift*/
void interpolatedDDF(grid& g,double u_shift){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        if((u_shift > 0.0) ){//&& (g.domain[j].flag ==0)
            for(int q = 0; q < 27; ++q){
                g.domain[numberofNodes].fold[q] = g.domain[numberofNodes].fshift[q];
                g.domain[numberofNodes].gold[q] = g.domain[numberofNodes].gshift[q];
            }
        }
    }
}


/*Stream pull*/
void stream(grid& g){
    #pragma omp parallel for
    for(int k = 0; k < g.NZ; ++k){
        for(int j = 0; j < g.NY; ++j){
            for(int i = 0 ; i < g.NX; ++i){  //0 to nx-1
                if(g(k,j,i).flag == 0){
                    /*density population*/
                    g(k,j,i).fnew[0]   = g(k,j,i)    .fold[0];
                    g(k,j,i).fnew[1]   = g(k,j,i-1)  .fold[1];
                    g(k,j,i).fnew[2]   = g(k,j,i+1)  .fold[2];
                    g(k,j,i).fnew[3]   = g(k,j-1,i)  .fold[3];
                    g(k,j,i).fnew[4]   = g(k,j+1,i)  .fold[4];
                    g(k,j,i).fnew[5]   = g(k-1,j,i)  .fold[5];
                    g(k,j,i).fnew[6]   = g(k+1,j,i)  .fold[6];
                    g(k,j,i).fnew[7]   = g(k,j-1,i-1).fold[7];
                    g(k,j,i).fnew[8]   = g(k,j+1,i+1).fold[8];
                    // 
                    g(k,j,i).fnew[9]    = g(k-1,j,i-1).fold[9]   ;
                    g(k,j,i).fnew[10]   = g(k+1,j,i+1).fold[10]  ;
                    g(k,j,i).fnew[11]   = g(k-1,j-1,i).fold[11]  ;
                    g(k,j,i).fnew[12]   = g(k+1,j+1,i).fold[12]  ;
                    g(k,j,i).fnew[13]   = g(k,j+1,i-1).fold[13]  ;
                    g(k,j,i).fnew[14]   = g(k,j-1,i+1).fold[14]  ;
                    g(k,j,i).fnew[15]   = g(k+1,j,i-1).fold[15]  ;
                    g(k,j,i).fnew[16]   = g(k-1,j,i+1).fold[16]  ;
                    g(k,j,i).fnew[17]   = g(k+1,j-1,i).fold[17]  ;
                    // 
                    g(k,j,i).fnew[18]   = g(k-1,j+1,i)  .fold[18]; //new it should be old
                    g(k,j,i).fnew[19]   = g(k-1,j-1,i-1).fold[19];
                    g(k,j,i).fnew[20]   = g(k+1,j+1,i+1).fold[20];
                    g(k,j,i).fnew[21]   = g(k+1,j-1,i-1).fold[21];
                    g(k,j,i).fnew[22]   = g(k-1,j+1,i+1).fold[22];
                    g(k,j,i).fnew[23]   = g(k-1,j+1,i-1).fold[23];
                    g(k,j,i).fnew[24]   = g(k+1,j-1,i+1).fold[24];
                    g(k,j,i).fnew[25]   = g(k-1,j-1,i+1).fold[25];
                    g(k,j,i).fnew[26]   = g(k+1,j+1,i-1).fold[26];
                    /*energy population*/
                    g(k,j,i).gnew[0]   = g(k,j,i)    .gold[0];
                    g(k,j,i).gnew[1]   = g(k,j,i-1)  .gold[1];
                    g(k,j,i).gnew[2]   = g(k,j,i+1)  .gold[2];
                    g(k,j,i).gnew[3]   = g(k,j-1,i)  .gold[3];
                    g(k,j,i).gnew[4]   = g(k,j+1,i)  .gold[4];
                    g(k,j,i).gnew[5]   = g(k-1,j,i)  .gold[5];
                    g(k,j,i).gnew[6]   = g(k+1,j,i)  .gold[6];
                    g(k,j,i).gnew[7]   = g(k,j-1,i-1).gold[7];
                    g(k,j,i).gnew[8]   = g(k,j+1,i+1).gold[8];
                    // 
                    g(k,j,i).gnew[9]    = g(k-1,j,i-1).gold[9]   ;
                    g(k,j,i).gnew[10]   = g(k+1,j,i+1).gold[10]  ;
                    g(k,j,i).gnew[11]   = g(k-1,j-1,i).gold[11]  ;
                    g(k,j,i).gnew[12]   = g(k+1,j+1,i).gold[12]  ;
                    g(k,j,i).gnew[13]   = g(k,j+1,i-1).gold[13]  ;
                    g(k,j,i).gnew[14]   = g(k,j-1,i+1).gold[14]  ;
                    g(k,j,i).gnew[15]   = g(k+1,j,i-1).gold[15]  ;
                    g(k,j,i).gnew[16]   = g(k-1,j,i+1).gold[16]  ;
                    g(k,j,i).gnew[17]   = g(k+1,j-1,i).gold[17]  ;
                    // 
                    g(k,j,i).gnew[18]   = g(k-1,j+1,i)  .gold[18];
                    g(k,j,i).gnew[19]   = g(k-1,j-1,i-1).gold[19];
                    g(k,j,i).gnew[20]   = g(k+1,j+1,i+1).gold[20];
                    g(k,j,i).gnew[21]   = g(k+1,j-1,i-1).gold[21];
                    g(k,j,i).gnew[22]   = g(k-1,j+1,i+1).gold[22];
                    g(k,j,i).gnew[23]   = g(k-1,j+1,i-1).gold[23];
                    g(k,j,i).gnew[24]   = g(k+1,j-1,i+1).gold[24];
                    g(k,j,i).gnew[25]   = g(k-1,j-1,i+1).gold[25];
                    g(k,j,i).gnew[26]   = g(k+1,j+1,i-1).gold[26];
                }
            }
        }
    }
}
/*inlet outlet boundary condition we normally use it if we have inlet flow scenarios*/
void boundaryConditionInletOutlet(grid& g,grid& g_bc, double uin=0.0,double rhoout=1.0){
    // int endpX = g.NX - 1;
    // int endpY = g.NY - 1;
    // int endpZ = g.NZ - 1;
    #pragma omp parallel for
    for(int k = 0; k < g.NZ ; ++k )
        for(int j = 0; j < g.NY ; ++j ){
            /*inlet*/
            if(g(k,j,0).flag == 20){
                for(int q = 0; q < 27; ++q){

                    g(k,j,0).fold[q]  =  g_bc(k,j,1).feq[q];
                    g(k,j,0).gold[q]  =  g_bc(k,j,1).geq[q];
 
                }
            }
            /*outlet*/
            if(g(k,j,0).flag == 30){
                for(int q = 0; q < 27; ++q){
                    g(k,j,(g.NX - 1)).fold[q] = g(k,j,(g.NX - 1)-1).fold[q];
                    g(k,j,(g.NX - 1)).gold[q] = g(k,j,(g.NX - 1)-1).gold[q];

                    
                }
            } 
        }
}
/*left right boundary conditions*/
void boundaryConditiontWallLR(grid& g,grid& bg){
    // int endpX = (g.NX - 1);
    // int endpY = g.NY - 1;
    // int endpZ = g.NZ - 1;
    #pragma omp parallel for
    for(int k = 1; k < g.NZ-1; ++k){
        for(int j = 1; j < g.NY-1; ++j){//it was 1 to g.NY-1
            // if(g.domain[j].flag == 0){
            /*periodic bc*/
            // /*left*/
                // g(k,j,0).  fold[1]  = g(k,j,endpX-1).fold[1] ;
                // g(k,j,0).  fold[7]  = g(k,j,endpX-1).fold[7] ;
                // g(k,j,0).  fold[9]  = g(k,j,endpX-1).fold[9] ;
                // g(k,j,0).  fold[13] = g(k,j,endpX-1).fold[13];
                // g(k,j,0).  fold[15] = g(k,j,endpX-1).fold[15];
                // g(k,j,0).  fold[19] = g(k,j,endpX-1).fold[19];
                // g(k,j,0).  fold[21] = g(k,j,endpX-1).fold[21];
                // g(k,j,0).  fold[23] = g(k,j,endpX-1).fold[23];
                // g(k,j,0).  fold[26] = g(k,j,endpX-1).fold[26];
                // // 
                // g(k,j,0).  gold[1]  = g(k,j,endpX-1).gold[1] ;
                // g(k,j,0).  gold[7]  = g(k,j,endpX-1).gold[7] ;
                // g(k,j,0).  gold[9]  = g(k,j,endpX-1).gold[9] ;
                // g(k,j,0).  gold[13] = g(k,j,endpX-1).gold[13];
                // g(k,j,0).  gold[15] = g(k,j,endpX-1).gold[15];
                // g(k,j,0).  gold[19] = g(k,j,endpX-1).gold[19];
                // g(k,j,0).  gold[21] = g(k,j,endpX-1).gold[21];
                // g(k,j,0).  gold[23] = g(k,j,endpX-1).gold[23];
                // g(k,j,0).  gold[26] = g(k,j,endpX-1).gold[26];

            /*right*/
                // g(k,j,endpX).  fold[2]  = g(k,j,1).fold[2] ;
                // g(k,j,endpX).  fold[8]  = g(k,j,1).fold[8] ;
                // g(k,j,endpX).  fold[10] = g(k,j,1).fold[10];
                // g(k,j,endpX).  fold[14] = g(k,j,1).fold[14];
                // g(k,j,endpX).  fold[16] = g(k,j,1).fold[16];
                // g(k,j,endpX).  fold[20] = g(k,j,1).fold[20];
                // g(k,j,endpX).  fold[22] = g(k,j,1).fold[22];
                // g(k,j,endpX).  fold[24] = g(k,j,1).fold[24];
                // g(k,j,endpX).  fold[25] = g(k,j,1).fold[25];
                // // 
                // g(k,j,endpX).  gold[2]  = g(k,j,1).gold[2] ;
                // g(k,j,endpX).  gold[8]  = g(k,j,1).gold[8] ;
                // g(k,j,endpX).  gold[10] = g(k,j,1).gold[10];
                // g(k,j,endpX).  gold[14] = g(k,j,1).gold[14];
                // g(k,j,endpX).  gold[16] = g(k,j,1).gold[16];
                // g(k,j,endpX).  gold[20] = g(k,j,1).gold[20];
                // g(k,j,endpX).  gold[22] = g(k,j,1).gold[22];
                // g(k,j,endpX).  gold[24] = g(k,j,1).gold[24];
                // g(k,j,endpX).  gold[25] = g(k,j,1).gold[25];
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
                for(int q = 0; q < 27; ++q){
                    // Distributiion for mass and moementun
                    g(k,j,0).fold[q]     = g(k,j,1).fold[q];
                    g(k,j,(g.NX - 1)).fold[q] = g(k,j,(g.NX - 1)-1).fold[q];
                    // Distribution for energy
                    g(k,j,0).gold[q]     = g(k,j,1).gold[q];
                    g(k,j,(g.NX - 1)).gold[q] = g(k,j,(g.NX - 1)-1).gold[q]; 
                    // g(j,1).fold[q]       = bg(j,1).feq[q];
                    // g(j,endpX-1).fold[q] = bg(j,endpX-1).feq[q];
                    // g(j,1).gold[q]       = bg(j,1).geq[q];
                    // g(j,endpX-1).gold[q] = bg(j,endpX-1).geq[q];
                }
        }
    }

        
        
    // }
}
/*top bottom boundary*/
void boundaryConditiontWallTP(grid& g){
    // int endpX = g.NX - 1;
    // int endpY = (g.NY - 1);
    // int endpZ = g.NZ - 1;
    #pragma omp parallel for
    for(int k = 0; k < g.NZ; ++k){
        for(int i = 0; i < g.NX; ++i){
            /*First order Neumann*/
            // for(int q = 0; q < 27; ++q){
            //         // double FEQW1     =  feqForBoundary(g(1,1).w[q], 1.0, 0.01, 0, g(1,1).cx[q], g(1,1).cy[q]);
            //         g(k,endpY,i).  fold[q] = g(k,endpY-1,i).  fold[q];//FEQW1;
            //         g(k,0,i).      fold[q] = g(k,1,i).        fold[q];      //FEQW1;//FEQW1;
            //         g(k,endpY,i).  gold[q] = g(k,endpY-1,i).  gold[q];//FEQW1;
            //         g(k,0,i).      gold[q] = g(k,1,i).        gold[q]; 

            // }
            /*Periodic boundary condition*/
            /*top*/
            // for mass and momentum distribution
            g(k,(g.NY - 1),i).  fold[4]   = g(k,1,i).  fold[4] ; 
            g(k,(g.NY - 1),i).  fold[8]   = g(k,1,i).  fold[8] ; 
            g(k,(g.NY - 1),i).  fold[12]  = g(k,1,i).  fold[12]; 
            g(k,(g.NY - 1),i).  fold[13]  = g(k,1,i).  fold[13];
            g(k,(g.NY - 1),i).  fold[18]  = g(k,1,i).  fold[18];
            g(k,(g.NY - 1),i).  fold[20]  = g(k,1,i).  fold[20];
            g(k,(g.NY - 1),i).  fold[22]  = g(k,1,i).  fold[22];
            g(k,(g.NY - 1),i).  fold[23]  = g(k,1,i).  fold[23];
            g(k,(g.NY - 1),i).  fold[26]  = g(k,1,i).  fold[26];
            // for energy distribution
            g(k,(g.NY - 1),i).  gold[4]   = g(k,1,i).  gold[4] ; 
            g(k,(g.NY - 1),i).  gold[8]   = g(k,1,i).  gold[8] ; 
            g(k,(g.NY - 1),i).  gold[12]  = g(k,1,i).  gold[12]; 
            g(k,(g.NY - 1),i).  gold[13]  = g(k,1,i).  gold[13];
            g(k,(g.NY - 1),i).  gold[18]  = g(k,1,i).  gold[18];
            g(k,(g.NY - 1),i).  gold[20]  = g(k,1,i).  gold[20];
            g(k,(g.NY - 1),i).  gold[22]  = g(k,1,i).  gold[22];
            g(k,(g.NY - 1),i).  gold[23]  = g(k,1,i).  gold[23];
            g(k,(g.NY - 1),i).  gold[26]  = g(k,1,i).  gold[26];
            /*bottom*/
            // for mass and momentum distribution
            g(k,0,i).    fold[3]  = g(k,(g.NY - 1)-1,i)  .fold[3] ;
            g(k,0,i).    fold[7]  = g(k,(g.NY - 1)-1,i)  .fold[7] ;
            g(k,0,i).    fold[11] = g(k,(g.NY - 1)-1,i)  .fold[11];
            g(k,0,i).    fold[14] = g(k,(g.NY - 1)-1,i)  .fold[14];
            g(k,0,i).    fold[17] = g(k,(g.NY - 1)-1,i)  .fold[17];
            g(k,0,i).    fold[19] = g(k,(g.NY - 1)-1,i)  .fold[19];
            g(k,0,i).    fold[21] = g(k,(g.NY - 1)-1,i)  .fold[21];
            g(k,0,i).    fold[24] = g(k,(g.NY - 1)-1,i)  .fold[24];
            g(k,0,i).    fold[25] = g(k,(g.NY - 1)-1,i)  .fold[25];
            // for energy distribution
            g(k,0,i).    gold[3]  = g(k,(g.NY - 1)-1,i)  .gold[3] ;
            g(k,0,i).    gold[7]  = g(k,(g.NY - 1)-1,i)  .gold[7] ;
            g(k,0,i).    gold[11] = g(k,(g.NY - 1)-1,i)  .gold[11];
            g(k,0,i).    gold[14] = g(k,(g.NY - 1)-1,i)  .gold[14];
            g(k,0,i).    gold[17] = g(k,(g.NY - 1)-1,i)  .gold[17];
            g(k,0,i).    gold[19] = g(k,(g.NY - 1)-1,i)  .gold[19];
            g(k,0,i).    gold[21] = g(k,(g.NY - 1)-1,i)  .gold[21];
            g(k,0,i).    gold[24] = g(k,(g.NY - 1)-1,i)  .gold[24];
            g(k,0,i).    gold[25] = g(k,(g.NY - 1)-1,i)  .gold[25];
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
}
/*front back boundary conditions*/
void boundaryConditiontWallFB(grid& g){
    // int endpX = g.NX - 1;
    // int endpY = g.NY - 1;
    // int endpZ = (g.NZ - 1);
    #pragma omp parallel for
    for(int j = 0; j < g.NY; ++j){
        for(int i = 0; i < g.NX; ++i){
            /*First order Neumann*/
            // for(int q = 0; q < 27; ++q){
            //         // double FEQW1     =  feqForBoundary(g(1,1).w[q], 1.0, 0.01, 0, g(1,1).cx[q], g(1,1).cy[q]);
            //         g(endpZ,j,i).  fold[q] = g(endpZ-1,j,i).  fold[q];//FEQW1;
            //         g(0,j,i).      fold[q] = g(1,j,i).        fold[q];      //FEQW1;//FEQW1;
            //         g(endpZ,j,i).  gold[q] = g(endpZ-1,j,i).  gold[q];//FEQW1;
            //         g(0,j,i).      gold[q] = g(1,j,i).        gold[q]; 

            // }
            /*Periodic boundary condition*/
            /*front*/
            // for mass and momentum distribution
            g((g.NZ - 1),j,i).  fold[6]   = g(1,j,i).  fold[6] ; 
            g((g.NZ - 1),j,i).  fold[10]  = g(1,j,i).  fold[10]; 
            g((g.NZ - 1),j,i).  fold[12]  = g(1,j,i).  fold[12]; 
            g((g.NZ - 1),j,i).  fold[15]  = g(1,j,i).  fold[15];
            g((g.NZ - 1),j,i).  fold[17]  = g(1,j,i).  fold[17];
            g((g.NZ - 1),j,i).  fold[20]  = g(1,j,i).  fold[20];
            g((g.NZ - 1),j,i).  fold[21]  = g(1,j,i).  fold[21];
            g((g.NZ - 1),j,i).  fold[24]  = g(1,j,i).  fold[24];
            g((g.NZ - 1),j,i).  fold[26]  = g(1,j,i).  fold[26];
            // for energy distribution
            g((g.NZ - 1),j,i).  gold[6]   = g(1,j,i).  gold[6] ; 
            g((g.NZ - 1),j,i).  gold[10]  = g(1,j,i).  gold[10]; 
            g((g.NZ - 1),j,i).  gold[12]  = g(1,j,i).  gold[12]; 
            g((g.NZ - 1),j,i).  gold[15]  = g(1,j,i).  gold[15];
            g((g.NZ - 1),j,i).  gold[17]  = g(1,j,i).  gold[17];
            g((g.NZ - 1),j,i).  gold[20]  = g(1,j,i).  gold[20];
            g((g.NZ - 1),j,i).  gold[21]  = g(1,j,i).  gold[21];
            g((g.NZ - 1),j,i).  gold[24]  = g(1,j,i).  gold[24];
            g((g.NZ - 1),j,i).  gold[26]  = g(1,j,i).  gold[26];
            /*back*/
            // for mass and momentum distribution
            g(0,j,i).    fold[5]  = g((g.NZ - 1)-1,j,i)  .fold[5] ;
            g(0,j,i).    fold[9]  = g((g.NZ - 1)-1,j,i)  .fold[9] ;
            g(0,j,i).    fold[11] = g((g.NZ - 1)-1,j,i)  .fold[11];
            g(0,j,i).    fold[16] = g((g.NZ - 1)-1,j,i)  .fold[16];
            g(0,j,i).    fold[18] = g((g.NZ - 1)-1,j,i)  .fold[18];
            g(0,j,i).    fold[19] = g((g.NZ - 1)-1,j,i)  .fold[19];
            g(0,j,i).    fold[22] = g((g.NZ - 1)-1,j,i)  .fold[22];
            g(0,j,i).    fold[23] = g((g.NZ - 1)-1,j,i)  .fold[23];
            g(0,j,i).    fold[25] = g((g.NZ - 1)-1,j,i)  .fold[25];
            // for energy distribution
            g(0,j,i).    gold[5]  = g((g.NZ - 1)-1,j,i)  .gold[5] ;
            g(0,j,i).    gold[9]  = g((g.NZ - 1)-1,j,i)  .gold[9] ;
            g(0,j,i).    gold[11] = g((g.NZ - 1)-1,j,i)  .gold[11];
            g(0,j,i).    gold[16] = g((g.NZ - 1)-1,j,i)  .gold[16];
            g(0,j,i).    gold[18] = g((g.NZ - 1)-1,j,i)  .gold[18];
            g(0,j,i).    gold[19] = g((g.NZ - 1)-1,j,i)  .gold[19];
            g(0,j,i).    gold[22] = g((g.NZ - 1)-1,j,i)  .gold[22];
            g(0,j,i).    gold[23] = g((g.NZ - 1)-1,j,i)  .gold[23];
            g(0,j,i).    gold[25] = g((g.NZ - 1)-1,j,i)  .gold[25];
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
}
/*Boundary conditions for corners in the domain
// void bcCorners(grid& g){
//     int endpX = g.NX - 1;
//     int endpY = g.NY - 1;
//     g(0,0)        .fold[5] = g(1,1)      .fold[7];
//     g(endpY,0)    .fold[8] = g(endpY-1,1).fold[6];
//     g(0,endpX)    .fold[6] = g(1,endpX-1).fold[8];
//     g(endpY,endpX).fold[7] = g(endpY-1,endpX-1).fold[5];

//     g(0,0)        .gold[5] = g(1,1)      .gold[7];
//     g(endpY,0)    .gold[8] = g(endpY-1,1).gold[6];
//     g(0,endpX)    .gold[6] = g(1,endpX-1).gold[8];
//     g(endpY,endpX).gold[7] = g(endpY-1,endpX-1).gold[5];
// }
*/
/*Knudsen number dependent stabilization*/
void KnudsenNumberStability(grid& g){
    #pragma omp parallel for
    for(int numberofNodes = 0 ; numberofNodes < g.NZ * g.NY * g.NX; ++numberofNodes){
        double diffSumf = 0.0;
        double diffSumg = 0.0;
        // double tauf     = 0.0;
        // double taug     = 0.0;
        if(g.domain[numberofNodes].flag == 0){
            for(int q=0; q < 27; ++q){
                diffSumf = diffSumf + abs(g.domain[numberofNodes].fold[q] - g.domain[numberofNodes].feq[q])/ g.domain[numberofNodes].feq[q];
                diffSumg = diffSumg + abs(g.domain[numberofNodes].gold[q] - g.domain[numberofNodes].geq[q])/ g.domain[numberofNodes].geq[q];
            }
            // tauf = 1.0 / g.domain[j].omegaf;
            // taug = 1.0 / g.domain[j].omegag;
            g.domain[numberofNodes].epsilonf = (1.0/27.0) * diffSumf;  //(1.0/grid[i].Q) //9 IS NUMBER OF VELOCOTY SETS
            g.domain[numberofNodes].epsilong = (1.0/27.0) * diffSumg;
            if(g.domain[numberofNodes].epsilonf < 0.01){
                g.domain[numberofNodes].alphaf = 1.0;
                g.domain[numberofNodes].alphag = 1.0;}
            else if((g.domain[numberofNodes].epsilonf >= 0.01) && (g.domain[numberofNodes].epsilonf < 0.1)){
                g.domain[numberofNodes].alphaf = 1.05;
                g.domain[numberofNodes].alphag = 1.05;}
            else if((g.domain[numberofNodes].epsilonf >= 0.1) && (g.domain[numberofNodes].epsilonf < 1.0)){
                g.domain[numberofNodes].alphaf = 1.35;
                g.domain[numberofNodes].alphag = 1.35;}
            else if(g.domain[numberofNodes].epsilonf >= 1.0){
                g.domain[numberofNodes].alphaf = 1.0/g.domain[numberofNodes].tauf;
                g.domain[numberofNodes].alphag = 1.0/g.domain[numberofNodes].taug;}
            g.domain[numberofNodes].tauf = g.domain[numberofNodes].tauf * g.domain[numberofNodes].alphaf;
            g.domain[numberofNodes].taug = g.domain[numberofNodes].taug * g.domain[numberofNodes].alphag;
            // std::cout <<  g.domain[j].tauf << "<>" <<  g.domain[j].taug <<std::endl; 
            g.domain[numberofNodes].omegaf = 1.0/g.domain[numberofNodes].tauf;
            g.domain[numberofNodes].omegag = 1.0/g.domain[numberofNodes].taug;        
        }
    }
}
/*swapping old and new populations*/
void swap(grid& g){
    #pragma omp parallel for
    for(int numberofNodes = 0; numberofNodes < g.NZ * g.NY * g.NX ; ++numberofNodes){
        for(int q = 0; q < 27; ++q){
            g.domain[numberofNodes].fold[q] = g.domain[numberofNodes].fnew[q];
            g.domain[numberofNodes].gold[q] = g.domain[numberofNodes].gnew[q];
        }
    } 
}

