#pragma once

#include <iostream>
#include "node.h"

class grid{

public: 
    int NX;
    int NY;
    int NZ;
    node* domain;//[1500*5]; //relates to nx FIG3: 3000

    grid(){}
    grid(int nx,int ny,int nz):NX(nx),NY(ny),NZ(nz){
        domain = nullptr;
        domain = (node*) malloc(sizeof(node)*(NX*NY*NZ));
        for(int i =0; i < NX*NY*NZ; ++i){
            domain[i].alphaf   = 0.0;
            domain[i].nr       = 0.0;
            domain[i].alphag   = 0.0;
            domain[i].epsilonf = 0.0;
            domain[i].epsilong = 0.0;
            domain[i].omegaf   = 0.0;
            domain[i].omegag   = 0.0;
            domain[i].tauf     = 0.0;
            domain[i].taug     = 0.0;
            domain[i].omegag2  = 0.0;
            domain[i].rho      = 0.0;
            domain[i].T        = 0.0;
            domain[i].ushift   = 0.0;
            domain[i].ux       = 0.0;
            domain[i].uy       = 0.0;
            domain[i].uz       = 0.0;
            for(int k = 0; k < 9; ++k){ //in 3d the pressure tensor is 3x3 = 9 elements
                domain[i].P[k]        = 0.0;
                domain[i].Peq[k]      = 0.0;
            }
            domain[i].D        = 3;
            domain[i].Q        = 27;
            domain[i].flag     = 0;
            domain[i].x        = 0;
            domain[i].y        = 0;
            domain[i].z        = 0;

            domain[i].cx[0] = 0.0 ; domain[i].cx[1] = 1.0; domain[i].cx[2] = -1.0; domain[i].cx[3] =  0.0; domain[i].cx[4] =  0.0; domain[i].cx[5] = 0.0; domain[i].cx[6] =  0.0; domain[i].cx[7] =  1.0; domain[i].cx[8] = -1.0;
            domain[i].cy[0] = 0.0 ; domain[i].cy[1] = 0.0; domain[i].cy[2] =  0.0; domain[i].cy[3] =  1.0; domain[i].cy[4] = -1.0; domain[i].cy[5] = 0.0; domain[i].cy[6] =  0.0; domain[i].cy[7] =  1.0; domain[i].cy[8] = -1.0; 
            domain[i].cz[0] = 0.0 ; domain[i].cz[1] = 0.0; domain[i].cz[2] =  0.0; domain[i].cz[3] =  0.0; domain[i].cz[4] =  0.0; domain[i].cz[5] = 1.0; domain[i].cz[6] = -1.0; domain[i].cz[7] =  0.0; domain[i].cz[8] =  0.0;
// 
            domain[i].cx[9] = 1.0 ; domain[i].cx[10] =-1.0; domain[i].cx[11] =  0.0; domain[i].cx[12] =  0.0; domain[i].cx[13] =  1.0; domain[i].cx[14] =-1.0; domain[i].cx[15] =  1.0; domain[i].cx[16] = -1.0; domain[i].cx[17] =  0.0;
            domain[i].cy[9] = 0.0 ; domain[i].cy[10] = 0.0; domain[i].cy[11] =  1.0; domain[i].cy[12] = -1.0; domain[i].cy[13] = -1.0; domain[i].cy[14] = 1.0; domain[i].cy[15] =  0.0; domain[i].cy[16] =  0.0; domain[i].cy[17] =  1.0; 
            domain[i].cz[9] = 1.0 ; domain[i].cz[10] =-1.0; domain[i].cz[11] =  1.0; domain[i].cz[12] = -1.0; domain[i].cz[13] =  0.0; domain[i].cz[14] = 0.0; domain[i].cz[15] = -1.0; domain[i].cz[16] =  1.0; domain[i].cz[17] = -1.0;
// 
            domain[i].cx[18] = 0.0 ; domain[i].cx[19] = 1.0; domain[i].cx[20] = -1.0; domain[i].cx[21] =  1.0; domain[i].cx[22] = -1.0; domain[i].cx[23] = 1.0; domain[i].cx[24] = -1.0; domain[i].cx[25] = -1.0; domain[i].cx[26] =  1.0;
            domain[i].cy[18] =-1.0 ; domain[i].cy[19] = 1.0; domain[i].cy[20] = -1.0; domain[i].cy[21] =  1.0; domain[i].cy[22] = -1.0; domain[i].cy[23] =-1.0; domain[i].cy[24] =  1.0; domain[i].cy[25] =  1.0; domain[i].cy[26] = -1.0; 
            domain[i].cz[18] = 1.0 ; domain[i].cz[19] = 1.0; domain[i].cz[20] = -1.0; domain[i].cz[21] = -1.0; domain[i].cz[22] =  1.0; domain[i].cz[23] = 1.0; domain[i].cz[24] = -1.0; domain[i].cz[25] =  1.0; domain[i].cz[26] = -1.0;

            domain[i].LM[0] = 0.5;
            domain[i].LM[1] = 0.5;
            domain[i].LM[2] = 0.5;
            domain[i].LM[3] = 0.5;

            for(int j = 0; j < 27; ++j){
                domain[i]. feq[j]   = 0.0;
                domain[i].fold[j]   = 0.0;
                domain[i].fnew[j]   = 0.0;
                domain[i]. geq[j]   = 0.0;
                domain[i].gold[j]   = 0.0;
                domain[i].gnew[j]   = 0.0;
                domain[i].qstr[j]   = 0.0;
                domain[i].   w[j]   = 0.0;
                domain[i].fshift[j] = 0.0;
                domain[i].gshift[j] = 0.0;
            }
        }
 


    }
    grid(const grid& src)
    :NX(src.NX),NY(src.NY),NZ(src.NZ){
        domain = nullptr;
        domain = (node*)malloc(sizeof(node)*(NX*NY*NZ));
        for(int i =0; i < NX*NY*NZ; ++i){
            domain[i].nr       = src.domain[i].nr   ;
            domain[i].alphaf   = src.domain[i].alphaf   ;
            domain[i].alphag   = src.domain[i].alphag   ;
            domain[i].epsilonf = src.domain[i].epsilonf ;
            domain[i].epsilong = src.domain[i].epsilong ;
            domain[i].omegaf   = src.domain[i].omegaf   ;
            domain[i].omegag   = src.domain[i].omegag   ;
            domain[i].tauf     = src.domain[i].tauf     ;
            domain[i].taug     = src.domain[i].taug     ;
            domain[i].omegag2  = src.domain[i].omegag2  ;
            domain[i].rho      = src.domain[i].rho      ;
            domain[i].T        = src.domain[i].T        ;
            domain[i].ushift   = src.domain[i].ushift   ;
            domain[i].ux       = src.domain[i].ux       ;
            domain[i].uy       = src.domain[i].uy       ;
            domain[i].uz       = src.domain[i].uz       ;
            for(int k =0; k < 9; ++k){
                domain[i].P[k]        = src.domain[i].P[k]        ;
                domain[i].Peq[k]      = src.domain[i].Peq[k]      ;
            }
            domain[i].D        = src.domain[i].D        ;
            domain[i].Q        = src.domain[i].Q        ;
            domain[i].flag     = src.domain[i].flag     ;
            domain[i].x        = src.domain[i].x;
            domain[i].y        = src.domain[i].y;
            domain[i].z        = src.domain[i].z;            
            // domain[i].cx[0]    = src.domain[i].cx[0];
            // domain[i].cx[1]    = src.domain[i].cx[1];
            // domain[i].cx[2]    = src.domain[i].cx[2];

            for(int j = 0; j < 27; ++j){
                domain[i].  cx[j]   = src.domain[i].  cx[j];
                domain[i].  cy[j]   = src.domain[i].  cy[j];
                domain[i]. feq[j]   = src.domain[i]. feq[j];
                domain[i].fold[j]   = src.domain[i].fold[j];
                domain[i].fnew[j]   = src.domain[i].fnew[j];
                domain[i]. geq[j]   = src.domain[i]. geq[j];
                domain[i].gold[j]   = src.domain[i].gold[j];
                domain[i].gnew[j]   = src.domain[i].gnew[j];
                domain[i].qstr[j]   = src.domain[i].qstr[j];
                domain[i].   w[j]   = src.domain[i].   w[j];
                domain[i].fshift[j] = src.domain[i].fshift[j];
                domain[i].gshift[j] = src.domain[i].gshift[j];
            }

            domain[i].LM[0] = src.domain[i].LM[0];
            domain[i].LM[1] = src.domain[i].LM[1];
            domain[i].LM[2] = src.domain[i].LM[2];
            domain[i].LM[3] = src.domain[i].LM[3];
        }
 
    }
    grid& operator=(const grid& src){
        if(this == &src){
            return *this;
        }
        NX = src.NX;
        NY = src.NY;
        NZ = src.NZ;
        domain = nullptr;
        domain = (node*)malloc(sizeof(node)*(NX*NY*NZ));
        for(int i =0; i < NX*NY*NZ; ++i){
            domain[i].nr       = src.domain[i].nr   ;
            domain[i].alphaf   = src.domain[i].alphaf   ;
            domain[i].alphag   = src.domain[i].alphag   ;
            domain[i].epsilonf = src.domain[i].epsilonf ;
            domain[i].epsilong = src.domain[i].epsilong ;
            domain[i].omegaf   = src.domain[i].omegaf   ;
            domain[i].omegag   = src.domain[i].omegag   ;
            domain[i].tauf     = src.domain[i].tauf     ;
            domain[i].taug     = src.domain[i].taug     ;
            domain[i].omegag2  = src.domain[i].omegag2  ;
            domain[i].rho      = src.domain[i].rho      ;
            domain[i].T        = src.domain[i].T        ;
            domain[i].ushift   = src.domain[i].ushift   ;
            domain[i].ux       = src.domain[i].ux       ;
            domain[i].uy       = src.domain[i].uy       ;
            domain[i].uz       = src.domain[i].uz       ;
            for(int k =0; k < 9; ++k){
                domain[i].P[k]        = src.domain[i].P[k]        ;
                domain[i].Peq[k]      = src.domain[i].Peq[k]      ;
            }
            domain[i].D        = src.domain[i].D        ;
            domain[i].Q        = src.domain[i].Q        ;
            domain[i].flag     = src.domain[i].flag     ;
            domain[i].x        = src.domain[i].x;
            domain[i].y        = src.domain[i].y;
            domain[i].z        = src.domain[i].z;
            // domain[i].cx[0]    = src.domain[i].cx[0];
            // domain[i].cx[1]    = src.domain[i].cx[1];
            // domain[i].cx[2]    = src.domain[i].cx[2];

            for(int j = 0; j < 27; ++j){
                domain[i].  cx[j]   = src.domain[i].  cx[j];
                domain[i].  cy[j]   = src.domain[i].  cy[j];
                domain[i]. feq[j]   = src.domain[i]. feq[j];
                domain[i].fold[j]   = src.domain[i].fold[j];
                domain[i].fnew[j]   = src.domain[i].fnew[j];
                domain[i]. geq[j]   = src.domain[i]. geq[j];
                domain[i].gold[j]   = src.domain[i].gold[j];
                domain[i].gnew[j]   = src.domain[i].gnew[j];
                domain[i].qstr[j]   = src.domain[i].qstr[j];
                domain[i].   w[j]   = src.domain[i].   w[j];
                domain[i].fshift[j] = src.domain[i].fshift[j];
                domain[i].gshift[j] = src.domain[i].gshift[j];
            }

            domain[i].LM[0] = src.domain[i].LM[0];
            domain[i].LM[1] = src.domain[i].LM[1];
            domain[i].LM[2] = src.domain[i].LM[2];
            domain[i].LM[3] = src.domain[i].LM[3];
        }

        return *this; 
    }

    node &operator()(int k, int j, int i){
        return domain[ NX * NY * k + NX * j + i];
    }

    ~grid(){
        free(domain);
    }

};