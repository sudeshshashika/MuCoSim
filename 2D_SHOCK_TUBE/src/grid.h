#pragma once

#include <iostream>
#include "node.h"

class grid{

public: 
    int NX;
    int NY;
    node* domain;//[1500*5]; //relates to nx FIG3: 3000

    grid(){}
    grid(int nx,int ny):NX(nx),NY(ny){
        domain = nullptr;
        domain = (node*) malloc(sizeof(node)*(NX*NY));
        for(int i =0; i < NX*NY; ++i){
            domain[i].nr       = 0.0;
            domain[i].alphaf   = 0.0;
            domain[i].alphag   = 0.0;
            domain[i].epsilonf = 0.0;
            domain[i].epsilong = 0.0;
            domain[i].omegaf   = 0.0;
            domain[i].omegag   = 0.0;
            domain[i].tauf   = 0.0;
            domain[i].taug   = 0.0;
            domain[i].omegag2  = 0.0;
            domain[i].rho      = 0.0;
            domain[i].T        = 0.0;
            domain[i].ushift   = 0.0;
            domain[i].ux       = 0.0;
            domain[i].uy       = 0.0;
            for(int k = 0; k < 4; ++k){
                domain[i].P[k]        = 0.0;
                domain[i].Peq[k]      = 0.0;
            }
            domain[i].D        = 2;
            domain[i].Q        = 9;
            domain[i].flag     = 0;

            domain[i].cx[0] = 0.0; domain[i].cx[1] = 1.0; domain[i].cx[2] = 0.0; domain[i].cx[3]= -1.0; domain[i].cx[4]=  0.0; 
            domain[i].cy[0] = 0.0; domain[i].cy[1] = 0.0; domain[i].cy[2] = 1.0; domain[i].cy[3]=  0.0; domain[i].cy[4]= -1.0; 

            domain[i].cx[5] = 1.0; domain[i].cx[6] = -1.0; domain[i].cx[7] = -1.0; domain[i].cx[8]=  1.0;
            domain[i].cy[5] = 1.0; domain[i].cy[6] =  1.0; domain[i].cy[7] = -1.0; domain[i].cy[8]= -1.0;

            domain[i].LM[0] = 0.5;
            domain[i].LM[1] = 0.5;
            domain[i].LM[2] = 0.5;
            for(int j = 0; j < 9; ++j){
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
    :NX(src.NX),NY(src.NY){
        domain = nullptr;
        domain = (node*)malloc(sizeof(node)*(NX*NY));
        for(int i =0; i < NX*NY; ++i){
            domain[i].nr       = src.domain[i].nr   ;
            domain[i].alphaf   = src.domain[i].alphaf   ;
            domain[i].alphag   = src.domain[i].alphag   ;
            domain[i].epsilonf = src.domain[i].epsilonf ;
            domain[i].epsilong = src.domain[i].epsilong ;
            domain[i].omegaf   = src.domain[i].omegaf   ;
            domain[i].omegag   = src.domain[i].omegag   ;
            domain[i].tauf   = src.domain[i].tauf   ;
            domain[i].taug   = src.domain[i].taug   ;
            domain[i].omegag2  = src.domain[i].omegag2  ;
            domain[i].rho      = src.domain[i].rho      ;
            domain[i].T        = src.domain[i].T        ;
            domain[i].ushift   = src.domain[i].ushift   ;
            domain[i].ux       = src.domain[i].ux       ;
            domain[i].uy       = src.domain[i].uy       ;
            for(int k =0; k < 4; ++k){
                domain[i].P[k]        = src.domain[i].P[k]        ;
                domain[i].Peq[k]      = src.domain[i].Peq[k]      ;
            }
            domain[i].D        = src.domain[i].D        ;
            domain[i].Q        = src.domain[i].Q        ;
            domain[i].flag     = src.domain[i].flag     ;
            // domain[i].cx[0]    = src.domain[i].cx[0];
            // domain[i].cx[1]    = src.domain[i].cx[1];
            // domain[i].cx[2]    = src.domain[i].cx[2];

            for(int j = 0; j < 9; ++j){
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
        }
 
    }
    grid& operator=(const grid& src){
        if(this == &src){
            return *this;
        }
        NX = src.NX;
        NY = src.NY;
        domain = nullptr;
        domain = (node*)malloc(sizeof(node)*(NX*NY));
        for(int i =0; i < NX*NY; ++i){
            domain[i].nr       = src.domain[i].nr   ;
            domain[i].alphaf   = src.domain[i].alphaf   ;
            domain[i].alphag   = src.domain[i].alphag   ;
            domain[i].epsilonf = src.domain[i].epsilonf ;
            domain[i].epsilong = src.domain[i].epsilong ;
            domain[i].omegaf   = src.domain[i].omegaf   ;
            domain[i].omegag   = src.domain[i].omegag   ;
            domain[i].tauf   = src.domain[i].tauf   ;
            domain[i].taug   = src.domain[i].taug   ;
            domain[i].omegag2  = src.domain[i].omegag2  ;
            domain[i].rho      = src.domain[i].rho      ;
            domain[i].T        = src.domain[i].T        ;
            domain[i].ushift   = src.domain[i].ushift   ;
            domain[i].ux       = src.domain[i].ux       ;
            domain[i].uy       = src.domain[i].uy       ;
            for(int k =0; k < 4; ++k){
                domain[i].P[k]        = src.domain[i].P[k]        ;
                domain[i].Peq[k]      = src.domain[i].Peq[k]      ;
            }
            domain[i].D        = src.domain[i].D        ;
            domain[i].Q        = src.domain[i].Q        ;
            domain[i].flag     = src.domain[i].flag     ;
            // domain[i].cx[0]    = src.domain[i].cx[0];
            // domain[i].cx[1]    = src.domain[i].cx[1];
            // domain[i].cx[2]    = src.domain[i].cx[2];

            for(int j = 0; j < 9; ++j){
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
        }

        return *this; 
    }

    node &operator()(int j, int i){
        return domain[NX * j + i];
    }

    ~grid(){
        free(domain);
    }

};