#pragma once

#include <iostream>

class node{

public:
    // latice parameters
    // int NX;
    int D;
    int Q;
    int nr;
    //macroscopic parameters
    int flag;
    double rho;
    double ux;
    double uy;
    double T;
    double ushift;
    double epsilonf;
    double alphaf  ;
    double epsilong;
    double alphag  ;
    //node velocities
    double cx[9];
    double cy[9];
    //node weights
    double  w[9];
    //relaxation factors dependon T
    double omegaf;
    double omegag;
    double tauf;
    double taug;
    double omegag2;
    // f population
    double fold[9];
    double fnew[9];
    double  feq[9];
    double  fshift[9];
    // g population
    double gold[9];
    double gnew[9];
    double  geq[9];
    double  gshift[9];
    double qstr[9];
    //pressure tensor since this is one dimensional we have 4 value per node
    double   P[4];
    double Peq[4];
    /*Lagrangian Multipliers considering ELBM minimization of ENTROPY function for 1D it's 2(dim+1 ;
      because of 2 constraints and one has to consider the macroscopic vewlocity
    */
    double LM[3];

    node(){}
    node(int dim, int qset)
        :D(dim),Q(qset){
            nr       = 0;
            flag     = 0;
            rho      = 0.0;
            ux       = 0.0;
            uy       = 0.0;
            T        = 0.0;
            ushift   = 0.0; 
            omegaf   = 0.0;
            omegag   = 0.0;
            tauf = 0.0;
            taug = 0.0;

            omegag2  = 0.0;
            epsilonf = 0.0;
            alphaf   = 0.0;
            epsilong = 0.0;
            alphag   = 0.0;
            for(int i =0; i < 4; ++i){ 
                P[i]      = 0.0;
                Peq[i]    = 0.0;
            }
            //node speeds in standard node
            cx[0] = 0.0 ; cx[1] = 1.0; cx[2] = 0.0; cx[3]= -1.0; cx[4]=  0.0; cx[5] = 1.0; cx[6] = -1.0; cx[7] = -1.0; cx[8]=  1.0;
            cy[0] = 0.0 ; cy[1] = 0.0; cy[2] = 1.0; cy[3]=  0.0; cy[4]= -1.0; cy[5] = 1.0; cy[6] =  1.0; cy[7] = -1.0; cy[8]= -1.0; 
            
            for(int i =0; i < Q; ++i){
                w[i]     = 0.0; 
                fold [i] = 0.0;
                fnew [i] = 0.0;
                feq  [i] = 0.0; 
                gold [i] = 0.0;
                gnew [i] = 0.0;
                geq  [i] = 0.0;
                qstr [i] = 0.0; 
                fshift[i] = 0.0;
                gshift[i] = 0.0;  
            }
            /*3 lagranginans*/
            for (int i =0 ; i < 3; ++i){
                LM[i] = 0.5;
            }
    }
    node(const node& src)
        :D(src.D),Q(src.Q){
            nr      = src.nr;
            flag    = src.flag;
            rho     = src.rho     ;
            ux       = src.ux     ;
            uy       = src.uy     ;
            T       = src.T       ;
            ushift  = src.ushift  ;
            epsilonf= src.epsilonf; 
            alphaf  = src.alphaf  ; 
            epsilong= src.epsilong; 
            alphag  = src.alphag  ; 
            omegaf  = src.omegaf  ;
            tauf = src.tauf;
            taug = src.taug;
            omegag  = src.omegag  ;
            omegag2 = src.omegag2 ; 
            for(int i =0; i < 4; ++i){
                P[i]       = src.P[i]     ;
                Peq[i]     = src.Peq[i]    ;
            }

            for(int i =0; i < 9; ++i){
                cx[i]    = src.cx[i];
                cy[i]    = src.cy[i];
                w[i]     = src.w[i];    
                fold [i] = src.fold [i];
                fnew [i] = src.fnew [i];
                feq  [i] = src.feq  [i];
                gold [i] = src.gold [i];
                gnew [i] = src.gnew [i];
                geq  [i] = src.geq  [i];
                qstr [i] = src.qstr [i];
                fshift[i] = src.fshift[i];
                gshift[i] = src.gshift[i];
            }
            for(int i = 0; i < 3;++i){
                LM[i]  = src.LM[i];
            }
        }
    node& operator=(const node& src){
        if(this == &src)
            return *this;
        nr      = src.nr;
        D       = src.D;
        Q       = src.Q;
        flag    = src.flag;
        rho     = src.rho     ;
        ux      = src.ux     ;
        uy      = src.uy     ;
        T       = src.T       ;
        epsilonf= src.epsilonf;
        alphaf  = src.alphaf  ;
        epsilong= src.epsilong;
        alphag  = src.alphag  ;
        ushift  = src.ushift  ; 
        omegaf  = src.omegaf  ;
        omegag  = src.omegag  ;
        tauf = src.tauf;
        taug = src.taug;
        omegag2 = src.omegag2 ;
        for(int i =0; i < 4; ++i){ 
            P[i]      = src.P[i]      ;
            Peq[i]    = src.Peq[i]   ;
        }
        for(int i =0; i < 9; ++i){
            cx[i]    = src.cx[i];
            cy[i]    = src.cy[i];
            w[i]     = src.w[i];    
            fold [i] = src.fold [i];
            fnew [i] = src.fnew [i];
            feq  [i] = src.feq  [i];
            gold [i] = src.gold [i];
            gnew [i] = src.gnew [i];
            geq  [i] = src.geq  [i];
            qstr [i] = src.qstr [i];
            fshift[i] = src.fshift[i];
            gshift[i] = src.gshift[i];
        }
        for(int i = 0; i < 3;++i){
            LM[i]  = src.LM[i];
        }

        return *this;
    }

    // node &operator()(int i, int j){
    //     return gridPoint[gridSize_X * j + i];
    // }

    ~node(){}
};

