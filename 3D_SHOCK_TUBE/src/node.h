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
    double uz;
    double T;
    double ushift;
    double epsilonf;
    double alphaf  ;
    double epsilong;
    double alphag  ;
    //distance for 3d case
    int x;
    int y;
    int z;
    //node velocities
    double cx[27];
    double cy[27];
    double cz[27];
    //node weights
    double  w[27];
    //relaxation factors dependon T
    double omegaf;
    double omegag;
    double tauf;
    double taug;
    double omegag2;
    // f population
    double fold[27];
    double fnew[27];
    double  feq[27];
    double  fshift[27];
    // g population
    double gold[27];
    double gnew[27];
    double  geq[27];
    double  gshift[27];
    double qstr[27];
    //pressure tensor since this is one dimensional we have 4 value per node
    double   P[9];
    double Peq[9];
    /*Lagrangian Multipliers considering ELBM minimization of ENTROPY function for 1D it's 2(dim+1 ;
      because of 2 constraints and one has to consider the macroscopic vewlocity
    */
    double LM[4];

    node(){}
    node(int dim, int qset)
        :D(dim),Q(qset){
            flag     = 0;
            nr       = 0;
            x        = 0;
            y        = 0;
            z        = 0;
            rho      = 0.0;
            ux       = 0.0;
            uy       = 0.0;
            uz       = 0.0;
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
            for(int i =0; i < 9; ++i){ 
                P[i]      = 0.0;
                Peq[i]    = 0.0;
            }
            //node speeds in standard node
            cx[0] = 0.0 ; cx[1] = 1.0; cx[2] = -1.0; cx[3] =  0.0; cx[4] =  0.0; cx[5] = 0.0; cx[6] =  0.0; cx[7] =  1.0; cx[8] = -1.0;
            cy[0] = 0.0 ; cy[1] = 0.0; cy[2] =  0.0; cy[3] =  1.0; cy[4] = -1.0; cy[5] = 0.0; cy[6] =  0.0; cy[7] =  1.0; cy[8] = -1.0; 
            cz[0] = 0.0 ; cz[1] = 0.0; cz[2] =  0.0; cz[3] =  0.0; cz[4] =  0.0; cz[5] = 1.0; cz[6] = -1.0; cz[7] =  0.0; cz[8] =  0.0;
// 
            cx[9] = 1.0 ; cx[10] =-1.0; cx[11] =  0.0; cx[12] =  0.0; cx[13] =  1.0; cx[14] =-1.0; cx[15] =  1.0; cx[16] = -1.0; cx[17] =  0.0;
            cy[9] = 0.0 ; cy[10] = 0.0; cy[11] =  1.0; cy[12] = -1.0; cy[13] = -1.0; cy[14] = 1.0; cy[15] =  0.0; cy[16] =  0.0; cy[17] =  1.0; 
            cz[9] = 1.0 ; cz[10] =-1.0; cz[11] =  1.0; cz[12] = -1.0; cz[13] =  0.0; cz[14] = 0.0; cz[15] = -1.0; cz[16] =  1.0; cz[17] = -1.0;
// 
            cx[18] = 0.0 ; cx[19] = 1.0; cx[20] = -1.0; cx[21] =  1.0; cx[22] = -1.0; cx[23] = 1.0; cx[24] = -1.0; cx[25] = -1.0; cx[26] =  1.0;
            cy[18] =-1.0 ; cy[19] = 1.0; cy[20] = -1.0; cy[21] =  1.0; cy[22] = -1.0; cy[23] =-1.0; cy[24] =  1.0; cy[25] =  1.0; cy[26] = -1.0; 
            cz[18] = 1.0 ; cz[19] = 1.0; cz[20] = -1.0; cz[21] = -1.0; cz[22] =  1.0; cz[23] = 1.0; cz[24] = -1.0; cz[25] =  1.0; cz[26] = -1.0;

            for(int i =0; i < 27; ++i){
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
            for (int i =0 ; i < 4; ++i){
                LM[i] = 0.5;
            }
    }
    node(const node& src)
        :D(src.D),Q(src.Q){
            flag    = src.flag;
            nr      = src.nr;
            x       = src.x;
            y       = src.y;
            z       = src.z;
            rho     = src.rho     ;
            ux       = src.ux     ;
            uy       = src.uy     ;
            uz       = src.uz     ;
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
            for(int i =0; i < 9; ++i){
                P[i]       = src.P[i]     ;
                Peq[i]     = src.Peq[i]    ;
            }

            for(int i =0; i < 27; ++i){
                cx[i]    = src.cx[i];
                cy[i]    = src.cy[i];
                cz[i]    = src.cz[i];
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
            for(int i = 0; i < 4;++i){
                LM[i]  = src.LM[i];
            }
        }
    node& operator=(const node& src){
        if(this == &src)
            return *this;

        D       = src.D;
        Q       = src.Q;
        nr      = src.nr;
        flag    = src.flag;
        x       = src.x;
        y       = src.y;
        z       = src.z;
        rho     = src.rho     ;
        ux      = src.ux     ;
        uy      = src.uy     ;
        uz      = src.uz     ;
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
        for(int i =0; i < 9; ++i){ 
            P[i]      = src.P[i]      ;
            Peq[i]    = src.Peq[i]   ;
        }
        for(int i =0; i < 27; ++i){
            cx[i]    = src.cx[i];
            cy[i]    = src.cy[i];
            cz[i]    = src.cz[i];
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
        for(int i = 0; i < 4;++i){
            LM[i]  = src.LM[i];
        }

        return *this;
    }

    // node &operator()(int i, int j){
    //     return gridPoint[gridSize_X * j + i];
    // }

    ~node(){}
};

