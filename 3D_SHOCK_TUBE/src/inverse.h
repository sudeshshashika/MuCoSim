#pragma once

#include <iostream>

class Inverse{
    friend std::ostream &operator << (std::ostream &os, const Inverse &rhs){
        // for(int i =0; i<rhs.N;++i){
        os << rhs.INV[0] << "          " <<  rhs.INV[1] << "          " <<  rhs.INV[2] << "\n";
        os << rhs.INV[3] << "          " <<  rhs.INV[4] << "          " <<  rhs.INV[5] << "\n";
        os << rhs.INV[6] << "          " <<  rhs.INV[7] << "          " <<  rhs.INV[8] << "\n";
        return os;
        // }
    }

public:
    int N; //create with total number of elements in 3*3 matrix which is 9
    double DET;
    double ADJ[9];
    double INV[9];
    double  TR[9];
    double   J[9];
    double   I[9];

    // double ADJ[4];
    // double INV[4];
    // double  TR[4];
    // double   J[4];
    // double   I[4];

    Inverse(){}
    Inverse(double N)
        :N(N){
            DET = 0.0;
            for(int i =0; i<N;++i){
                ADJ[i] = 0.0;
                INV[i] = 0.0;
                TR[i]  = 0.0;
                J[i]   = 0.0;
                I[i]   = 0.0;
            }
            // I[0] = 1.0; I[4] = 1.0; I[8] = 1.0;
        }

    Inverse(const Inverse& src)
        :N(src.N){
            DET = src.DET;
            for(int i =0; i<N;++i){
                ADJ[i] = src.ADJ[i];
                INV[i] = src.INV[i];
                TR[i]  = src.TR[i];
                J[i]   = src.J[i];
                I[i]   = src.I[i];
            }
    }

    Inverse& operator=(const Inverse& src){
        if(this == &src){
            return *this;
        }
        DET = src.DET;
        for(int i =0; i < N; ++i){
            ADJ[i] = src.ADJ[i];
            INV[i] = src.INV[i];
            TR[i]  = src.TR[i];
            J[i]   = src.J[i];
            I[i]   = src.I[i];
        }
    }

    void assign(double* array){
        for(int i =0; i < N;++i){
            J[i] = array[i];
        }
    }

    void adjoint(){
        ADJ[0] =  1.0 * (J[4]*J[8] - J[5]*J[7]);
        ADJ[1] = -1.0 * (J[3]*J[8] - J[5]*J[6]);
        ADJ[2] =  1.0 * (J[3]*J[7] - J[4]*J[6]);
        ADJ[3] = -1.0 * (J[1]*J[8] - J[2]*J[7]);
        ADJ[4] =  1.0 * (J[0]*J[8] - J[2]*J[6]);
        ADJ[5] = -1.0 * (J[0]*J[7] - J[1]*J[6]);
        ADJ[6] =  1.0 * (J[1]*J[5] - J[2]*J[4]);
        ADJ[7] = -1.0 * (J[0]*J[5] - J[2]*J[3]);
        ADJ[8] =  1.0 * (J[0]*J[4] - J[1]*J[3]);
    }

    void Det(){
        /*3d*/
        DET =    J[0]*(J[4]*J[8] - J[5]*J[7]) -
                 J[1]*(J[3]*J[8] - J[5]*J[6]) +
                 J[2]*(J[3]*J[7] - J[4]*J[6]);
        /*2d*/
        // DET =   (J[0]*J[3] - J[1]*J[2]);
    } 


    void transpose(){
        /*3d*/
        for(int i=0; i < N; ++i){
            TR[i] = ADJ[i];
        }
        TR[1] = ADJ[3];
        TR[3] = ADJ[1];
        TR[2] = ADJ[6];
        TR[6] = ADJ[2];
        TR[5] = ADJ[7];
        TR[7] = ADJ[5];
        /*2d*/
        // for(int i=0; i < N; ++i){
        //     TR[i] = J[i];
        // }
        // TR[0] = J[3];
        // TR[1] = -1.0*J[1];
        // TR[2] = -1.0*J[2];
        // TR[3] = J[0];

   
    }

    void Inv(){
        adjoint(); //for3D uncomment and  for 2D comment
        transpose();

        Det();
        for(int i=0; i < N; ++i){
            INV[i] = TR[i] / DET;
        }
    }

    ~Inverse(){
    }
};



// class Inverse{
//     friend std::ostream &operator << (std::ostream &os, const Inverse &rhs){
//         // for(int i =0; i<rhs.N;++i){
//         os << rhs.INV[0] << "          " <<  rhs.INV[1] << "          " <<  rhs.INV[2] << "\n";
//         os << rhs.INV[3] << "          " <<  rhs.INV[4] << "          " <<  rhs.INV[5] << "\n";
//         os << rhs.INV[6] << "          " <<  rhs.INV[7] << "          " <<  rhs.INV[8] << "\n";
//         return os;
//         // }
//     }

// public:
//     int N;
//     double DET;
//     double* ADJ;
//     double* INV;
//     double* TR;
//     double* J;
//     double* I;

//     // double ADJ[4];
//     // double INV[4];
//     // double  TR[4];
//     // double   J[4];
//     // double   I[4];

//     Inverse(){}
//     Inverse(double N)
//         :N(N){
//             DET = 0.0;
//             ADJ = (double*)malloc(sizeof(double) * (N));
//             INV = (double*)malloc(sizeof(double) * (N));
//             TR  = (double*)malloc(sizeof(double) * (N));
//             J   = (double*)malloc(sizeof(double) * (N));
//             I   = (double*)malloc(sizeof(double) * (N));
//             for(int i =0; i<N;++i){
//                 ADJ[i] = 0.0;
//                 INV[i] = 0.0;
//                 TR[i]  = 0.0;
//                 J[i]   = 0.0;
//                 I[i]   = 0.0;
//             }
//             // I[0] = 1.0; I[4] = 1.0; I[8] = 1.0;
//         }

//     Inverse(const Inverse& src)
//         :N(src.N){
//             DET = src.DET;
//             ADJ = nullptr;
//             INV = nullptr;
//             J   = nullptr;
//             I   = nullptr;
//             ADJ = (double*)malloc(sizeof(double) * (N));
//             INV = (double*)malloc(sizeof(double) * (N));
//             TR  = (double*)malloc(sizeof(double) * (N));
//             J   = (double*)malloc(sizeof(double) * (N));
//             I   = (double*)malloc(sizeof(double) * (N));
//             for(int i =0; i<N;++i){
//                 ADJ[i] = src.ADJ[i];
//                 INV[i] = src.INV[i];
//                 TR[i]  = src.TR[i];
//                 J[i]   = src.J[i];
//                 I[i]   = src.I[i];
//             }
//     }

//     Inverse& operator=(const Inverse& src){
//         if(this == &src){
//             return *this;
//         }
//         DET = src.DET;
//         ADJ = nullptr;
//         INV = nullptr;
//         TR  = nullptr;
//         J   = nullptr;
//         I   = nullptr;
//         ADJ = (double*)malloc(sizeof(double) * (N));
//         INV = (double*)malloc(sizeof(double) * (N));
//         TR  = (double*)malloc(sizeof(double) * (N));
//         J   = (double*)malloc(sizeof(double) * (N));
//         I   = (double*)malloc(sizeof(double) * (N));
//         for(int i =0; i<N;++i){
//             ADJ[i] = src.ADJ[i];
//             INV[i] = src.INV[i];
//             TR[i]  = src.TR[i];
//             J[i]   = src.J[i];
//             I[i]   = src.I[i];
//         }
//     }

//     void assign(double* array){
//         for(int i =0; i<N;++i){
//             J[i] = array[i];
//         }
//     }

//     void adjoint(){
//         ADJ[0] =  1.0 * (J[4]*J[8] - J[5]*J[7]);
//         ADJ[1] = -1.0 * (J[3]*J[8] - J[5]*J[6]);
//         ADJ[2] =  1.0 * (J[3]*J[7] - J[4]*J[6]);
//         ADJ[3] = -1.0 * (J[1]*J[8] - J[2]*J[7]);
//         ADJ[4] =  1.0 * (J[0]*J[8] - J[2]*J[6]);
//         ADJ[5] = -1.0 * (J[0]*J[7] - J[1]*J[6]);
//         ADJ[6] =  1.0 * (J[1]*J[5] - J[2]*J[4]);
//         ADJ[7] = -1.0 * (J[0]*J[5] - J[2]*J[3]);
//         ADJ[8] =  1.0 * (J[0]*J[4] - J[1]*J[3]);
//     }

//     void Det(){
//         /*3d*/
//         DET =    J[0]*(J[4]*J[8] - J[5]*J[7]) -
//                  J[1]*(J[3]*J[8] - J[5]*J[6]) +
//                  J[2]*(J[3]*J[7] - J[4]*J[6]);
//         /*2d*/
//         // DET =   (J[0]*J[3] - J[1]*J[2]);
//     } 


//     void transpose(){
//         /*3d*/
//         for(int i=0; i < N; ++i){
//             TR[i] = ADJ[i];
//         }
//         TR[1] = ADJ[3];
//         TR[3] = ADJ[1];
//         TR[2] = ADJ[6];
//         TR[6] = ADJ[2];
//         TR[5] = ADJ[7];
//         TR[7] = ADJ[5];
//         /*2d*/
//         // for(int i=0; i < N; ++i){
//         //     TR[i] = J[i];
//         // }
//         // TR[0] = J[3];
//         // TR[1] = -1.0*J[1];
//         // TR[2] = -1.0*J[2];
//         // TR[3] = J[0];

   
//     }

//     void Inv(){
//         adjoint(); //for3D uncomment and  for 2D comment
//         transpose();

//         Det();
//         for(int i=0; i < N; ++i){
//             INV[i] = TR[i] / DET;
//         }
//     }

//     ~Inverse(){
//         free(ADJ);
//         free(INV);
//         free(TR);
//         free(I);
//         free(J);
//     }
// };