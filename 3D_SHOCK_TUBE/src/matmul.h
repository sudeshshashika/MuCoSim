#pragma once

#include <iostream>

class Matmul{

public:
    int n;
    double   f[4];  // 4 functions
    double   M[16]; // 4*4 is matrix of 16 elements
    double res[4];  // 4 residuals

    Matmul(){}

    Matmul(int n)
        :n(n){
            for(int i = 0; i < (n * n);++i){
                M[i] = 0.0;
            }
            for(int i = 0; i < n;++i){
                  f[i] = 0.0;
                res[i] = 0.0;
            }
        }

    Matmul(const Matmul& src)
        :n(src.n){
            for(int i = 0; i < (n * n);++i){
                M[i] = src.M[i];
            }
            for(int i = 0; i < n;++i){
                  f[i] = src.f[i];
                res[i] = src.res[i];
            } 
        }

    Matmul& operator=(const Matmul& src){
        if(this == &src){
            return *this;
        }
        n = src.n;
        for(int i = 0; i < (n * n);++i){
            M[i] = src.M[i];
        }
        for(int i = 0; i < n;++i){
            f[i]   = src.f[i];
            res[i] = src.res[i];
        } 
        return *this;        
            
    }

    void Mul(double* matrix , double* array ){
        /*3d*/
        res[0] = matrix[0] *array[0] + matrix[1] *array[1] + matrix[2] *array[2]  + matrix[3] *array[3];
        res[1] = matrix[4] *array[0] + matrix[5] *array[1] + matrix[6] *array[2]  + matrix[7] *array[3];
        res[2] = matrix[8] *array[0] + matrix[9] *array[1] + matrix[10]*array[2]  + matrix[11]*array[3];
        res[3] = matrix[12]*array[0] + matrix[13]*array[1] + matrix[14]*array[2]  + matrix[15]*array[3];
        /*2*/
        // res[0] = matrix[0]*array[0] + matrix[1]*array[1];
        // res[1] = matrix[2]*array[0] + matrix[3]*array[1];
    }

    ~Matmul(){
        // free(f);
        // free(M);
        // free(res);
    }
        
};

// class Matmul{

// public:
//     int n;
//     double* f;
//     double* M;
//     double* res;

//     Matmul(){}

//     Matmul(int n)
//         :n(n){
//             f   = (double*)malloc(sizeof(double) * n);
//             res = (double*)malloc(sizeof(double) * n);
//             M   = (double*)malloc(sizeof(double) * (n * n));
//             for(int i = 0; i < (n * n);++i){
//                 M[i] = 0.0;
//             }
//             for(int i = 0; i < n;++i){
//                 f[i] = 0.0;
//                 res[i] = 0.0;
//             }
//         }

//     Matmul(const Matmul& src)
//         :n(src.n){
//             f = M = res = nullptr;
//             f   = (double*)malloc(sizeof(double) * n);
//             res = (double*)malloc(sizeof(double) * n);
//             M   = (double*)malloc(sizeof(double) * (n * n));
//             for(int i = 0; i < (n * n);++i){
//                 M[i] = src.M[i];
//             }
//             for(int i = 0; i < n;++i){
//                 f[i] = src.f[i];
//                 res[i] = 0.0;
//             } 
//         }

//     Matmul& operator=(const Matmul& src){
//         if(this == &src){
//             return *this;
//         }
//         n = src.n;
//         f = M = res = nullptr;
//         f   = (double*)malloc(sizeof(double) * n);
//         res = (double*)malloc(sizeof(double) * n);
//         M   = (double*)malloc(sizeof(double) * (n * n));
//         for(int i = 0; i < (n * n);++i){
//             M[i] = src.M[i];
//         }
//         for(int i = 0; i < n;++i){
//             f[i]   = src.f[i];
//             res[i] = src.res[i];
//         }         
            
//     }

//     void Mul(double* matrix , double* array ){
//         /*3d*/
//         res[0] = matrix[0]*array[0] + matrix[1]*array[1] + matrix[2]*array[2];
//         res[1] = matrix[3]*array[0] + matrix[4]*array[1] + matrix[5]*array[2];
//         res[2] = matrix[6]*array[0] + matrix[7]*array[1] + matrix[8]*array[2];
//         /*2*/
//         // res[0] = matrix[0]*array[0] + matrix[1]*array[1];
//         // res[1] = matrix[2]*array[0] + matrix[3]*array[1];
//     }

//     ~Matmul(){
//         free(f);
//         free(M);
//         free(res);
//     }
        
// };