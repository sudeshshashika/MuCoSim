#pragma once
#include <iostream>
#include <omp.h>
#include <math.h>
#include "node.h"
#include "grid.h"


double radius(int x,int y,int z, int i, int j, int k){

    double b = pow(x-i,2) + pow(y-j,2) + pow(z-k,2); 
    double r = pow(b,0.5);
    return r;
}


void circleWall(grid& g,int x, int y, int z,double Radius){
    // #pragma omp parallel for
    for(int k = z-Radius-5; k < z+Radius+5; ++k){
        for(int j =y-Radius-5;j < y+Radius+5; ++j){
            for(int i=x-Radius-5; i < x+Radius+5; ++i){
                
                if(radius(x,y,z,i,j,k) < Radius){ // do not use less than or equal then there will be some edges
                    g(k,j,i).flag = 5;
                }
                
            }
        }
    }
}

void bccircleWall(grid& g,int x, int y,int z, double Radius){

    for(int k = z-Radius-5; k < z+Radius+5; ++k){
        for(int j =y-Radius-5;j < y+Radius+5; ++j){
            for(int i=x-Radius-5; i < x+Radius+5; ++i){
                if(g(k,j,i).flag == 5){
                    // g(k,j,i).fnew[0]   = g(k,j,i)    .fold[0];
                    g(k,j,i).fold[1]   = g(k,j,i+1)  .fold[2];   
                    g(k,j,i).fold[2]   = g(k,j,i-1)  .fold[1];
                    g(k,j,i).fold[3]   = g(k,j+1,i)  .fold[4];
                    g(k,j,i).fold[4]   = g(k,j-1,i)  .fold[3];
                    g(k,j,i).fold[5]   = g(k+1,j,i)  .fold[6];
                    g(k,j,i).fold[6]   = g(k-1,j,i)  .fold[5];
                    g(k,j,i).fold[7]   = g(k,j+1,i+1).fold[8];
                    g(k,j,i).fold[8]   = g(k,j-1,i-1).fold[7];
                    // 
                    g(k,j,i).fold[9]    = g(k+1,j,i+1).fold[10]   ;
                    g(k,j,i).fold[10]   = g(k-1,j,i-1).fold[9]  ;
                    g(k,j,i).fold[11]   = g(k+1,j+1,i).fold[12]  ;
                    g(k,j,i).fold[12]   = g(k-1,j-1,i).fold[11]  ;
                    g(k,j,i).fold[13]   = g(k,j-1,i+1).fold[14]  ;
                    g(k,j,i).fold[14]   = g(k,j+1,i-1).fold[13]  ;
                    g(k,j,i).fold[15]   = g(k-1,j,i+1).fold[16]  ;
                    g(k,j,i).fold[16]   = g(k+1,j,i-1).fold[15]  ;
                    g(k,j,i).fold[17]   = g(k-1,j+1,i).fold[18]  ;
                    // 
                    g(k,j,i).fold[18]   = g(k+1,j-1,i)  .fold[17]; //new it should be old
                    g(k,j,i).fold[19]   = g(k+1,j+1,i+1).fold[20];
                    g(k,j,i).fold[20]   = g(k-1,j-1,i-1).fold[19];
                    g(k,j,i).fold[21]   = g(k-1,j+1,i+1).fold[22];
                    g(k,j,i).fold[22]   = g(k+1,j-1,i-1).fold[21];
                    g(k,j,i).fold[23]   = g(k+1,j-1,i+1).fold[24];
                    g(k,j,i).fold[24]   = g(k-1,j+1,i-1).fold[23];
                    g(k,j,i).fold[25]   = g(k+1,j+1,i-1).fold[26];
                    g(k,j,i).fold[26]   = g(k-1,j-1,i+1).fold[25];
                    /*energy population*/
                    // g(k,j,i).gnew[0]   = g(k,j,i)    .gold[0];
                    g(k,j,i).gold[1]   = g(k,j,i+1)  .gold[2];   
                    g(k,j,i).gold[2]   = g(k,j,i-1)  .gold[1];
                    g(k,j,i).gold[3]   = g(k,j+1,i)  .gold[4];
                    g(k,j,i).gold[4]   = g(k,j-1,i)  .gold[3];
                    g(k,j,i).gold[5]   = g(k+1,j,i)  .gold[6];
                    g(k,j,i).gold[6]   = g(k-1,j,i)  .gold[5];
                    g(k,j,i).gold[7]   = g(k,j+1,i+1).gold[8];
                    g(k,j,i).gold[8]   = g(k,j-1,i-1).gold[7];
                    // 
                    g(k,j,i).gold[9]    = g(k+1,j,i+1).gold[10]   ;
                    g(k,j,i).gold[10]   = g(k-1,j,i-1).gold[9]  ;
                    g(k,j,i).gold[11]   = g(k+1,j+1,i).gold[12]  ;
                    g(k,j,i).gold[12]   = g(k-1,j-1,i).gold[11]  ;
                    g(k,j,i).gold[13]   = g(k,j-1,i+1).gold[14]  ;
                    g(k,j,i).gold[14]   = g(k,j+1,i-1).gold[13]  ;
                    g(k,j,i).gold[15]   = g(k-1,j,i+1).gold[16]  ;
                    g(k,j,i).gold[16]   = g(k+1,j,i-1).gold[15]  ;
                    g(k,j,i).gold[17]   = g(k-1,j+1,i).gold[18]  ;
                    // 
                    g(k,j,i).gold[18]   = g(k+1,j-1,i)  .gold[17]; //new it should be old
                    g(k,j,i).gold[19]   = g(k+1,j+1,i+1).gold[20];
                    g(k,j,i).gold[20]   = g(k-1,j-1,i-1).gold[19];
                    g(k,j,i).gold[21]   = g(k-1,j+1,i+1).gold[22];
                    g(k,j,i).gold[22]   = g(k+1,j-1,i-1).gold[21];
                    g(k,j,i).gold[23]   = g(k+1,j-1,i+1).gold[24];
                    g(k,j,i).gold[24]   = g(k-1,j+1,i-1).gold[23];
                    g(k,j,i).gold[25]   = g(k+1,j+1,i-1).gold[26];
                    g(k,j,i).gold[26]   = g(k-1,j-1,i+1).gold[25];
                }
            }  
        }
    }
}
