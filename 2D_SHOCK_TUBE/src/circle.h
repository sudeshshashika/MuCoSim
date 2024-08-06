#pragma once
#include <iostream>
#include <omp.h>
#include <math.h>
#include "node.h"
#include "grid.h"


double radius(int x,int y, int i, int j){

    double b = pow(x-i,2) + pow(y-j,2); 
    double r = pow(b,0.5);
    return r;
}


void circleWall(grid& g,int x, int y,double Radius){
    #pragma omp parallel for
    for(size_t j =0;j<g.NY;++j){
        for(size_t i=0; i< g.NX;++i){
            
            if(radius(x,y,i,j)<Radius){
                g(j,i).flag = 5;
            }
            
        }
    }
}

void bccircleWall(grid& g,int x, int y,double Radius){

    
    for(size_t j =y-Radius-5;j<y+Radius+5;++j){
        for(size_t i=x-Radius-5; i< x+Radius+5;++i){
        
            if(g(j,i).flag == 5){
                   g(j,i).fold[1] = g(j,i+1).fold[3];
                   g(j,i).fold[2] = g(j+1,i).fold[4];
                   g(j,i).fold[3] = g(j,i-1).fold[1];
                   g(j,i).fold[4] = g(j-1,i).fold[2];
                   g(j,i).fold[5] = g(j+1,i+1).fold[7];
                   g(j,i).fold[6] = g(j+1,i-1).fold[8];
                   g(j,i).fold[7] = g(j-1,i-1).fold[5];
                   g(j,i).fold[8] = g(j-1,i+1).fold[6];

            }
            
        }
    }
}
