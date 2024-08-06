#pragma once
#include <iostream>
#include "grid.h"

void zOuterWall(grid& g){
    //this for front and back
    for (int j = 0; j < g.NY; ++j){
        for(int i = 0; i < g.NX; ++i){
            g(0,j,i).flag      = 10;
            g(g.NZ-1,j,i).flag = 10;
        }
    } 
}

void yOuterWall(grid& g){
    //this is for top and bottom
    for (int k = 0; k < g.NZ; ++k){
        for(int i = 0; i < g.NX; ++i){
            g(k,0,i).flag      = 10;
            g(k,g.NY-1,i).flag = 10;
        }
    } 
}

void xOuterWall(grid& g){
    //this is for left and right
    for (int k = 0; k < g.NZ; ++k){
        for(int j = 0; j < g.NY; ++j){
            g(k,j,0).flag      = 10;//20;
            g(k,j,g.NX-1).flag = 10;//30;
        }
    } 
}
