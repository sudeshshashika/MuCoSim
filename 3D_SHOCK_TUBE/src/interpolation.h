#pragma once

#include <iostream>
#include "math.h"
#include "node.h"
#include "grid.h"


void LPMiddleFluidnode(grid& g){
    // int endpx = g.NX-1;
    for(size_t j =0; j < g.NY; ++j){
        for(size_t i =2; i < g.NX-2; ++i){
            if( (g(j,i).ushift > 0.0) && (g(j,i).flag ==0) ){
                for(size_t q =0; q < 9; ++q){
                    g(j,i).fshift[q] = (1.0 - g(j,i).ushift)*(1.0 + g(j,i).ushift) * g(j,i).ushift * ( 0.5 * g(j,i-1).fold[q] + g(j,i).fold[q] + 0.5*g(j,i+1).fold[q]) ;
                    g(j,i).gshift[q] = (1.0 - g(j,i).ushift)*(1.0 + g(j,i).ushift) * g(j,i).ushift * ( 0.5 * g(j,i-1).gold[q] + g(j,i).gold[q] + 0.5*g(j,i+1).gold[q]) ;;
                }
            }
        }
    }
}

void LPExtrapolation(grid& g){
    for(size_t j =0; j < g.NY; ++j){
            int i = 1;
        // for(size_t i =2; i < g.NX-2; ++i){
            if( (g(j,i).ushift > 0.0) && (g(j,i).flag ==0) ){
                for(size_t q =0; q < 9; ++q){
                    g(j,i).fshift[q] = (1.0 + g(j,i).ushift)*(2.0 + g(j,i).ushift) * g(j,i).ushift * ( 0.5 * g(j,i).fold[q] + g(j,i+1).fold[q] + g(j,i+2).fold[q] );
                    g(j,i).gshift[q] = (1.0 + g(j,i).ushift)*(2.0 + g(j,i).ushift) * g(j,i).ushift * ( 0.5 * g(j,i).gold[q] + g(j,i+1).gold[q] + g(j,i+2).gold[q] );
                } 
            }
        // }
    }
}

void LPEndfluidnode(grid& g){
    for(size_t j =0; j < g.NY; ++j){
            int i = g.NX - 2;
            if( (g(j,i).ushift > 0.0) && (g(j,i).flag == 0) ){
                for(size_t q =0; q < 9; ++q){
                    g(j,i).fshift[q] = (1.0 - g(j,i).ushift)*(2.0 - g(j,i).ushift) * g(j,i).ushift * ( 0.5 * g(j,i).fold[q] + g(j,i-1).fold[q] + g(j,i-2).fold[q] );
                    g(j,i).gshift[q] = (1.0 - g(j,i).ushift)*(2.0 - g(j,i).ushift) * g(j,i).ushift * ( 0.5 * g(j,i).gold[q] + g(j,i-1).gold[q] + g(j,i-2).gold[q] );       
                }
            }
    }
}