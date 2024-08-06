#pragma once
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <fstream>
#include "node.h"
#include "state.h"
#include "grid.h"

void reading(std::string &path, grid& g){
    std::ifstream infile{path};
    for(int i = 0; i < g.NX*g.NY; ++i){
            infile >> g.domain[i].flag;
    }
}


void createFlags(grid& d){
    #pragma omp parallel for schedule(static)
    for(size_t j = 0; j < d.NY; ++j){
        for(size_t i = 0 ; i <d.NX; ++i){

            // outfile << circle(i,j,10,100,40)<<std::endl;

            if((j>=1) && (j<=d.NY-2) && (i<1)){//inlet
                d(j,i).flag =  10 ;

            }else if ((j>=1) && (j<=d.NY-2) && (i>d.NX-2))//outlet
            {
                d(j,i).flag =  10 ;
                
            }else if ((i>=0) && (i<=d.NX) && (j>d.NY-2))//top
            {
                d(j,i).flag =  10 ;

            }else if ((i>=0) && (i<=d.NX) && (j<1))//bottom
            {
                d(j,i).flag =  10 ;
            }
            else{//fluid node
                d(j,i).flag =  0 ;
            }
        }
    }
}


void WriteOutVTK(grid& g,const int& ITERATION,std::string& path,state d){
    
    
    if(ITERATION%d.vtk_out_freq==0){
    std::string part_out_PATH = path;                //change folder path
    std::string file_extension     = ".vtk";
    std::string vtk_out_name_base  = d.vtk_out_name_base;
    // if( (ITERATION % part_out_freq)==0){
    std::ofstream outfile{path+vtk_out_name_base+std::to_string(ITERATION)+ file_extension};
        // outfile << N <<"\n";
    outfile << "# vtk DataFile Version 2.0\n";
    outfile << "CLBM\n";
    outfile << "ASCII\n";
    outfile << "DATASET STRUCTURED_POINTS\n";
    outfile << "DIMENSIONS "<<g.NX << " "<<g.NY<<" "<<"1\n";
    outfile << "ORIGIN 0 0 0\n";
    outfile << "SPACING 1 1 1\n";
    outfile << "POINT_DATA "<< (g.NX)*(g.NY)<<"\n";

    
        
    outfile << "SCALARS flag int 1\n";
    outfile << "LOOKUP_TABLE default\n";
    for(int i = 0; i < g.NX*g.NY; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
        outfile << std::fixed << g.domain[i].flag << "\n";
        // }
    }

    outfile << "SCALARS density double 1\n";
    outfile << "LOOKUP_TABLE default\n";

    for(int i = 0; i < g.NX*g.NY; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
            outfile << std::fixed << g.domain[i].rho << "\n";
        // }
    }

    outfile << "SCALARS Temp double 1\n";
    outfile << "LOOKUP_TABLE default\n";

    for(int i = 0; i < g.NX*g.NY; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
            outfile << std::fixed << g.domain[i].T << "\n";
        // }
    }

    outfile << "SCALARS Mach_x double 1\n";
    outfile << "LOOKUP_TABLE default\n";

    for(int i = 0; i < g.NX*g.NY; ++i){
        if(g.domain[i].ux == 0 || g.domain[i].T ==0){
            outfile << std::fixed << 0.0 << "\n";
        }else{
            outfile << std::fixed << g.domain[i].ux / sqrt(d.sphr * g.domain[i].T) << "\n";
        }
    }

    // outfile << "SCALARS NRIter double 1\n";
    // outfile << "LOOKUP_TABLE default\n";
    // for(int i = 0; i < g.NX*g.NY; ++i){
    //     // if(g.domain[i].ux == 0 || g.domain[i].T ==0){
    //     //     outfile << std::fixed << 0.0 << "\n";
    //     // }else{
    //         outfile << std::fixed << g.domain[i]. << "\n";
    //     // }
    // }
    // outfile << "SCALARS Ma double 1\n";
    // outfile << "LOOKUP_TABLE default\n";

    // for(int i = 0; i < g.NX*g.NY; ++i){
    //     // if(g.domain[i].flag == 0){
    //         outfile << std::fixed << g.domain[i].T << "\n";
    //     // }
    // }
    outfile << "VECTORS velocity double\n";
    for(int i = 0; i < g.NX*g.NY; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
            outfile << std::fixed << g.domain[i].ux <<" " << g.domain[i].uy <<" "<< "0"<< "\n";
        // }
    }
        outfile.close(); 
    }  
}

void WriteOut(grid& g,const int& ITERATION,std::string& path,state d){
    
    std::string part_out_PATH = path;                //change folder path
    std::string file_extension     = ".out";
    std::string part_out_name_base     = d.part_out_name_base;
    // if( (ITERATION % part_out_freq)==0){
    std::ofstream outfile{path+part_out_name_base+std::to_string(ITERATION)+ file_extension};
        // outfile << N <<"\n";
    outfile << "Density" <<  std::setw(15) << "x vel" << std::setw(15) << "y vel" << std::setw(15) << "Temp" << std::setw(15) << "NR" << "\n";
    for(int i = 0; i < g.NX*g.NY; ++i){
        if(g.domain[i].flag == 0){
            
            outfile << std::fixed << g.domain[i].rho<< std::setw(15)<< g.domain[i].ux << std::setw(15)
                                  << g.domain[i].uy << std::setw(15)<< g.domain[i].T  << std::setw(15)<< g.domain[i].nr  << "\n";
        }
    }
        outfile.close();   
    // }
}
// 
