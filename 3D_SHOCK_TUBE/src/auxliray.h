#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include "node.h"
#include "grid.h"
#include "state.h"

// 
void reading(std::string &path, grid& g){
    std::ifstream infile{path};
    for(int i = 0; i < g.NX*g.NY*g.NZ; ++i){
            infile >> g.domain[i].x >> g.domain[i].y >> g.domain[i].z;
    }
}
// 
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
    outfile << "DIMENSIONS "<<g.NX << " "<<g.NY<<" "<<g.NZ<<"\n";
    outfile << "ORIGIN 0 0 0\n";
    outfile << "SPACING 1 1 1\n";
    outfile << "POINT_DATA "<< (g.NX)*(g.NY)*(g.NZ)<<"\n";

    
        
    outfile << "SCALARS flag int 1\n";
    outfile << "LOOKUP_TABLE default\n";
    for(int i = 0; i < g.NX*g.NY*g.NZ; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
        outfile << std::fixed << g.domain[i].flag << "\n";
        // }
    }

    outfile << "SCALARS density double 1\n";
    outfile << "LOOKUP_TABLE default\n";

    for(int i = 0; i < g.NX*g.NY*g.NZ; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
            outfile << std::fixed << g.domain[i].rho << "\n";
        // }
    }

    outfile << "SCALARS Temp double 1\n";
    outfile << "LOOKUP_TABLE default\n";

    for(int i = 0; i < g.NX*g.NY*g.NZ; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
            outfile << std::fixed << g.domain[i].T << "\n";
        // }
    }
    outfile << "SCALARS alpha double 1\n";
    outfile << "LOOKUP_TABLE default\n";

    for(int i = 0; i < g.NX*g.NY*g.NZ; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
            outfile << std::fixed << g.domain[i].alphaf << "\n";
        // }
    }

    outfile << "SCALARS Mach_x double 1\n";
    outfile << "LOOKUP_TABLE default\n";

    for(int i = 0; i < g.NX*g.NY*g.NZ; ++i){
        if(g.domain[i].ux == 0 || g.domain[i].T <= 0){
            outfile << std::fixed << 0.0 << "\n";
        }else{
            outfile << std::fixed << g.domain[i].ux / sqrt(d.sphr * g.domain[i].T) << "\n";
        }
    }

    outfile << "VECTORS velocity double\n";
    for(int i = 0; i < g.NX*g.NY*g.NZ; ++i){
        // if(g.domain[i].flag == 0 || g.domain[i].flag == 5){
            outfile << std::fixed << g.domain[i].ux <<" " << g.domain[i].uy <<" "<< g.domain[i].uz<< "\n";
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
    outfile << "Density" <<  std::setw(15) << "x vel" << std::setw(15) << "y vel" << std::setw(15) << "Temp" << "\n";
    for(int i = 0; i < g.NX*g.NY*g.NZ; ++i){
        // if(g.domain[i].flag == 0){
            
            outfile << std::fixed << g.domain[i].rho<< std::setw(15)<< g.domain[i].ux <<std::setw(15)
                                  << g.domain[i].uy << std::setw(15)<< g.domain[i].uz <<std::setw(15)<< g.domain[i].T << "\n";
        // }
    }
        outfile.close();   
    // }
}