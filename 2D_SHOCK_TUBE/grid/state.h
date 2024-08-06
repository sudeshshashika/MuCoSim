// This Class define the input parameters for the 
#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

class state{

public:
       
    int NX;
    int NY;
    int NZ;

    std::string  part_input_file;
    std::string  vtk_out_name_base; 
    std::string  part_out_name_base;                   
    int          vtk_out_freq;   
    int          part_out_freq;           
    int          domain_size_x; 
    int          domain_size_y; 
    int          domain_size_z;   
    int          helper_cells;   
    int          total_time;   
    int          fluid_flag;    
    double       Pr;                   
    double       sphr;                               
    double       u_x;                 
    double       u_y;                 
    double       u_z;                 
    double       mu;                 
    double       u_shift;   
    double       CV;                 

void setState(std::string infile1Name){
    std::string temp; 
    std::fstream infile1(infile1Name);
    infile1 >> temp >> part_input_file;
    infile1 >> temp >> vtk_out_name_base;
    infile1 >> temp >> part_out_name_base;
    infile1 >> temp >> vtk_out_freq;  
    infile1 >> temp >> part_out_freq; 
    infile1 >> temp >> domain_size_x;     
    infile1 >> temp >> domain_size_y;  
    infile1 >> temp >> domain_size_z;
    infile1 >> temp >> helper_cells;
    infile1 >> temp >> total_time;
    infile1 >> temp >> fluid_flag;
    infile1 >> temp >> Pr;
    infile1 >> temp >> sphr;
    infile1 >> temp >> u_x;
    infile1 >> temp >> u_y;
    infile1 >> temp >> u_z;
    infile1 >> temp >> mu;
    infile1 >> temp >> u_shift;
    infile1.close();
}
void setDomain(){
    NX = domain_size_x + helper_cells;
    NY = domain_size_y + helper_cells;
    NZ = domain_size_z + helper_cells;
}

void setCV(){
    CV = 1.0 / (sphr - 1.0);
}

void printState(){
    std::cout << std::setfill('-') << std::setw(45) << "" << std::setfill(' ') << std::endl;
    std::cout <<  std::left << "Input file"            <<" "<<std::internal<< std::setw(40)   << part_input_file    << std::endl;    
    std::cout <<  std::left << ".*vtk output name"     <<" "<<std::internal<< std::setw(17)   << vtk_out_name_base  << std::endl;
    std::cout <<  std::left << ".*out output name"     <<" "<<std::internal<< std::setw(17)   << part_out_name_base << std::endl;
    std::cout <<  std::left << "Total Lattice nodes"   <<" "<<std::internal<< std::setw(16)   << (NX*NY*1)          << std::endl;     
    std::cout <<  std::left << "Total Time"            <<" "<<std::internal<< std::setw(24)   << total_time         << std::endl;
    std::cout <<  std::left << "Prandtl number"        <<" "<<std::internal<< std::setw(20)   << Pr                 << std::endl;
    std::cout <<  std::left << "Specifif heat ratio"   <<" "<<std::internal<< std::setw(12)   << sphr               << std::endl;                
    std::cout <<  std::left << "kinematic viscosity"   <<" "<<std::internal<< std::setw(16)   << mu                 << std::endl;                  
    std::cout << std::setfill('-') << std::setw(45) << "" << std::setfill(' ') << std::endl;            
}

};

