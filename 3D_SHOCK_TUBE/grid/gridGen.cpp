#include <iostream>
#include <omp.h>
#include <iomanip>
#include <fstream>
#include <math.h> 
#include "state.h"
/*This code setup the grid for the the compressible Lattice Boltzmann Solver. Here the 
boundary conditions is also defined by using a flag field. The respective flag field for thr relevent 
boundary condition has been stated as follws.

no-slip    =  10;
inlet      =  20;
outlet     =  30;
periodic   =  40;
fluid cell =  0;

Helper cells is to handle the boundary layer. Thus, the initial domain is extened by extra 2 cells in each direction.  

*/
int main(){
    /*rectangle*/
    std::string prmIn{"../input_parameters_ST_2.txt"};
    state d;
    d.setState(prmIn);
    d.setDomain();
    // int fluidDomainx = 400;
    // int fluidDomainy = 80;
    // int fluidDomainz = 80;
    // int helperCells  = 2;
    // int NX = fluidDomainx+helperCells;
    // int NY = fluidDomainy+helperCells;
    // int NZ = fluidDomainy+helperCells;
    
    std::ofstream outfile{"gridGen.in"};

    // outfile << NX;
    // outfile << NY;
    #pragma omp parallel for  
    for(int k = 0; k < d.NZ; ++k){
        for(int j = 0; j < d.NY; ++j){
            for(int i = 0 ; i < d.NX; ++i){

                // outfile << circle(i,j,10,100,40)<<std::endl;
                outfile << i << std::setw(5) << j << std::setw(5) << k <<std::endl;

            }
        }
    }
    outfile.close();
}


/*
                if((j>=1) && (j<=NY-2) && (i<1)){//inlet
                    outfile<< std::fixed<< 20 << "\n";

                }else if ((j>=1) && (j<=NY-2) && (i>NX-2))//outlet
                {
                    outfile<< std::fixed<< 30 << "\n";
                    
                }else if ((i>=0) && (i<=NX) && (j>NY-2))//top
                {
                    outfile<< std::fixed<< 10 << "\n";

                }else if ((i>=0) && (i<=NX) && (j<1))//bottom
                {
                    outfile<< std::fixed<< 10 << "\n";
                }
                else{//fluid node
                    outfile<< std::fixed<< 0 << "\n";
                }

*/