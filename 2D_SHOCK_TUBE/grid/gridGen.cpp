#include <iostream>
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

    // int fluidDomainx = 3000;
    // int fluidDomainy = 4;
    // int helperCells  = 2;
    // int NX = fluidDomainx+helperCells;
    // int NY = fluidDomainy+helperCells;
    
    std::ofstream outfile{"gridGen.in"};

    // outfile << NX;
    // outfile << NY;
    for(size_t j = 0; j < d.NY; ++j){
        for(size_t i = 0 ; i <d.NX; ++i){

            // outfile << circle(i,j,10,100,40)<<std::endl;

            if((j>=1) && (j<=d.NY-2) && (i<1)){//inlet
                outfile<< std::fixed<< 10 << "\n";

            }else if ((j>=1) && (j<=d.NY-2) && (i>d.NX-2))//outlet
            {
                outfile<< std::fixed<< 10 << "\n";
                
            }else if ((i>=0) && (i<=d.NX) && (j>d.NY-2))//top
            {
                outfile<< std::fixed<< 10 << "\n";

            }else if ((i>=0) && (i<=d.NX) && (j<1))//bottom
            {
                outfile<< std::fixed<< 10 << "\n";
            }
            else{//fluid node
                outfile<< std::fixed<< 0 << "\n";
            }
        }
    }
    outfile.close();
}
