#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <likwid-marker.h>
#include "math.h"
#include "newtonRaphson.h"
#include "lbmkernerls.h"
#include "node.h"
#include "grid.h"
#include "inverse.h"
#include "auxliray.h"
#include "interpolation.h"
// #include "timing.h"
#include "state.h"

//#define CHUNK_SIZE 1 

int main(){
    LIKWID_MARKER_INIT;
    #pragma omp parallel
    {
      LIKWID_MARKER_REGISTER("main");

    }
    std::string prmIn{"../input_parameters_ST_2.txt"};       //specify correct path
    // std::string pathin{"../build/gridGen.in"};           //specify correct path
    std::string pathvtk{"../vtk/"};                     //specify correct path
    std::string path{"../out/"};                        //specify correct path
    state prm;
    prm.setState(prmIn);
    prm.setDomain();
    prm.setCV();

    grid G(prm.NX,prm.NY);
    grid BC_G(prm.NX,prm.NY);
    // reading(pathin,G);
    // reading(pathin,BC_G);
    createFlags(G);
    createFlags(BC_G);
    shiftAsign(G,prm.u_shift);
    shiftAsign(BC_G,prm.u_shift);

    initializeLagrangians(G);
    initializationMacro(G,prm,prm.fluid_flag);
    setWeights(G);
    calcFeq(G,prm.fluid_flag);
    newtonRaphson(G,prm);
    calcGeq(G,prm.fluid_flag);
    initilizeDDF(G);
    calcRelaxationFactors(G,prm.sphr,prm.mu,prm.CV,prm.Pr);
    calc_parameters_from_PDF(G,prm.CV,prm.fluid_flag);
    // WriteOut(G,0,path,prm);
    // WriteOutVTK(G,0,pathvtk,prm);
    int latticeTime{0};
    // double START = getTimeStamp();
    const auto START = std::chrono::steady_clock::now();
    #pragma omp parallel
    {
      LIKWID_MARKER_START("main");

    }
    while(latticeTime <= prm.total_time){
        // std::cout << "+--ITER | " << latticeTime << std::endl;
        calc_parameters_from_PDF(G,prm.CV,prm.fluid_flag);
        calcRelaxationFactors(G,prm.sphr,prm.mu,prm.CV,prm.Pr);
        //WriteOutVTK(G,latticeTime,pathvtk,prm);
        setWeights(G);
        calcFeq(G,prm.fluid_flag);
        newtonRaphson(G,prm);
        calcGeq(G,prm.fluid_flag);
        pTensor(G);
        pEqTensor(G);
        calcQuasiEqG(G);
        //KnudsenNumberStability(G);
        collision(G);
        resetDDFshitS(G);
        IDW(G,prm.u_shift);
        IDWE(G,prm.u_shift);
        interpolatedDDF(G,prm.u_shift); 
        boundaryConditiontWallLR(G,BC_G);  
        boundaryConditiontWallTP(G);
        stream(G);
        swap(G);
        // if(latticeTime == prm.total_time){
        //     WriteOut(G,latticeTime,path,prm);
        //     WriteOutVTK(G,latticeTime,pathvtk,prm);
        // }
        latticeTime+=1;
    }
    #pragma omp parallel
    {
      LIKWID_MARKER_STOP("main");

    }
    const auto END = std::chrono::steady_clock::now();
    const auto DURATION = std::chrono::duration<double>(END - START);
    // double END = getTimeStamp();
    // std::cout << std::setfill('-') << std::setw(45) << "" << std::setfill(' ') << std::endl;
    // std::cout <<"Simulation END"<<std::endl;
    // std::cout <<"Simulation Time: "<< END-START << " s"<<std::endl;
    // std::cout <<"Simulation Time: "<<DURATION.count() << " s"<<std::endl;
    std::cout <<prm.NX << "    "<<prm.NY<<"    "<<DURATION.count() << "    "<< (prm.total_time * prm.NX * prm.NY)/(DURATION.count() * 1000000.0) <<std::endl;
    //prm.printState();
    LIKWID_MARKER_CLOSE;
    return 0;
}





