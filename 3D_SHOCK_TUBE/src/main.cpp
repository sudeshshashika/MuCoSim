#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <chrono>
#include <likwid-marker.h>
#include "boundaryStup.h"
#include "circle.h"
#include "math.h"
#include "newtonRaphson.h"
#include "lbmkernerls.h"
#include "node.h"
#include "grid.h"
#include "auxliray.h"
//#include "timing.h"
#include "state.h"

int main(){
    LIKWID_MARKER_INIT;
    #pragma omp parallel
    {
      LIKWID_MARKER_REGISTER("MAIN");

    }
    std::string prmIn{"../input_parameters_ST_2.txt"};
    std::string path{"../out/"};                     //specify correct path
    // std::string pathin{"../build/gridGen.in"};         //specify correct path
    std::string pathvtk{"../vtk/"};                  //specify correct path
    state prm;
    prm.setState(prmIn);
    prm.setDomain();
    prm.setCV();
    prm.setUshift();

    grid G(prm.NX,prm.NY,prm.NZ);
    grid BC_G(prm.NX,prm.NY,prm.NZ);
    // reading(pathin,G);
    // reading(pathin,BC_G);
    xOuterWall(G);
    yOuterWall(G);
    zOuterWall(G);
    xOuterWall(BC_G);
    yOuterWall(BC_G);
    zOuterWall(BC_G);
    shiftAsign(G,prm.u_shift);
    shiftAsign(BC_G,prm.u_shift);

    initializeLagrangians(G);
    initializationMacro(G,prm,prm.fluid_flag);
    setWeights(G);
    calcFeq(G,prm.fluid_flag);
    newtonRaphson(G,prm.CV,prm.fluid_flag);
    calcGeq(G,prm.fluid_flag);
    initilizeDDF(G);
    calcRelaxationFactors(G,prm.sphr,prm.mu,prm.CV,prm.Pr);
    density(G,prm.fluid_flag);
    velocity(G,prm.fluid_flag);
    temperature(G,prm.CV,prm.fluid_flag);
    // WriteOut(G,0,path,prm);
    // WriteOutVTK(G,0,pathvtk,prm);

    int latticeTime = 0;
  //  double START = getTimeStamp();
    const auto START = std::chrono::steady_clock::now();
//     LIKWID_MARKER_START("MAIN");
    // #pragma omp parallel
    // {
    //   LIKWID_MARKER_START("MAIN");

    // }
    while(latticeTime <= prm.total_time){
        //std::cout << "ITER = " << latticeTime << std::endl;
        density(G,prm.fluid_flag);
        velocity(G,prm.fluid_flag);
        temperature(G,prm.CV,prm.fluid_flag);
        calcRelaxationFactors(G,prm.sphr,prm.mu,prm.CV,prm.Pr);
        // WriteOutVTK(G,latticeTime,pathvtk,prm);
        setWeights(G);
        calcFeq(G,prm.fluid_flag);
        // LIKWID_MARKER_START("newtonRaphson");
        newtonRaphson(G,prm.CV,prm.fluid_flag);
        // LIKWID_MARKER_START("newtonRaphson");
        calcGeq(G,prm.fluid_flag);
        pTensor(G);
        calcQuasiEqG(G);
        // KnudsenNumberStability(G);
        collision(G);
        resetDDFshitS(G);
        IDW(G,prm.u_shift);
        IDWE(G,prm.u_shift);
        interpolatedDDF(G,prm.u_shift); 
        boundaryConditiontWallLR(G,BC_G); 
        boundaryConditiontWallTP(G);
        boundaryConditiontWallFB(G);
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
      LIKWID_MARKER_STOP("MAIN");

    }
//     LIKWID_MARKER_STOP("MAIN");

    //double END = getTimeStamp();
    const auto END = std::chrono::steady_clock::now();
    const auto DURATION = std::chrono::duration<double>(END - START);
    std::cout <<prm.NX << "    "<<prm.NY<<"    "<<prm.NZ<<"    "<<DURATION.count() << "    "<< (prm.total_time * prm.NX * prm.NY * prm.NZ)/(DURATION.count() * 1000000.0) <<std::endl;
//     std::cout << std::setfill('-') << std::setw(45) << "" << std::setfill(' ') << std::endl;
//     std::cout <<"Simulation END"<<std::endl;    
//     std::cout <<"Simulation Time: "<<DURATION.count() << " s"<<std::endl;
//     prm.printState();
    //std::cout <<"Calculation Time: "<< END-START <<std::endl;
    LIKWID_MARKER_CLOSE;
    return 0;
}

