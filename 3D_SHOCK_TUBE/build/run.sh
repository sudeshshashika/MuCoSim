#!/bin/bash
cd .. 
cd out
rm *.out

cd .. 
cd vtk
rm *.vtk
cd .. 

cd build
make clean
cmake ..
make
# ./GRID
export OMP_NUM_THREADS=2 ./CLBM





#************If running for the first time uncomment following otherwise keep it
# as it is for 3D casees since it takes some time to generate the mesh
# cd grid 
# rm gridGen gridGen.in
# g++ -o gridGen gridGen.cpp
# ./gridGen
# cd ..



# cd src
# rm a.out main.o timing.o timing.h.gch
 
# g++ -O3 -Wall -I/usr/local/include -c main.cpp timing.c
# g++ -L/usr/local/lib main.o timing.o -lgsl -lgslcblas -lm
# ./a.out
