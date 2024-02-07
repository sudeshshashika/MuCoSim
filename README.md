
<!-----------------------------------------------------------------------------


This document should be written based on the Github flavored markdown specs:
https://github.github.com/gfm/
It can be converted to html or pdf with pandoc:
pandoc -s -o logbook.html  -f gfm -t html logbook.md
pandoc test.txt -o test.pdf
or with the kramdown converter:
kramdown --template document  -i GFM  -o html logbook.md

Optional: Document how much time was spent. A simple python command line tool
for time tracking is [Watson](http://tailordev.github.io/Watson/).

------------------------------------------------------------------------------>


# MuCoSim: Analysis of a Compressible Lattice Boltzmann Solver (CLBM)


[TOC]


## Project Description



### Customer Info


* Name: Thomas Gruber
* E-Mail: [thomas.gruber@fau.de](mailto:thomas.gruber@fau.de)


### Application Info

* Code: TeaLeaf
* URL: https://github.com/UK-MAC/TeaLeaf


Tealeaf is a miniapp solving linear heat conduction equation on a spatially decomposed regularly grid using a 5 point stencil with implicit solvers. The benchmark was part of the Mantevo Benchmark Suite that contains selected benchmarks of the National Labratories in the USA. The TeaLeaf mini-app itself was created and is maintained by researchers from the UK (https://ieeexplore.ieee.org/document/8049027)


> <media-tag title="Example of heat conduction. The arrow shows the direction of heat flow in the rod. By MikeRun, CC BY-SA 4.0 &lt;https://creativecommons.org/licenses/by-sa/4.0&gt;, via Wikimedia Commons" src="https://files.cryptpad.fr/blob/6a/6a7ab991aa406f811167e25701a769b8d999a68160ac30e4" data-crypto-key="cryptpad:PDJgMmn8suVztlmHa1sBk+ZyyBRXaGu3qz2IhMX1zEI="></media-tag>
> _Symbolic picture of linear heat conduction._ ([source, CC license](https://commons.wikimedia.org/wiki/File:Heat-conduction.svg))


But the heat conduction calculations are only a placeholder. The idea of the TeaLeaf mini-app is to explore the design space of new, scalable iterative sparse linear solvers which can target next generation architectures.



### Testsystem

* Host/Clustername:
* Cluster Info URL:
* CPU type:
* Memory capacity:
* Number of cores per node:
* Interconnect:

### Software Environment

* Compiler: 
* MPI: 
* Tools:
  * 

### How to build software

<!-----------------------------------------------------------------------------

Commands to build the software including flags, modules and required changes

------------------------------------------------------------------------------>

### Testcase description

<!-----------------------------------------------------------------------------

Describe the test case like domain size, number of iterations, ...

------------------------------------------------------------------------------>

### How to run software


<!-----------------------------------------------------------------------------

Commands to run the software including environment settings

------------------------------------------------------------------------------>

## Task1: Scaling runs


<!-----------------------------------------------------------------------------
For scaling runs, we find a runtime/performance number in the output of TeaLeaf (file tea.out). The last line always look like this:

 Wall clock    2.03252601623535

We use this value to plot the runtime of 10 iterations (see end-step in testcase definition) for 1 to 72 OpenMP threads inside a single MPI process (for i in {1..72}; do likwid-mpirun -np 1 -t $i ./tea_leaf; tail -n 1 tea.out | awk '{print $3}'; done).

We plot it with gnuplot:

set terminal png
set output 'openmp-scaling.png'
set title 'Weak scaling of TeaLeaf on Intel IcelakeSP 8360Y, Intel 19.0.5, default flags, ppcg method'
set xlabel '#Threads'
set xrange [1:72]
set ylabel 'Runtime [s]'
plot 'openmp-scaling.dat' w linespoints title 'Runtime'

------------------------------------------------------------------------------>


## Task2: Whole application measurements


<!-----------------------------------------------------------------------------

In order to get a first impression about the runtime behavior of the code, we use a basic set of LIKWID performance groups:

    DATA: Load-store-ratio
    FLOPS_DP/SP: Flops rate
    L2: L1 <-> L2 memory traffic
    L3: L2 <-> L3 memory traffic
    MEM: Memory traffic
    ENERGY: Energy consumption

Running: likwid-mpirun -np X -t Y -g GROUP ./tea_leaf

Analyze the results and come up with some main observations.

Note: In case of strange MPI related errors, try --mpiopts "-mpi=pmi2".

------------------------------------------------------------------------------>

## Task3: Runtime profile

<!-----------------------------------------------------------------------------

For the runtime profile, we commonly use gprof due to the already available support by the compilers. In order to get the sampling data, we add the compiler flag -pg and run with a single thread!

Identify the hotspots of the application.

------------------------------------------------------------------------------>
## Task4: Instrument kernels with MarkerAPI

<!-----------------------------------------------------------------------------

Put MarkerAPI calls around the hotspots.

------------------------------------------------------------------------------>

## Task5: Measurements of the selected hot spots


<!-----------------------------------------------------------------------------

Measure the hotspots with LIKWID


Running: likwid-mpirun -np X -t Y -m -g GROUP ./tea_leaf

------------------------------------------------------------------------------>


## Task6: Discussion of hot spot measurements


<!-----------------------------------------------------------------------------

Analyze the results of the measurements.

- Do you have ideas how to optimize the code?
- Any compile flag you would change based on the results?
- etc.
------------------------------------------------------------------------------>
