#!/bin/bash -l
#SBATCH --job-name=CLBM
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=/home/hpc/muco/muco114h/App/phase02/mucosim/2D_SHOCK_TUBE/build/%j_%x.out

unset SLURM_EXPORT_ENV 
#module load gsl

# This pins openMP threads to physical cores on the socket
export OMP_PLACES=cores

# This places threads as close as possible to one another
export OMP_PROC_BIND=close


# export OMP_NUM_THREADS=1

#srun --cpu-freq=2400000-2400000:performance ./CLBM
# srun --cpu-freq=2400000-2400000:performance likwid-pin -C S0:0 ./CLBM
for i in {1..72};
 do
  export OMP_NUM_THREADS=$i
     
#     TODO: put benchmark stuff in here
#     srun --cpu-freq=2400000-2400000:performance likwid-pin -C S0:0-$i ./CLBM
  srun --cpu-freq=2400000-2400000:performance ./CLBM
done
# 
# srun --cpu-freq=2400000-2400000:performance likwid-pin -C S0:0-35@S1:0 ./CLB
# for i in {1..35};
#  do
#   export OMP_NUM_THREADS=$i
     
#     TODO: put benchmark stuff in here
#     srun --cpu-freq=2400000-2400000:performance likwid-pin -C S0:0-35@S1:0-$i ./CLBM
#   srun --cpu-freq=2400000-2400000:performance ./CLBM
# done
