#!/usr/bin/env bash
#SBATCH --job-name=test_TDP2
#SBATCH --output=job_test.out
#SBATCH --error=job_test.err
#SBATCH -p defq
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --nodes=4 --ntasks-per-node=2

module load intel/mkl/64/11.2/2015.3.187
module load compiler/gcc/5.1.0
module load slurm/14.03.0
module load hardware/hwloc/1.11.0
module load mpi/openmpi/gcc/1.8.6-tm

cd /home/pichon/

sleep 50
echo "TEST"
