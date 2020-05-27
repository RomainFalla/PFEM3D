#!/bin/bash
# Submission script for Dragon2
#SBATCH --account=fevrier
#SBATCH --job-name=KoshizukaDamBreak
#SBATCH --time=05:00:00 # hh:mm:ss
#
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2625 # megabytes
#SBATCH --partition=defq
#
#SBATCH --comment=PFEM
#
#SBATCH --output=PFEMTestOut.txt

module load gcc/9.3.0 boost/1.78.0 gmp/6.0.2 mpfr/4.0.2 cmake/current CGAL/4.14.3

export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12

cd $HOME/PFEM/build/bin
srun ./pfem ../../params/3D/damBreakKoshizuka.json ../../geometry/3D/damBreakKoshizuka.msh