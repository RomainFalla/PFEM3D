#!/bin/bash
# Submission script for Dragon2
#SBATCH --job-name=PFEM test
#SBATCH --time=5:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2625 # megabytes
#SBATCH --partition=batch,long
#
#SBATCH --comment=PFEM
#
#SBATCH --outfile=PFEMTestOut.txt

module load CMake/3.12.1-GCCcore-7.3.0
module load CGAL/4.11.1-foss-2018b-Python-2.7.15

export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=32

cd $HOME/PFEM-master/build/bin
srun ./pfem ../../params/3D/damBreakKoshizuka.json ../../geometry/3D/damBreakKoshizuka.msh results.msh