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
#SBATCH --mail-user=simon.fevrier@student.uliege.be
#SBATCH --mail-type=ALL
#
#SBATCH --comment=PFEM
#
#SBATCH --output=PFEMTestOut.txt

export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12

cd $HOME/PFEM/build/bin
srun ./pfem ../../params/3D/damBreakKoshizuka.json ../../geometry/3D/damBreakKoshizuka.msh results.msh