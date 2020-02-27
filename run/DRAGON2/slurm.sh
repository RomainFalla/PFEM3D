#!/bin/bash
# Submission script for Dragon2
#SBATCH --job-name=PFEM test
#SBATCH --time=05:00:00 # hh:mm:ss
#
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2625 # megabytes
#SBATCH --partition=batch,long
#
#SBATCH --mail-user=simon.fevrier@student.uliege.be
#SBATCH --mail-type=ALL
#
#SBATCH --comment=PFEM
#
#SBATCH --outfile=PFEMTestOut.txt

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

cd $HOME/PFEM-master/build/bin
srun ./main params/params.json geometry/square.msh results.msh