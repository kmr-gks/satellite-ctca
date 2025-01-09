#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr20001a
#SBATCH -t 00:10:00
#SBATCH --rsc p=13:t=1:c=1:m=1070M
#SBATCH -o %x.%j.out

#============ Shell Script ============
set -x

export I_MPI_FABRICS=ofi

date
rm *.txt
srun -n13 -c1 -l --multi-prog multi.conf
echo ...done
date

