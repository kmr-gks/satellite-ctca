#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr10451a
#SBATCH -t 10:00
#SBATCH --rsc p=898:t=1:c=2:m=2140M
#SBATCH -o %x.%j.out

#============ Shell Script ============

#標準エラー出力を標準出力にリダイレクト
exec 1>&2

export EMSES_DEBUG=no

#set -x
module load hdf5/1.12.2_intel-2022.3-impi
module load fftw/3.3.10_intel-2022.3-impi

## -nの値は --rscで書いたpの値と同じにする
## -cの値は --rscで書いたcの値と同じにする 
date
srun -n898 -c2 -l --multi-prog multi.conf
#srun -l ./requester/bin/mpiemses3D plasma.inp

echo ...done

date

