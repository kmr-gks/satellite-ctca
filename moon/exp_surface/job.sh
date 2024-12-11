#!/bin/bash
#SBATCH -p gr20001a
#SBATCH --rsc p=130:t=1:c=1
#SBATCH -t 168:00:00
#SBATCH -o %x.%j.out

# 標準出力と標準エラー出力をリダイレクト
exec 2>&1

#set -x
#export LD_LIBRARY_PATH="/LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/lib:/LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/lib:$LD_LIBRARY_PATH"

module load fftw
module load hdf5/1.12.2_intel-2022.3-impi

export EMSES_DEBUG=no

date

rm *_0000.h5
#dont use sbatch, use mysbatch and njob.sh instead
srun -n130 -l --multi-prog multi.conf
date

# Postprocessing(visualization code, etc.)

echo ...done

date
