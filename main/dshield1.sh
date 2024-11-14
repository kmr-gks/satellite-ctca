#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr10451a
#SBATCH -t 10:00
#SBATCH --rsc p=114:t=1:c=1
#SBATCH -o %x.%j.out

#============ Shell Script ============
#gr10451a
#標準エラー出力を標準出力にリダイレクト
exec 1>&2

set -x
export LD_LIBRARY_PATH="/home/b/b38257/kobeb4/satellite-ctca/main/requester:/LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/lib:/LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/lib:$LD_LIBRARY_PATH"

## -nの値は --rscで書いたpの値と同じにする
## -cの値は --rscで書いたcの値と同じにする 
date
srun -n114 -l --multi-prog dshield1.conf
echo ...done
#/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python generate_xdmf3.py ./requester-data/nd1p00_0000.h5 ./requester-data/nd2p00_0000.h5 ./requester-data/phisp00_0000.h5
#/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python generate_xdmf3.py ./requester-data/ex00_0000.h5 ./requester-data/ey00_0000.h5 ./requester-data/ez00_0000.h5
#/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python generate_xdmf3.py ./requester-data/j1x00_0000.h5 ./requester-data/j1y00_0000.h5 ./requester-data/j1z00_0000.h5
#/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python generate_xdmf3.py ./requester-data/j2x00_0000.h5 ./requester-data/j2y00_0000.h5 ./requester-data/j2z00_0000.h5
#/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python generate_xdmf3.py ./requester-data/rho00_0000.h5
#/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python3 generate_xdmf3.py j3x00_0000.h5 j3y00_0000.h5 j3z00_0000.h5
#/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python3 generate_xdmf3.py j4x00_0000.h5 j4y00_0000.h5 j4z00_0000.h5
#/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python3 generate_xdmf3.py j5x00_0000.h5 j5y00_0000.h5 j5z00_0000.h5
date

