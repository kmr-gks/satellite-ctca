#!/bin/bash
#SBATCH -p gr20001a
#SBATCH --rsc p=128:t=1:c=1
#SBATCH -t 10:00
#SBATCH -o %x.%j.out

# 標準出力と標準エラー出力をリダイレクト
exec 2>&1

#環境変数を表示
echo "LD_LIBRARY_PATH="
echo $LD_LIBRARY_PATH
echo "PATH="
echo $PATH

# 自分自身のスクリプトファイルの内容を出力
echo "===== ジョブスクリプトの内容 ====="
cat "$0"
echo "===== 実行開始 ====="

#set -x
#export LD_LIBRARY_PATH="/LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/lib:/LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/lib:$LD_LIBRARY_PATH"

module load fftw
module load hdf5/1.12.2_intel-2022.3-impi

export EMSES_DEBUG=no

date

rm *_0000.h5
#dont use sbatch, use mysbatch and njob.sh instead
cp ../requester/bin/mpiemses3D .
srun ./mpiemses3D plasma.inp

date

# Postprocessing(visualization code, etc.)

echo ...done

python generate_xdmf3.py nd*.h5 rhobk00_0000.h5
python generate_xdmf3.py rho00_0000.h5
python generate_xdmf3.py phisp00_0000.h5
python generate_xdmf3.py ex00_0000.h5 ey00_0000.h5 ez00_0000.h5
python generate_xdmf3.py j1x00_0000.h5 j1y00_0000.h5 j1z00_0000.h5 j2x00_0000.h5 j2y00_0000.h5 j2z00_0000.h5 j3x00_0000.h5 j3y00_0000.h5 j3z00_0000.h5

date
