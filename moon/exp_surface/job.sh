#!/bin/bash
#SBATCH -p gr10451a
#SBATCH --rsc p=130:t=1:c=1
#SBATCH -t 168:00:00
#SBATCH -o %x.%j.out

#gr10451a or gr20001a

# 標準出力と標準エラー出力をリダイレクト
exec 2>&1

cat $0
cat plasma.inp

#set -x
#export LD_LIBRARY_PATH="/LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/lib:/LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/lib:$LD_LIBRARY_PATH"

module load fftw
module load hdf5/1.12.2_intel-2022.3-impi

#set environment variables
export EMSES_DEBUG=no

export SHIP_X_FROM=1
export SHIP_X_TO=10
export SHIP_Y_FROM=16
export SHIP_Y_TO=16
export SHIP_Z_FROM=256
export SHIP_Z_TO=256
export NEIGHBOUR_THR=10
#for python script after simulation
export OUTPUT_DIR_NAME="(${SHIP_X_FROM},${SHIP_Y_FROM},${SHIP_Z_FROM})-(${SHIP_X_TO},${SHIP_Y_TO},${SHIP_Z_TO})t${NEIGHBOUR_THR}.${SLURM_JOB_ID}.out"
export OUTPUT_FILE_NAME="output.csv"
export JOB_OUT_FILE="job.sh.${SLURM_JOB_ID}.out"

echo "output file: $OUTPUT_DIR_NAME"
mkdir $OUTPUT_DIR_NAME
cd $OUTPUT_DIR_NAME

date

rm *_0000.h5
srun -l --multi-prog ../multi.conf
date

# Postprocessing(visualization code, etc.)

echo ...done

#move other output files
mkdir -p ${SLURM_JOB_ID}_output
mv *.h5 chgacm1 chgacm2 chgmov currnt energy energy1 energy2 ewave icur influx isflux nesc noflux ocur oltime pbody pbodyd pbodyr plasma.inp plasma.out seyield SNAPSHOT1 volt ${SLURM_JOB_ID}_output/.

echo "Running python script with $OUTPUT_DIR_NAME"
python ../histogram.py 

date
