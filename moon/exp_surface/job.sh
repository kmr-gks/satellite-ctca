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

#set ship position and neighbour threshold [m]
export SHIP_COORD="(0,8,30)-(16,8,30)"

export GRID_LENGTH=0.5
export NEIGHBOUR_THR=1

#step range (step unit)
export STEP_FROM=1
export STEP_TO=10000

#set histogram parameters (1 for yes, 0 for no)
export CORRECT_BY_BIN_WIDTH=0

#for python script after simulation
export OUTPUT_DIR_NAME="${SLURM_JOB_ID}.${SHIP_COORD}.out"
export OUTPUT_FILE_NAME="${SHIP_COORD}.csv"

echo "output file: $OUTPUT_DIR_NAME"
mkdir $OUTPUT_DIR_NAME
cd $OUTPUT_DIR_NAME

mkdir -p ../output/${SLURM_JOB_ID}
cp $0 ../plasma.inp ../output/${SLURM_JOB_ID}/.
date

rm *_0000.h5
srun -l --multi-prog ../multi.conf

# Postprocessing(visualization code, etc.)

echo ...done

#move other output files
mv *.h5 chgacm1 chgacm2 chgmov currnt energy energy1 energy2 ewave icur influx isflux nesc noflux ocur oltime pbody pbodyd pbodyr plasma.out seyield SNAPSHOT1 volt ../output/${SLURM_JOB_ID}/.

echo "Running python script with $OUTPUT_DIR_NAME"
python ../histogram.py 

date
mv "../job.sh.${SLURM_JOB_ID}.out" .
