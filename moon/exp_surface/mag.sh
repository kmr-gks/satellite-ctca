#!/bin/bash
#SBATCH -p gr20001a
#SBATCH --rsc p=66:t=1:c=1:m=8G
#SBATCH -t 168:00:00
#SBATCH -o %x.%j.out

#gr10451a or gr20001a

# 標準出力と標準エラー出力をリダイレクト
exec 2>&1

cat $0
cat mag.inp

#set -x
#export LD_LIBRARY_PATH="/LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/lib:/LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/lib:$LD_LIBRARY_PATH"

module load fftw
module load hdf5/1.12.2_intel-2022.3-impi

#set environment variables
export EMSES_DEBUG=no

#set ship position and neighbour threshold [m]
export SHIP_NUM=34
export SHIP_COORD1="(64,0,25)-(64,128,25)"
export SHIP_COORD2="(64,0,26)-(64,128,26)"
export SHIP_COORD3="(64,0,27)-(64,128,27)"
export SHIP_COORD4="(64,0,28)-(64,128,28)"
export SHIP_COORD5="(64,0,29)-(64,128,29)"
export SHIP_COORD6="(64,0,29.2)-(64,128,29.2)"
export SHIP_COORD7="(64,0,29.4)-(64,128,29.4)"
export SHIP_COORD8="(64,0,29.6)-(64,128,29.6)"
export SHIP_COORD9="(64,0,29.8)-(64,128,29.8)"
export SHIP_COORD10="(64,0,30)-(64,128,30)"
export SHIP_COORD11="(64,0,30.2)-(64,128,30.2)"
export SHIP_COORD12="(64,0,30.4)-(64,128,30.4)"
export SHIP_COORD13="(64,0,30.6)-(64,128,30.6)"
export SHIP_COORD14="(64,0,30.8)-(64,128,30.8)"
export SHIP_COORD15="(64,0,31)-(64,128,31)"
export SHIP_COORD16="(64,0,32)-(64,128,32)"
export SHIP_COORD17="(64,0,33)-(64,128,33)"
export SHIP_COORD18="(64,0,34)-(64,128,34)"
export SHIP_COORD19="(64,0,35)-(64,128,35)"
export SHIP_COORD20="(64,0,40)-(64,128,40)"
export SHIP_COORD21="(64,0,45)-(64,128,45)"
export SHIP_COORD22="(64,0,50)-(64,128,50)"
export SHIP_COORD23="(64,0,55)-(64,128,55)"
export SHIP_COORD24="(64,0,60)-(64,128,60)"
export SHIP_COORD25="(256,64,30)-(0,64,30)"
export SHIP_COORD26="(256,64,31)-(0,64,31)"
export SHIP_COORD27="(256,64,32)-(0,64,32)"
export SHIP_COORD28="(256,64,33)-(0,64,33)"
export SHIP_COORD29="(256,64,34)-(0,64,34)"
export SHIP_COORD30="(256,64,35)-(0,64,35)"
export SHIP_COORD31="(256,64,40)-(0,64,40)"
export SHIP_COORD32="(256,64,45)-(0,64,45)"
export SHIP_COORD33="(256,64,50)-(0,64,50)"
export SHIP_COORD34="(256,64,55)-(0,64,55)"

export GRID_LENGTH=0.5
export NEIGHBOUR_THR=1

#step range (step unit)
export STEP_FROM=99000
export STEP_TO=100000


#set histogram parameters (1 for yes, 0 for no)
export CORRECT_BY_BIN_WIDTH=0

#for python script after simulation
export OUTPUT_DIR_NAME="${SLURM_JOB_ID}mag.out"

echo "output file: $OUTPUT_DIR_NAME"
mkdir $OUTPUT_DIR_NAME
cd $OUTPUT_DIR_NAME

mkdir -p ../output/${SLURM_JOB_ID}
cp $0 ../output/${SLURM_JOB_ID}/.
cp ../mag.inp ../output/${SLURM_JOB_ID}/plasma.inp
date

rm *_0000.h5
srun -l --multi-prog ../mag.conf

# Postprocessing(visualization code, etc.)

echo ...done

#move other output files
mkdir -p ../output/${SLURM_JOB_ID}
mv *.h5 chgacm1 chgacm2 chgmov currnt energy energy1 energy2 ewave icur influx isflux nesc noflux ocur oltime pbody pbodyd pbodyr plasma.out seyield SNAPSHOT1 volt ../output/${SLURM_JOB_ID}/.

echo "Running python script with $OUTPUT_DIR_NAME"
for i in $(seq 1 $SHIP_NUM); do
	eval "SHIP_COORD=\$SHIP_COORD$i"
	python ../histogram.py "${SHIP_COORD}.csv"
	echo python "${SHIP_COORD}.csv"
done

date
mv "../mag.sh.${SLURM_JOB_ID}.out" .
