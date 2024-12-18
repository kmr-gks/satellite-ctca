#!/bin/bash
#SBATCH -p gr10451a
#SBATCH --rsc p=130:t=1:c=1
#SBATCH -t 168:00:00
#SBATCH -o %x.%j.out

# 標準出力と標準エラー出力をリダイレクト
exec 2>&1

cat $0

#set -x
#export LD_LIBRARY_PATH="/LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/lib:/LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/lib:$LD_LIBRARY_PATH"

module load fftw
module load hdf5/1.12.2_intel-2022.3-impi

#set environment variables
export EMSES_DEBUG=no

export SHIPY=16
export SHIPZ=256
export NEIGHBOUR_THR=10
export OUTPUT_FILE_NAME="output,sy=${SHIPY},sz=${SHIPZ},nt=${NEIGHBOUR_THR}"
export EXTENTION=".csv"

# check if the output file exists
NEW_FILE_NAME="${OUTPUT_FILE_NAME}${EXTENTION}"
COUNTER=0
while [ -f "$NEW_FILE_NAME" ]; do
	COUNTER=$((COUNTER+1))
	NEW_FILE_NAME="${OUTPUT_FILE_NAME}_${COUNTER}${EXTENTION}"
done

OUTPUT_FILE_NAME="${NEW_FILE_NAME}"
echo "output file: $OUTPUT_FILE_NAME"

date

rm *_0000.h5
#dont use sbatch, use mysbatch and njob.sh instead
srun -l --multi-prog multi.conf
date

# Postprocessing(visualization code, etc.)

echo ...done

if [ -f "$OUTPUT_FILE_NAME" ]; then
	# ファイルサイズを取得 (バイト単位)
	FILE_SIZE=$(stat -c%s "$OUTPUT_FILE_NAME")
	if [ "$FILE_SIZE" -ge 1024 ]; then
		echo "Running python script..."
		python histogram.py 
	else
		echo "Size of file '$OUTPUT_FILE_NAME' is too small."
	fi
else
	echo "File '$OUTPUT_FILE_NAME' does not exist."
fi

date
