#!/bin/bash
#SBATCH -p gr10451a
#SBATCH --rsc p=1:t=1:c=1:m=40G
#SBATCH -t 168:00:00
#SBATCH -o %x.%j.out

#gr10451a or gr20001a

# 標準出力と標準エラー出力をリダイレクト
exec 2>&1

date
python plot.py
date
