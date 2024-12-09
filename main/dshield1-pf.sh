#!/bin/bash
#PBS -q large
#PBS -l select=4:ncpus=32:mpiprocs=32
#PBS -l walltime=12:00:00
#PBS -N jobemses
##PBS -o stdout
#PBS -j oe

source /etc/profile.d/modules.sh
module load compiler mpt python-3.8.8-gcc-8.3.1-xbx7faq
export MPICC_CC='icc'
export MPI_TYPE_DEPTH=20

cd ${PBS_O_WORKDIR}

date

mpiexec -n 112 ./requester/mpiemses3D dshield1.inp : -n 1 ./coupler/coupler : -n 1 ./worker/dout

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

