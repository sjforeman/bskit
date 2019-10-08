#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=02:00:00
#PBS -N cross_hiAGN_BAH_2048_dm_b_z0

###
# Example script (for CITA Sunnyvale cluster):
#
# Measuring cross power spectrum of two density grids in bigfile format.
###

export Z=0

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export ILL_DIR=/mnt/raid-cita/sforeman/bahamas/
export FILE_ROOT=Bahamas_hiAGN_delta_dm_2048_z${Z}.0
export CROSS_FILE_ROOT=Bahamas_hiAGN_delta_b_2048_z${Z}.0

export IN_FILE=${ILL_DIR}sjf_grids/${FILE_ROOT}.bigfile
export CROSS_IN_FILE=${ILL_DIR}sjf_grids/${CROSS_FILE_ROOT}.bigfile
export OUT_JSON_FILE=${ILL_DIR}ps_bahamas/${FILE_ROOT}_cross_${CROSS_FILE_ROOT}_kf_0.5kf.json
export OUT_DAT_FILE=${ILL_DIR}ps_bahamas/${FILE_ROOT}_cross_${CROSS_FILE_ROOT}_kf_0.5kf.dat
export NGRID=2048
export NGRIDCIC=2048
export BOX_SIZE=400.
export DK_IN_KF=1.
export KMIN_IN_KF=0.5

source activate nbodykit-clean-env
cd $PBS_O_WORKDIR
mpirun -np 128 python measure_ps_from_bigfile.py $IN_FILE $OUT_JSON_FILE $OUT_DAT_FILE $NGRID $NGRIDCIC $BOX_SIZE $DK_IN_KF $KMIN_IN_KF --cross $CROSS_IN_FILE