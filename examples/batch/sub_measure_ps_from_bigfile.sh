#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=01:00:00
#PBS -N ps_hiAGN_BAH_2048_m_z0

###
# Example script (for CITA Sunnyvale cluster):
#
# Measuring power spectrum from density grid in bigfile format.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export ILL_DIR=/mnt/raid-cita/sforeman/bahamas/
export FILE_ROOT=Bahamas_hiAGN_delta_m_2048_z0.0

export IN_FILE=${ILL_DIR}sjf_grids/${FILE_ROOT}.bigfile
export OUT_JSON_FILE=${ILL_DIR}ps_bahamas/${FILE_ROOT}_ps_kf_0.5kf.json
export OUT_DAT_FILE=${ILL_DIR}ps_bahamas/${FILE_ROOT}_ps_kf_0.5kf.dat
export NGRID=2048
export NGRIDCIC=2048
export BOX_SIZE=400.
export DK_IN_KF=1.
export KMIN_IN_KF=0.5

source activate nbodykit-clean-env
cd $PBS_O_WORKDIR
mpirun -np 128 python measure_ps_from_bigfile.py $IN_FILE $OUT_JSON_FILE $OUT_DAT_FILE $NGRID $NGRIDCIC $BOX_SIZE $DK_IN_KF $KMIN_IN_KF