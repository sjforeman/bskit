#!/bin/bash -l
#PBS -l nodes=2:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=03:00:00
#PBS -N ps_bf_subboxes_orig_b

###
# Example script (for CITA Sunnyvale cluster):
#
# Measuring power spectrum from subboxes of a density grid in bigfile format.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export ILL_DIR=/mnt/raid-cita/sforeman/original_illustris/
export FILE_ROOT=Illustris-1_delta_b_1024_z3.0
export SUBBOX_STR=_8subboxes

export IN_FILE=${ILL_DIR}sjf_grids/${FILE_ROOT}.bigfile
export OUT_DAT_FILE=${ILL_DIR}ps_orig/${FILE_ROOT}${SUBBOX_STR}_ps_kf_0.5kf.dat
export NGRID=1024
export NGRIDCIC=2048
export BOX_SIZE=75.
export DK_IN_KF=1.
export KMIN_IN_KF=0.5
export NSUB_PER_SIDE=2

source activate nbodykit-clean-env
cd $PBS_O_WORKDIR
source set_script_dir.sh
mpirun -np 32 python ${SCRIPT_DIR}measure_ps_from_bigfile_subboxes.py $IN_FILE $OUT_DAT_FILE $NGRID $NGRIDCIC $BOX_SIZE $DK_IN_KF $KMIN_IN_KF $NSUB_PER_SIDE