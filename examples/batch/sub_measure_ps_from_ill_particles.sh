#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=12:00:00
#PBS -N ps_part_Ill1_1024_z0.0_128c

###
# Example script (for CITA Sunnyvale cluster):
#
# Measuring power spectrum from particle snapshot in HDF5 format.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export ILL_DIR=/mnt/raid-cita/sforeman/original_illustris/
export SIM_DIR=Illustris-1
export FILE_ROOT=snap_135.hdf5.

export SNAP_PREFIX=${ILL_DIR}${SIM_DIR}/${FILE_ROOT}
export NMESH=1024
export OUT_JSON_FILE=${ILL_DIR}ps_orig/${SIM_DIR}_part_delta_m_ps_kf_0.5kf_${NMESH}.json
export OUT_DAT_FILE=${ILL_DIR}ps_orig/${SIM_DIR}_part_delta_m_ps_kf_0.5kf_${NMESH}.dat
export SIM_TYPE=hydro
export BOX_SIZE=75.
export DK_IN_KF=1.
export KMIN_IN_KF=0.5

source activate nbodykit-env
cd $PBS_O_WORKDIR
mpirun -np 128 python ../measure_ps_from_illustris_particles.py $SNAP_PREFIX $OUT_JSON_FILE $OUT_DAT_FILE $SIM_TYPE $NMESH $BOX_SIZE $DK_IN_KF $KMIN_IN_KF