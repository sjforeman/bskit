#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=02:00:00
#PBS -N ds_paco

###
# Example script (for CITA Sunnyvale cluster):
#
# Downsampling a bigfile density grid to a lower resolution (specified on the command line).
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

source activate nbodykit-clean-env
cd $PBS_O_WORKDIR

export N_MESH=${2}
export IN_FILE=/mnt/raid-cita/sforeman/paco_sims/sjf_grids/Paco${1}_DM_delta_m_2048_z0.0.bigfile
export OUT_FILE=/mnt/raid-cita/sforeman/paco_sims/sjf_grids/Paco${1}_DM_delta_m_${N_MESH}_z0.0.bigfile

mpirun -np 128 python ../downsample_bigfile_grid.py $IN_FILE $OUT_FILE $N_MESH


