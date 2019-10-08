#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=24:00:00
#PBS -N conv_ill-1_z3_b_to_delta_m_bigfile_2048

###
# Example script (for CITA Sunnyvale cluster):
#
# Converting an Illustris particle snapshot in HDF5 format into a density grid
# in bigfile format.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export ILL_DIR=/mnt/raid-cita/sforeman/original_illustris/
export SNAP_ROOT=Illustris-1/snap060/snap_060.

export IN_SNAP=${ILL_DIR}${SNAP_ROOT}
export OUT_FILE=${ILL_DIR}sjf_grids/Illustris-1_delta_b_2048_z3.0.bigfile
export SNAP_TYPE=baryons
export NMESH=2048
export BOX_SIZE=75.
export MASS_DM=0.000440896524361096  # Illustris-1

source activate nbodykit-env
cd $PBS_O_WORKDIR
mpirun -np 128 python ../convert_illustris_hdf5_particles_to_bigfile_grid.py  $IN_SNAP $OUT_FILE $SNAP_TYPE $NMESH $BOX_SIZE -mDM $MASS_DM