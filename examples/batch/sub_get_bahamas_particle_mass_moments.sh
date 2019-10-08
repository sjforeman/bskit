#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=01:00:00
#PBS -N BAH_lowAGN_mass_hydro_z1

###
# Example script (for CITA Sunnyvale cluster):
#
# Computing moments of particle mass distribution from BAHAMAS snapshot.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export SIM_DIR=/mnt/scratch-lustre/mccarthy/BAHAMAS_snapshots/
export SNAP_ROOT=BAHAMAS_low_AGN_nu0_L400N1024_WMAP9/snapshot_026/snap_026.
export IN_SNAP=${SIM_DIR}${SNAP_ROOT}

export OUT_DIR=/mnt/raid-cita/sforeman/bahamas/
export OUT_FILE=${OUT_DIR}sjf_grids/NONE
export SNAP_TYPE=hydro
export NMESH=2048
export BOX_SIZE=400.

source activate nbodykit-clean-env
cd $PBS_O_WORKDIR
mpirun -np 128 python ../convert_bahamas_hdf5_particles_to_bigfile_grid.py  $IN_SNAP $OUT_FILE $SNAP_TYPE $NMESH $BOX_SIZE --mass_moments_only