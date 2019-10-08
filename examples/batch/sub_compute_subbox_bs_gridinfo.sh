#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=24:00:00
#PBS -N bs_gridinfo_1024_75_subbox

###
# Example script (for CITA Sunnyvale cluster):
#
# Computing bispectrum binning information (average k_i and number of triangles in each
# (k_1,k_2,k_3) triangle bin) for a given binning scheme. IN THE ORIGINAL BOX,
# this file corresponds to
# a 1024^3 grid in a 75 h^-1 Mpc box, where each k_i is binned with dk = k_f (=2pi/Lbox)
# for the lowest 40 bins, and then dk = 6k_f for bins at higher k.
#
# This file computes the binning information corresponding to measurements with the same
# binning scheme, but in a subbox of the original box, with side length and grid resolution
# of 1/NSUB_PER_SIDE times the original box values.
#
# This uses the faster (but more memory-intensive) routine.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export SNAP_PREFIX=/mnt/raid-cita/sforeman/illustrisTNG/sjf_grids/TNG100-1_delta_m_1024_z0.0.bigfile
export NGRID=1024
export NGRIDCIC=2048
export LBOX=75.
export NSUB_PER_SIDE=2
export STARTI=0
export ENDI=200000
export START_SUBBOX_IND=0
export END_SUBBOX_IND=0
export OUTFILE_PREFIX=/mnt/raid-cita/sforeman/bs/ipy/bs_grid_info/Lbox75_1024_kf_6kf_40lowkbins

export DK=0.0837758
export DK_HIGH=0.50265
export KMIN=0.0418879
export KMAX=23.5
export NUM_LOWK_BINS=40

source activate nbodykit-clean-env
cd $PBS_O_WORKDIR
source set_script_dir.sh
mpirun -np 128 python ${SCRIPT_DIR}measure_subbox_bs_fast.py bigfile_grid $SNAP_PREFIX $NGRID $NGRIDCIC $LBOX $NSUB_PER_SIDE $STARTI $ENDI $START_SUBBOX_IND $END_SUBBOX_IND $OUTFILE_PREFIX --k_bin_info $KMIN $KMAX $DK --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH  --triangle_type all --meas_type grid_info