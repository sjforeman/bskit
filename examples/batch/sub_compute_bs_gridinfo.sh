#!/bin/bash -l
#PBS -l nodes=2:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=02:00:00
#PBS -N bs_gridinfo_1024_test

###
# Example script (for CITA Sunnyvale cluster):
#
# Computing bispectrum binning information (average k_i and number of triangles in each
# (k_1,k_2,k_3) triangle bin) for a given binning scheme. This file corresponds to
# a 1024^3 grid in a 75 h^-1 Mpc box, where each k_i is binned with dk = k_f (=2pi/Lbox)
# for the lowest 40 bins, and then dk = 6k_f for bins at higher k.
#
# This uses the faster (but more memory-intensive) routine.
###

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export SNAP_PREFIX=/mnt/raid-cita/sforeman/illustrisTNG/sjf_grids/TNG100-1_delta_m_1024_z0.0.bigfile
export NGRID=1024
export NGRIDCIC=2048
export LBOX=75.
export STARTI=0
export ENDI=25000
export OUTFILE=/mnt/raid-cita/sforeman/bs/ipy/bs_grid_info/Lbox75_1024_kf_6kf_40lowkbins.dat


export DK=0.0837758
export DK_HIGH=0.50265
export KMIN=0.0418879
export KMAX=23.5
export NUM_LOWK_BINS=40

source activate nbodykit-py3-env
source set_script_dir.sh
mpirun -np 128 python ${SCRIPT_DIR}measure/measure_bs_fast.py bigfile_grid $SNAP_PREFIX $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE --k_bin_info $KMIN $KMAX $DK -cic --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH  --triangle_type all --meas_type grid_info
