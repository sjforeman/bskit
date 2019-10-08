#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=04:00:00
#PBS -N fftb_ill100-1_z3

###
# Example script (for CITA Sunnyvale cluster):
#
# Measuring bispectrum from specified density grid in bigfile format, using faster
# (but more memory-intensive) routine.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export Z=3
export TYPE=b

export NGRID=1024
export NGRIDCIC=2048

export SNAP_PREFIX=/mnt/raid-cita/sforeman/illustrisTNG/sjf_grids/TNG100-1_delta_${TYPE}_${NGRID}_z${Z}.0.bigfile
export LBOX=75.
export STARTI=0
export ENDI=2000000
export OUTFILE=/mnt/raid-cita/sforeman/illustrisTNG/bs/TNG100-1_delta_${TYPE}_${NGRID}_z${Z}.0_bs_kf_6kf_40lowbins.dat

export DK=0.0837758
export DK_HIGH=0.50265
export KMIN=0.0418879
export KMAX=23.5
export NUM_LOWK_BINS=40

source activate nbodykit-clean-env
cd $PBS_O_WORKDIR
source set_script_dir.sh
mpirun -np 128 python ${SCRIPT_DIR}measure_bs_fast.py bigfile_grid $SNAP_PREFIX $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE --k_bin_info $KMIN $KMAX $DK -cic --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH  --triangle_type all

