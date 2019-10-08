#!/bin/bash -l
#PBS -l nodes=8:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=12:00:00
#PBS -N bs_subboxes_ill_dm_paper_z0

###
# Example script (for CITA Sunnyvale cluster):
#
# Measuring bispectrum from subboxes of specified density grid in bigfile format, using faster
# (but more memory-intensive) routine.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export Z=0
export PART_TYPE=dm

export SNAP_PREFIX=/mnt/raid-cita/sforeman/illustrisTNG/sjf_grids/TNG100-1_delta_${PART_TYPE}_1024_z${Z}.0.bigfile
export NGRID=1024
export NGRIDCIC=2048
export LBOX=75.
export NSUB_PER_SIDE=2
export STARTI=0
export ENDI=2000000
export START_SUBBOX_IND=0
export END_SUBBOX_IND=7
export OUTFILE_PREFIX=/mnt/raid-cita/sforeman/illustrisTNG/bs/TNG100-1_delta_${PART_TYPE}_1024_z${Z}.0_bs_kf_6kf_40lowbins

export DK=0.0837758
export DK_HIGH=0.50265
export KMIN=0.0418879
export KMAX=23.5
export NUM_LOWK_BINS=40

source activate nbodykit-clean-env
cd $PBS_O_WORKDIR
source set_script_dir.sh
mpirun -np 128 python ${SCRIPT_DIR}measure_subbox_bs_fast.py bigfile_grid $SNAP_PREFIX $NGRID $NGRIDCIC $LBOX $NSUB_PER_SIDE $STARTI $ENDI $START_SUBBOX_IND $END_SUBBOX_IND $OUTFILE_PREFIX --k_bin_info $KMIN $KMAX $DK -cic --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH  --triangle_type all --meas_type unnorm_b_value