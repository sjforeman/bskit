#!/bin/bash
###
# Example script (for command line):
#
# Measuring bispectrum from specified density grid in bigfile format, using faster
# (but more memory-intensive) routine.
###

export OMP_NUM_THREADS=1

export FILE_ROOT=test_grid_512_0
export SNAP=${GRID_DIR}${FILE_ROOT}.bigfile
export OUTFILE_PREFIX=${OUT_DIR}${FILE_ROOT}_unnormbs_kf_3kf_3lowbins

export NGRID=512
export NGRIDCIC=512

export LBOX=1000.
export NSUB_PER_SIDE=2
export STARTI=0
export ENDI=100
export START_SUBBOX_IND=0
export END_SUBBOX_IND=7

export DK=0.00628
export DK_HIGH=0.01884
export KMIN=0.00314
export KMAX=0.1
export NUM_LOWK_BINS=3

mpirun -np 8 python ${SCRIPT_DIR}measure/measure_subbox_bs_fast.py bigfile_grid $SNAP $NGRID $NGRIDCIC $LBOX $NSUB_PER_SIDE $STARTI $ENDI $START_SUBBOX_IND $END_SUBBOX_IND $OUTFILE_PREFIX --k_bin_info $KMIN $KMAX $DK --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH  --triangle_type all --meas_type unnorm_b_value
