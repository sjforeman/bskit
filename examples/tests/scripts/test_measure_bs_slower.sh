#!/bin/bash
###
# Example script (for command line):
#
# Measuring bispectrum from specified density grid in bigfile format, using faster
# (but more memory-intensive) routine.
###

export OMP_NUM_THREADS=1

export FILE_ROOT=test_grid_512_0
export SNAP=${BSKIT_GRID_DIR}${FILE_ROOT}.bigfile
export OUTFILE=${BSKIT_OUT_DIR}${FILE_ROOT}_bs_kf_3kf_3lowbins_slow.dat

export NGRID=512
export NGRIDCIC=512

export LBOX=1000.
export STARTI=0
export ENDI=100

export DK=0.00628
export DK_HIGH=0.01884
export KMIN=0.00314
export KMAX=0.1
export NUM_LOWK_BINS=3

mpirun -np 8 python ${BSKIT_SCRIPT_DIR}measure/measure_bs_slow.py bigfile_grid $SNAP $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE --k_bin_info $KMIN $KMAX $DK --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH  --triangle_type all

