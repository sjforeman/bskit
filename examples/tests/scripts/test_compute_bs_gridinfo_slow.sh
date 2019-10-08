#!/bin/bash

###
# Example script (for command line):
#
# Computing bispectrum binning information (average k_i and number of triangles in each
# (k_1,k_2,k_3) triangle bin) for a given binning scheme.
#
# This uses the slower (but less memory-intensive) routine.
###

export OMP_NUM_THREADS=1

export SNAP=${GRID_DIR}test_grid_512_0.bigfile
export NGRID=512
export NGRIDCIC=512
export LBOX=1000.
export STARTI=0
export ENDI=200
export OUTFILE=${OUT_DIR}Lbox1000_512_kf_3kf_3lowkbins_slow.dat

export DK=0.00628
export DK_HIGH=0.01884
export KMIN=0.00314
export KMAX=0.1
export NUM_LOWK_BINS=3

mpirun -np 8 python ${SCRIPT_DIR}measure/measure_bs_slow.py bigfile_grid $SNAP $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE --k_bin_info $KMIN $KMAX $DK --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH --triangle_type all --meas_type grid_info
