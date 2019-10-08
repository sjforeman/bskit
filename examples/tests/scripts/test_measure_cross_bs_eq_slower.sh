#!/bin/bash
###
# Example script (for command line):
#
# Measuring bispectrum from specified density grid in bigfile format, using faster
# (but more memory-intensive) routine.
###

export OMP_NUM_THREADS=1

export FILE_ROOT1=test_grid_512_0
export FILE_ROOT2=test_grid_512_1
export SNAP1=${GRID_DIR}${FILE_ROOT1}.bigfile
export SNAP2=${GRID_DIR}${FILE_ROOT2}.bigfile

export OUTFILE112=${OUT_DIR}${FILE_ROOT1}_cross_${FILE_ROOT2}_bs_eq_kf_3kf_3lowbins_slow.dat
export OUTFILE221=${OUT_DIR}${FILE_ROOT2}_cross_${FILE_ROOT1}_bs_eq_kf_3kf_3lowbins_slow.dat

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

mpirun -np 8 python ${SCRIPT_DIR}measure/measure_bs_slow.py bigfile_grid $SNAP1 $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE112 --k_bin_info $KMIN $KMAX $DK --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH  --triangle_type equilateral --snap_prefix2 $SNAP2

mpirun -np 8 python ${SCRIPT_DIR}measure/measure_bs_slow.py bigfile_grid $SNAP2 $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE221 --k_bin_info $KMIN $KMAX $DK --high_k_bin_info $NUM_LOWK_BINS $DK_HIGH  --triangle_type equilateral --snap_prefix2 $SNAP1
