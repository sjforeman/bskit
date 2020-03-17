#!/bin/bash
###
# Example script (for command line):
#
# Measuring power spectrum from subboxes of a density grid in bigfile format.
###

export OMP_NUM_THREADS=1

export FILE_ROOT=test_grid_512_0
export SUBBOX_STR=_8subboxes

export IN_FILE=${BSKIT_GRID_DIR}${FILE_ROOT}.bigfile
export OUT_DAT_FILE=${BSKIT_OUT_DIR}${FILE_ROOT}${SUBBOX_STR}_ps_kf_0.5kf.dat

export NGRID=512
export NGRIDCIC=512
export BOX_SIZE=1000.
export DK_IN_KF=1.
export KMIN_IN_KF=0.5
export NSUB_PER_SIDE=2

mpirun -np 8 python ${BSKIT_SCRIPT_DIR}measure/measure_ps_from_bigfile_subboxes.py $IN_FILE $OUT_DAT_FILE $NGRID $NGRIDCIC $BOX_SIZE $DK_IN_KF $KMIN_IN_KF $NSUB_PER_SIDE
