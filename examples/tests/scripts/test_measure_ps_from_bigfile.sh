#!/bin/bash
###
# Example script (for command line):
#
# Measuring power spectrum from density grid in bigfile format.
###

export OMP_NUM_THREADS=1

export FILE_ROOT=test_grid_512_0

export IN_FILE=${GRID_DIR}${FILE_ROOT}.bigfile
export NGRID=512
export NGRIDCIC=512
export BOX_SIZE=1000.
export DK_IN_KF=1.
export KMIN_IN_KF=0.5
export OUT_JSON_FILE=${OUT_DIR}${FILE_ROOT}_ps_kf_0.5kf.json
export OUT_DAT_FILE=${OUT_DIR}${FILE_ROOT}_ps_kf_0.5kf.dat

mpirun -np 8 python ${SCRIPT_DIR}measure/measure_ps_from_bigfile.py $IN_FILE $OUT_JSON_FILE $OUT_DAT_FILE $NGRID $NGRIDCIC $BOX_SIZE $DK_IN_KF $KMIN_IN_KF