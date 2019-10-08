#!/bin/bash

###
# Example script (for command line):
#
# Combining bispectrum binning information and unnormalized values into single file.
###

export BIN_INFO_FILE=${OUT_DIR}Lbox1000_512_kf_3kf_3lowkbins.dat
export B_UNNORM_FILE=${OUT_DIR}test_grid_512_0_unnormbs_kf_3kf_3lowbins.dat
export KMAX=0.085
export OUT_FILE=${OUT_DIR}test_grid_512_0_bs_comb_kf_3kf_3lowbins.dat

python ${SCRIPT_DIR}process/process_fast_bs_measurement.py $BIN_INFO_FILE $B_UNNORM_FILE $KMAX $OUT_FILE