#!/bin/bash

###
# Example script (for command line):
#
# Generating 2 Gaussian density grids for testing.
#
###

export OMP_NUM_THREADS=1

mpirun -n 8 python ${BSKIT_GRID_DIR}generate_gaussian_test_grids.py
