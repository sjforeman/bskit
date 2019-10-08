"""Convert numpy-format grid of \rho/\rho_mean to density grid in bigfile format.

Output format is a grid of \delta(\vec{x}), *NOT* corrected for the CIC mass
assignment (this needs to be corrected separately before any statistics are measured).
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py
import time

import nbodykit.lab as nbk
from nbodykit.source.catalog.file import Gadget1Catalog,HDFCatalog
from nbodykit import CurrentMPIComm

# Routine to print script status to command line, with elapsed time
def print_status(start_time,message):
    elapsed_time = time.clock() - start_time
    print('%d\ts: %s' % (elapsed_time,message))

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('in_file',
                    help='input hdf5 file')
parser.add_argument('out_file',
                    help='output bigfile')
parser.add_argument('box_size',type=float,
                    help='box size in Mpc/h')

# Parse arguments
args = parser.parse_args()
in_file = args.in_file
out_file = args.out_file
box_size = args.box_size

# Set start time, and print first status message
start_time = time.clock()
print_status(start_time,'Starting conversion')

# Read in grid values from npy file
in_grid = np.load(in_file)
print_status(start_time,'Read input grid from %s' % in_file)

# Subtract 1 grid values, since they are actually \rho/\rho_mean
in_grid -= 1.
print_status(start_time,'Subtracted 1 from input grid')

# Since the full array will be hosted in memory on a single rank, we can use
# ArrayMesh. After the ArrayMesh is initialized, the data will be spread across
# multiple ranks
mesh = nbk.ArrayMesh(in_grid,BoxSize=box_size)
print_status(start_time,'Converted grid to ArrayMesh')

# Save mesh to disk, as bigfile
mesh.save(out_file,dataset='Field',mode='real')
print_status(start_time,'Wrote mesh to %s' % out_file)
