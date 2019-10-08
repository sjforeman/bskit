"""Convert density grid in bigfile format to a lower grid resolution.
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
def print_status(comm,start_time,message):
    if comm.rank == 0:
        elapsed_time = time.time() - start_time
        print('%d\ts: %s' % (elapsed_time,message))

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('in_file',
                    help='input bigfile')
parser.add_argument('out_file',
                    help='output bigfile')
parser.add_argument('n_mesh',type=int,
                    help='resolution of downsampled mesh')

# Parse arguments
args = parser.parse_args()
in_file = args.in_file
out_file = args.out_file
n_mesh = args.n_mesh

comm = CurrentMPIComm.get()

# Set start time, and print first status message
start_time = time.time()
print_status(comm,start_time,'Starting script')

# Import mesh from file
orig_mesh = nbk.BigFileMesh(in_file,'Field')
print_status(comm,start_time,'Imported mesh from %s with Nmesh %d' % (in_file,orig_mesh.attrs['Nmesh'][0]))

# Define new FieldMesh from downsampled version of original mesh
new_mesh = nbk.FieldMesh(orig_mesh.compute(mode='real',Nmesh=n_mesh))
print_status(comm,start_time,'Created new mesh with Nmesh %d' % new_mesh.attrs['Nmesh'][0])

# Save mesh to disk, as bigfile (actual computation of downsampling happens in this step)
new_mesh.save(out_file,dataset='Field',mode='real')
print_status(comm,start_time,'Wrote mesh to %s' % out_file)
