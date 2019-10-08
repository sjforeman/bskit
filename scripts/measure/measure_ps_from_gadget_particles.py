"""Measure matter power spectrum from particle snapshot in Gadget format.
"""

import numpy as np
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
parser.add_argument('snap_prefix',
                    help='input Gadget snapshot (up to final dot in filenames)')
parser.add_argument('out_json_file',
                    help='output power spectrum file, in nbodykit json format')
parser.add_argument('out_dat_file',
                    help='output power spectrum file, in plain text format')
parser.add_argument('nmesh',
                    help='Nmesh for measurement')
parser.add_argument('boxsize',type=float,
                    help='BoxSize in Mpc/h')
parser.add_argument('dk_in_kf',type=float,
                    help='dk for power spectrum binning, as a multiple of k_f')
parser.add_argument('kmin_in_kf',type=float,
                    help='lower edge of lowest power spectrum bin, as a multiple of k_f')

# Parse arguments
args = parser.parse_args()
snap_prefix = args.snap_prefix
out_json_file = args.out_json_file
out_dat_file = args.out_dat_file
nmesh = args.nmesh
box_size = args.boxsize
dk_in_kf = args.dk_in_kf
kmin_in_kf = args.kmin_in_kf

comm = CurrentMPIComm.get()
    
# Set start time, and print first status message
start_time = time.time()
print_status(comm,start_time,'Starting measurement')

# Import particle positions
in_cat = Gadget1Catalog(snap_prefix+'*')
    
# Create mesh for density grid
in_mesh = in_cat.to_mesh(Nmesh=nmesh,BoxSize=box_size,
                         window='cic',compensated=True)

# Save density field, if desired (for debugging)
# in_mesh.save(out_dat_file+'_FIELD',dataset='Field',mode='real')
 
print_status(comm,start_time,'Created mesh for density grid')

# Set binning parameters, based on multiples of k_f
kf = 2*np.pi/box_size
dk = dk_in_kf*kf
kmin = kmin_in_kf*kf

# Measure power spectrum using nbodykit's built-in routines
print_status(comm,start_time,'About to measure power spectrum')
r = nbk.FFTPower(in_mesh, mode='1d', dk=dk, kmin=kmin)
Pk = r.power
print_status(comm,start_time,'Measured power spectrum')

# Save power spectrum to json file
if comm.rank == 0:
    r.save(out_json_file)
print_status(comm,start_time,'Saved power spectrum to json file')

# Also save power spectrum to plain text file
if comm.rank == 0:
    np.savetxt(out_dat_file,np.real(np.transpose([Pk['k'],Pk['power'],Pk['modes']])),
                header='k_mean_over_bin [h Mpc^-1] P_bin [h^-3 Mpc^3] N_modes')
    
print_status(comm,start_time,'Saved power spectrum to plain text file')
