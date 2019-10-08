"""Measure matter power spectrum from density grid in bigfile format.

The -cic option controls whether the CIC window function deconvolution is
applied to the grid before the measurement is made. The script allows for
the original mass assignment to have occured at a different grid resolution
than the input grid: nmesh is the actual input grid resolution, while nmesh_cic
is the resolution at which the mass assignment took place. For example,
for the grids used for the bispectrum paper, we originally made 2048^3 grids and
then downsampled them to 1024^3, so the 2048^3 CIC window function needs to be
deconvolved.
"""

import numpy as np
import argparse
import h5py
import time

import nbodykit.lab as nbk
from nbodykit.source.catalog.file import Gadget1Catalog,HDFCatalog
from nbodykit import CurrentMPIComm,setup_logging
from nbodykit.source.mesh.catalog import CompensateCICShotnoise

# Routine to print script status to command line, with elapsed time
def print_status(comm,start_time,message):
    if comm.rank == 0:
        elapsed_time = time.time() - start_time
        print('%d\ts: %s' % (elapsed_time,message))

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('in_file',
                    help='input bigfile')
parser.add_argument('out_json_file',
                    help='output power spectrum file, in nbodykit json format')
parser.add_argument('out_dat_file',
                    help='output power spectrum file, in plain text format')
parser.add_argument('nmesh',type=int,
                    help='Nmesh for FFTs')
parser.add_argument('nmesh_cic',type=int,
                    help='Nmesh to be used for CIC correction')
parser.add_argument('box_size',type=float,
                    help='box size in Mpc/h')
parser.add_argument('dk_in_kf',type=float,
                    help='dk for power spectrum binning, as a multiple of k_f')
parser.add_argument('kmin_in_kf',type=float,
                    help='lower edge of lowest power spectrum bin, as a multiple of k_f')
parser.add_argument('--cross',
                    help='second input bigfile, for computing cross spectrum')
parser.add_argument('-cic','--cic_correction',action='store_true',
                    help='whether to apply extra CIC correction to input snapshot')

# Parse arguments
args = parser.parse_args()
in_file = args.in_file
out_json_file = args.out_json_file
out_dat_file = args.out_dat_file
Nmesh = args.nmesh
NmeshCIC = args.nmesh_cic
box_size = args.box_size
dk_mult = args.dk_in_kf
kmin_mult = args.kmin_in_kf
cross_in_file = args.cross
do_cic = args.cic_correction


def CompensateCICShotnoiseNgrid(w,v):
    """Compensate for CIC window function at original resolution.
    
    CIC window function to be applied to grid with resolution Nmesh, but corrected
    to correspond to a grid with resolution NmeshCIC. We need to do this because 
    the original density grid may have had a higher resolution than the (possibly
    subsampled) grid we operate on. See nbodykit.source.mesh.catalog.py for
    original syntax, and Jing et al. 2005 Eq. 20 for formula reference.
    """
    for i in range(3):
        wi = w[i]
        v = v / (1 - 2. / 3 * np.sin(0.5 * wi * Nmesh / NmeshCIC) ** 2) ** 0.5
    return v


comm = CurrentMPIComm.get()

# Set start time, and print first status message
start_time = time.time()
print_status(comm,start_time,'Starting script')

# Import mesh from file
mesh = nbk.BigFileMesh(in_file,'Field')
print_status(comm,start_time,'Read input grid from %s' % in_file)
if (cross_in_file is not None):
    cross_mesh = nbk.BigFileMesh(cross_in_file,'Field')
    print_status(comm,start_time,'Read input grid from %s' % cross_in_file)

if do_cic:
    # Apply CIC window compensation to mesh(es)
    mesh = mesh.apply(CompensateCICShotnoiseNgrid,kind='circular',mode='complex')
    print_status(comm,start_time,'Applied CIC compensation to mesh')
    if (cross_in_file is not None):
        cross_mesh = cross_mesh.apply(CompensateCICShotnoiseNgrid,kind='circular',mode='complex')
        print_status(comm,start_time,'Applied CIC compensation to cross_mesh')

# Set binning parameters, based on multiples of k_f
kf = 2*np.pi/box_size
dk = dk_mult*kf
kmin = kmin_mult*kf

# Set output string for measured spectrum
if (cross_in_file is not None):
    spec_string = 'cross spectrum'
else:
    spec_string = 'power spectrum'

# Measure power spectrum using nbodykit's built-in routines
print_status(comm,start_time,'About to measure '+spec_string)
if (cross_in_file is not None):
    r = nbk.FFTPower(mesh, second=cross_mesh, mode='1d', dk=dk, kmin=kmin)
else:
    r = nbk.FFTPower(mesh, mode='1d', dk=dk, kmin=kmin)
Pk = r.power
print_status(comm,start_time,'Measured '+spec_string)

# Save power spectrum to json file
if comm.rank == 0:
    r.save(out_json_file)
print_status(comm,start_time,'Saved '+spec_string+' to json file')

# Also save power spectrum to plain text file
if comm.rank == 0:
    np.savetxt(out_dat_file,np.real(np.transpose([Pk['k'],Pk['power'],Pk['modes']])),
               header='k_mean_over_bin [h Mpc^-1] P_bin [h^-3 Mpc^3] N_modes')
    
print_status(comm,start_time,'Saved '+spec_string+' to plain text file')
