"""Measure matter power spectrum from subregions of density grid in bigfile format.

The -cic option controls whether the CIC window function deconvolution is
applied to the grid before the measurement is made. The script allows for
the original mass assignment to have occured at a different grid resolution
than the input grid: nmesh is the actual input grid resolution, while nmesh_cic
is the resolution at which the mass assignment took place. For example,
for the grids used for the bispectrum paper, we originally made 2048^3 grids and
then downsampled them to 1024^3, so the 2048^3 CIC window function needs to be
deconvolved.

Note that subbox binning needs to be specified in terms of k_f of the full box.
"""

import numpy as np
import argparse
import h5py
import time

import nbodykit.lab as nbk
from nbodykit.source.catalog.file import Gadget1Catalog,HDFCatalog
from nbodykit import CurrentMPIComm,setup_logging
from nbodykit.source.mesh.catalog import CompensateCICShotnoise

import bskit

# Routine to print script status to command line, with elapsed time
def print_status(comm,start_time,message):
    if comm.rank == 0:
        elapsed_time = time.time() - start_time
        print('%d\ts: %s' % (elapsed_time,message))

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('in_file',
                    help='input bigfile')
parser.add_argument('out_dat_file',
                    help='output power spectrum file, in plain text format')
parser.add_argument('nmesh',type=int,
                    help='Nmesh for FFTs of full box')
parser.add_argument('nmesh_cic',type=int,
                    help='Nmesh to be used for CIC correction')
parser.add_argument('box_size',type=float,
                    help='box size in Mpc/h')
parser.add_argument('dk_in_kf',type=float,
                    help='dk for power spectrum binning, as a multiple of k_f')
parser.add_argument('kmin_in_kf',type=float,
                    help='lower edge of lowest power spectrum bin, as a multiple of k_f')
parser.add_argument('nsub_per_side',type=int,
                    help='cube root of number of subboxes to take')
parser.add_argument('--cross',
                    help='second input bigfile, for computing cross spectrum')
parser.add_argument('-cic','--cic_correction',action='store_true',
                    help='whether to apply extra CIC correction to input snapshot')

# Parse arguments
args = parser.parse_args()
in_file = args.in_file
out_dat_file = args.out_dat_file
Nmesh = args.nmesh
NmeshCIC = args.nmesh_cic
box_size = args.box_size
dk_mult = args.dk_in_kf
kmin_mult = args.kmin_in_kf
nsub_per_side = args.nsub_per_side
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
    # Apply CIC window compensation to mesh
    mesh = mesh.apply(CompensateCICShotnoiseNgrid,kind='circular',mode='complex')
    print_status(comm,start_time,'Applied CIC compensation to mesh')
    if (cross_in_file is not None):
        cross_mesh = cross_mesh.apply(CompensateCICShotnoiseNgrid,kind='circular',mode='complex')
        print_status(comm,start_time,'Applied CIC compensation to cross_mesh')

# Set binning parameters for subbox, based on multiples of k_f in full box,
# but adding an extra dk to kmin, since the longest mode in the full box
# will not be present in the subbox
kf = 2*np.pi/box_size
dk = dk_mult*kf
kmin = (kmin_mult+1)*kf

# Set output string for measured spectrum
if (cross_in_file is not None):
    spec_string = 'cross spectrum'
else:
    spec_string = 'power spectrum'

# Paint mesh (and cross-mesh, if present) to field(s)
field = mesh.compute(mode='real')
if (cross_in_file is not None):
    cross_field = cross_mesh.compute(mode='real')
    
# Construct string defining bounds for each subbox, to be included
# in header of output file
header_str = 'subbox_index k_mean_over_bin [h Mpc^-1] P_bin [h^-3 Mpc^3] N_modes\n'
for i in range(nsub_per_side):
    for j in range(nsub_per_side):
        for k in range(nsub_per_side):
            header_str += str(bskit.subbox_multiindex_to_index((i,j,k),nsub_per_side)) \
              + ': ' + str(np.array(((i,i+1),(j,j+1),(k,k+1)))*box_size/nsub_per_side) + '\n'

# Loop over subboxes
result = {}    
for i in range(nsub_per_side):
    for j in range(nsub_per_side):
        for k in range(nsub_per_side):
            # Convert indices representing subset of each coordinate axis
            # to single index for all Nsub^3 subboxes
            ind = bskit.subbox_multiindex_to_index((i,j,k),nsub_per_side)
            
            # Fetch mesh corresponding to subbox
            print_status(comm,start_time,'About to measure '+spec_string \
                         +' from subbox %d' % ind)
            sub_mesh = bskit.field_subbox_pm((i,j,k),nsub_per_side,field)
        
            # Measure subbox auto spectrum, or fetch cross subbox and measure cross spectrum
            if (cross_in_file is None):
                r = nbk.FFTPower(sub_mesh, mode='1d', dk=dk, kmin=kmin)
            else:
                sub_cross_mesh = bskit.field_subbox_pm((i,j,k),nsub_per_side,cross_field)
                r = nbk.FFTPower(sub_mesh, second=sub_cross_mesh, mode='1d', dk=dk, 
                                 kmin=kmin)
        
            # Add P(k) dict to results dict
            result[ind] = r.power
            print_status(comm,start_time,'Measured '+spec_string+' from subbox %d' % ind)
        
            # Save plain text file with all results so far (so that partial results
            # have been saved in case of crash or timeout)
            if comm.rank == 0:
                # Construct big array with all results, with first column as subbox index,
                # First, build a new dict with tables with (k,Pk,N_modes) columns
                out_array = {}
                for jj in range(ind+1):
                    out_array[jj] = np.transpose([np.ones(len(result[jj]['k']))*jj,
                                                 result[jj]['k'],
                                                 result[jj]['power'],
                                                 result[jj]['modes']])
           
                # If we've measured from more than one subbox, concatenate the results
                if ind == 0:
                    full_output = out_array[0]
                else:
                    full_output = np.concatenate([out_array[jj] for jj in range(ind+1)])
       
                # Save to text file, with header defining subbox bounds by index
                np.savetxt(out_dat_file,np.real(full_output),header=header_str,
                           fmt='%d %e %e %e')
                print_status(comm,start_time,'Saved '+spec_string+' to plain text file')

print_status(comm,start_time,'Finished all subbox measurements')