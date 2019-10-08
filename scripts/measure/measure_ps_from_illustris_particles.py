"""Measure matter power spectrum from particle snapshot in Illustris format.

Assumes coordinates are stored as kpc/h, and internally converts to Mpc/h
before computing power spectrum.
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
                    help='input Illustris snapshot (up to final dot in filenames)')
parser.add_argument('out_json_file',
                    help='output power spectrum file, in nbodykit json format')
parser.add_argument('out_dat_file',
                    help='output power spectrum file, in plain text format')
parser.add_argument('sim_type',
                    choices=['DM','hydro'],
                    help='type of snapshot')
parser.add_argument('nmesh',
                    help='Nmesh for measurement')
parser.add_argument('boxsize',type=float,
                    help='BoxSize in Mpc/h')
parser.add_argument('-mDM','--mass_DM',type=float,default=None,
                    help='DM particle mass, in 10^10 h^-1 Msun')
parser.add_argument('dk_in_kf',type=float,
                    help='dk for power spectrum binning, as a multiple of k_f')
parser.add_argument('kmin_in_kf',type=float,
                    help='lower edge of lowest power spectrum bin, as a multiple of k_f')

# Parse arguments
args = parser.parse_args()
snap_prefix = args.snap_prefix
out_json_file = args.out_json_file
out_dat_file = args.out_dat_file
sim_type = args.sim_type
nmesh = args.nmesh
box_size = args.boxsize
mass_DM = args.mass_DM
dk_in_kf = args.dk_in_kf
kmin_in_kf = args.kmin_in_kf

if sim_type == 'hydro' and mass_DM is None:
    raise RuntimeError('Must specify a DM particle mass for a full hydro snapshot!')

comm = CurrentMPIComm.get()
    
# Set start time, and print first status message
start_time = time.time()
print_status(comm,start_time,'Starting measurement')

# If reading DM-only simulation:
if sim_type == 'DM':
    
    # Import particle positions (PartType1 = DM particles)
    in_cat = HDFCatalog(snap_prefix+'*',dataset='PartType1')
    
    # Convert particle positions from kpc/h to Mpc/h
    in_cat['Coordinates'] /= 1e3
    
    # Create mesh for density grid
    in_mesh = in_cat.to_mesh(Nmesh=nmesh,BoxSize=box_size,
                             resampler='cic',compensated=True,
                             position='Coordinates')
    
    # Save density field in Fourier space, if desired (for debugging)
#     in_mesh.save(out_dat_file+'_FIELD',dataset='Field',mode='complex')
    
# If reading full hydro simulation:
elif sim_type == 'hydro':
    
    # For each particle type (gas, DM, stars, BHs):
    in_cat = {}
    for n in [0,1,4,5]:
        
        # Read in positions and masses of each particle type
        in_cat[n] = HDFCatalog(snap_prefix+'*',dataset='PartType'+str(n))
        
        # Convert particle positions from kpc/h to Mpc/h
        in_cat[n]['Coordinates'] /= 1e3
        
    # DM particles don't have a stored mass, so we create that column
    in_cat[1]['Masses'] = mass_DM
    
    # Make a combined catalog of all particles and masses
    comb_cat = nbk.MultipleSpeciesCatalog(['gas','DM','stars','BHs'],
                                          in_cat[0],in_cat[1],in_cat[4],in_cat[5])

    # Create mesh for density grid
    in_mesh = comb_cat.to_mesh(Nmesh=nmesh,BoxSize=box_size,compensated=True, 
                               resampler='cic',position='Coordinates',weight='Masses')

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
