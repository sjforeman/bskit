"""Convert Illustris (or TNG) HDF5 particle snapshot to density grid in bigfile format.

Output format is a grid of 1+\delta(\vec{x}), *NOT* corrected for the CIC mass
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
def print_status(comm,start_time,message):
    if comm.rank == 0:
        elapsed_time = time.time() - start_time
        print('%d\ts: %s' % (elapsed_time,message))

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('snap_prefix',
                    help='input Illustris snapshot (up to final dot in filenames)')
parser.add_argument('out_file',
                    help='output bigfile for density grid')
parser.add_argument('sim_type',
                    choices=['DM','hydro','baryons'],
                    help='type of density grid to make')
parser.add_argument('nmesh',
                    help='Nmesh for output grid')
parser.add_argument('boxsize',type=float,
                    help='BoxSize in Mpc/h')
parser.add_argument('-mDM','--mass_DM',type=float,
                    help='DM particle mass, in 10^10 h^-1 Msun')

# Parse arguments
args = parser.parse_args()
snap_prefix = args.snap_prefix
out_file = args.out_file
sim_type = args.sim_type
nmesh = args.nmesh
boxsize = args.boxsize
mass_DM = args.mass_DM

if sim_type == 'hydro' and not mass_DM:
    raise RuntimeError('Must specify a DM particle mass for a full hydro snapshot!')

comm = CurrentMPIComm.get()
    
# Set start time, and print first status message
start_time = time.time()
print_status(comm,start_time,'Starting conversion')

# If reading DM-only simulation, or making delta_DM grid from hydro snapshot:
if sim_type == 'DM':
    
    # Import particle positions (PartType1 = DM particles)
    in_cat = HDFCatalog(snap_prefix+'*',dataset='PartType1')
    
    # Convert particle positions from kpc/h to Mpc/h
    in_cat['Coordinates'] /= 1e3
    
    # Create mesh for density grid
    in_mesh = in_cat.to_mesh(Nmesh=nmesh,BoxSize=boxsize,
                             resampler='cic',compensated=False,
                             position='Coordinates')
    
# If reading full hydro simulation or making delta_b grid:
elif sim_type == 'hydro' or sim_type == 'baryons':
    
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
    if sim_type == 'hydro':
        comb_cat = nbk.MultipleSpeciesCatalog(['gas','DM','stars','BHs'],
                                              in_cat[0],in_cat[1],in_cat[4],in_cat[5])
    elif sim_type == 'baryons':
        comb_cat = nbk.MultipleSpeciesCatalog(['gas','stars','BHs'],
                                              in_cat[0],in_cat[4],in_cat[5])

    # Create mesh for density grid, weighting particles by mass
    in_mesh = comb_cat.to_mesh(Nmesh=nmesh,BoxSize=boxsize,compensated=False, 
                               resampler='cic',position='Coordinates',weight='Masses')

print_status(comm,start_time,'Created mesh for density grid')

# Paint grid to mesh, and save as bigfile.
# nbodykit may issue warnings here, because saving a mesh derived from a
# MultipleSpeciesCatalog will throw away information related to the individual
# species. As a blunt way to deal with this, we suppress warnings when saving
# the mesh to disk
print_status(comm,start_time,'Starting paint and file output')
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    in_mesh.save(out_file,dataset='Field',mode='real')
print_status(comm,start_time,'Wrote grid to %s' % out_file)
