"""Convert BAHAMAS HDF5 particle snapshot to density grid in bigfile format.

Output format is a grid of 1+\delta(\vec{x}), *NOT* corrected for the CIC mass
assignment (this needs to be corrected separately before any statistics are measured).
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py
import time
import warnings

import dask.array as da
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
                    help='input BAHAMAS snapshot (up to final dot in filenames)')
parser.add_argument('out_file',
                    help='output bigfile for density grid')
parser.add_argument('sim_type',
                    choices=['DM','DM-2fluid','hydro','baryons'],
                    help='type of density grid to make')
parser.add_argument('nmesh',
                    help='Nmesh for output grid')
parser.add_argument('boxsize',type=float,
                    help='BoxSize in Mpc/h')
parser.add_argument('--mass_moments_only',action='store_true',default=False,
                    help='Only compute <m>, <m^2>, <m^3> over all sim_type particles')

# Parse arguments
args = parser.parse_args()
snap_prefix = args.snap_prefix
out_file = args.out_file
sim_type = args.sim_type
nmesh = args.nmesh
boxsize = args.boxsize
do_mass_moments_only = args.mass_moments_only

comm = CurrentMPIComm.get()
    
# Set start time, and print first status message
start_time = time.time()
print_status(comm,start_time,'Starting conversion')

# If reading DM-only simulation, or making delta_DM grid from hydro snapshot:
if sim_type == 'DM':
    
    # Import particle positions (PartType1 = DM particles)
    in_cat = HDFCatalog(snap_prefix+'*',dataset='PartType1')
    
    # Create mesh for density grid
    in_mesh = in_cat.to_mesh(Nmesh=nmesh,BoxSize=boxsize,
                             resampler='cic',compensated=False,
                             position='Coordinates')

# If reading DM-only simulations with 2 DM fluids:
elif sim_type == 'DM-2fluid':
    
    # For each particle type (DM PartType1, DM PartType3):
    in_cat = {}
    for n in [1,3]:
        # Read in positions and masses of each particle type
        in_cat[n] = HDFCatalog(snap_prefix+'*',dataset='PartType'+str(n))
        
    # DM particles don't have a stored mass, so we create that column
    # using the values from the MassTable in the HDF5 file header
    f = h5py.File(snap_prefix+'0.hdf5')
    mass_table = f['Header'].attrs['MassTable']
    in_cat[1]['Masses'] = mass_table[1]
    in_cat[3]['Masses'] = mass_table[3]
    print_status(comm,start_time,'DM particle masses: %e %e' % \
                     (mass_table[1],mass_table[3]))
    
    # Make a combined catalog of all particles and masses
    comb_cat = nbk.MultipleSpeciesCatalog(['DM1','DM3'],in_cat[1],in_cat[3])

    # Create mesh for density grid
    in_mesh = comb_cat.to_mesh(Nmesh=nmesh,BoxSize=boxsize,compensated=False, 
                               resampler='cic',position='Coordinates',weight='Masses')
    
# If reading full hydro simulation or making delta_b grid:
elif sim_type == 'hydro' or sim_type == 'baryons':

    # For each particle type (gas, DM, stars, BHs):
    in_cat = {}
    for n in [0,1,4,5]:
        
        # Read in positions and masses of each particle type
        in_cat[n] = HDFCatalog(snap_prefix+'*',dataset='PartType'+str(n))
        
    # DM particles don't have a stored mass, so we create that column
    # using the values from the MassTable in the HDF5 file header
    f = h5py.File(snap_prefix+'0.hdf5')
    in_cat[1]['Mass'] = f['Header'].attrs['MassTable'][1]
    print_status(comm,start_time,'DM particle mass: %e' % \
                     f['Header'].attrs['MassTable'][1])
    
    # We also need to reassign the BH_Mass column to the Mass column for PartType5,
    # since Ian says we need to use BH_Mass for the black holes
    in_cat[5]['Mass'] = in_cat[5]['BH_Mass']
    
    # Make a combined catalog of all particles and masses
    if sim_type == 'hydro':
        part_types = ['gas','DM','stars','BHs']
        comb_cat = nbk.MultipleSpeciesCatalog(part_types,
                                              in_cat[0],in_cat[1],in_cat[4],in_cat[5])
    elif sim_type == 'baryons':
        part_types = ['gas','stars','BHs']
        comb_cat = nbk.MultipleSpeciesCatalog(part_types,
                                              in_cat[0],in_cat[4],in_cat[5])

    if do_mass_moments_only:
        # Compute moments of particle mass distribution. I don't know a one-line
        # command for taking the mean over multiple fields, so I do it manually.
        m_len = {pt: comb_cat.compute(da.count_nonzero(comb_cat[pt + '/Mass'])) for pt in part_types}
        m_sum = {pt: comb_cat.compute(comb_cat[pt + '/Mass'].sum()) for pt in part_types}
        m2_arr = {pt: comb_cat[pt + '/Mass']**2. for pt in part_types}
        m2_sum = {pt: comb_cat.compute(m2_arr[pt].sum()) for pt in part_types}
        m3_arr = {pt: comb_cat[pt + '/Mass']**3. for pt in part_types}
        m3_sum = {pt: comb_cat.compute(m3_arr[pt].sum()) for pt in part_types}
        
        m_len_all = np.sum(m_len.values())
        m_mean = np.sum(m_sum.values())/m_len_all
        m2_mean = np.sum(m2_sum.values())/m_len_all
        m3_mean = np.sum(m3_sum.values())/m_len_all

        print_status(comm,start_time,'Mean particle mass: %g' % m_mean)
        print_status(comm,start_time,'Mean squared particle mass: %g' % m2_mean)
        print_status(comm,start_time,'Mean cubed particle mass: %g' % m3_mean)
    else:
        # Create mesh for density grid
        in_mesh = comb_cat.to_mesh(Nmesh=nmesh,BoxSize=boxsize,compensated=False, 
                                   resampler='cic',position='Coordinates',weight='Mass')

# Paint grid to mesh, and save as bigfile
if not do_mass_moments_only:
    print_status(comm,start_time,'Created mesh for density grid')
    print_status(comm,start_time,'Starting paint and file output')

    # nbodykit may issue warnings here, because saving a mesh derived from a
    # MultipleSpeciesCatalog will throw away information related to the individual
    # species. As a blunt way to deal with this, we suppress warnings when saving
    # the mesh to disk
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        in_mesh.save(out_file,dataset='Field',mode='real')

    print_status(comm,start_time,'Wrote grid to %s' % out_file)
