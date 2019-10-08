"""Convert density grid in bigfile format to density grid in numpy format.

There are no command-line arguments for this - the input and output
directories and filenames need to be specified directly in the script.
"""

import numpy as np

import nbodykit.lab as nbk
from nbodykit import CurrentMPIComm

# Set directories for input and output files
in_dir = '/mnt/raid-cita/sforeman/illustrisTNG/sjf_grids/'
out_dir = '/mnt/raid-cita/sforeman/illustrisTNG/will_grids/'

# Routine to convert bigfile to numpy file
def convert_bigfile_to_npy(in_file,out_file):
    orig_mesh = nbk.BigFileMesh(in_file,'Field')
    orig_field = orig_mesh.compute(mode='real')
    print('Finished reading bigfile grid')
    # Note that Yu recommends not using preview for exporting computed grids
    # (see https://github.com/bccp/nbodykit/issues/599), but it actually
    # worked fine for our numerical tests (i.e. bispectra measured from 
    # the bigfile and numpy files agreed at close to numerical precision).
    orig_npy_array = orig_field.preview()
    np.save(out_file,orig_npy_array)
    print('Converted grid %s to npy file %s' % (in_file,out_file))
    
# Loop over files to convert - can loop over multiple Nmesh values or filenames
# by manually editing, as desired
for Nmesh in [2048]:
    convert_bigfile_to_npy(in_dir + 'TNG_DM300-1_delta_m_%d_z0.0.bigfile' % Nmesh,
                           out_dir + 'TNG_DM300-1_delta_m_%d_z0.0.npy' % Nmesh)
    convert_bigfile_to_npy(in_dir + 'TNG300-1_delta_m_%d_z0.0.bigfile' % Nmesh,
                           out_dir + 'TNG300-1_delta_m_%d_z0.0.npy' % Nmesh)
