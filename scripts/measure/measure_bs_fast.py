"""Measure auto or cross bispectrum from input snapshot(s) with a specified binning scheme.

The -cic option controls whether the CIC window function deconvolution is
applied to the grid before the measurement is made. The script allows for
the original mass assignment to have occured at a different grid resolution
than the input grid: nmesh is the actual input grid resolution, while nmesh_cic
is the resolution at which the mass assignment took place. For example,
for the grids used for the bispectrum paper, we originally made 2048^3 grids and
then downsampled them to 1024^3, so the 2048^3 CIC window function needs to be
deconvolved.

Binning can be specified in two ways:
    - a k_bin_file, containing an array of shape (6,N_tri), where rows of this array
    specify k1_min,k1_max,k2_min,k2_max,k3_min,k3_max for a set of triangle bins
    - the --k_bin_info argument, with numerical values for k_min, k_max, and dk
        - in addition to this, the --high_k_bin_info argument can supply a number of
        low-k bins with width dk, and a second value of dk for bins at higher k.
        EXAMPLE: passing "--k_bin_info 0.05 1. 0.1 --high_k_bin_info 10 0.2"
        specifies a binning scheme where the 10 lowest-k bins have width 0.1,
        and the remaining higher-k bins have width 0.2
        
The --meas_type argument specifies whether to measure the mean k1,k2,k3 values 
and number of triangles in each bin ('grid info') or bispectrum value 
(unnormalized by the number of triangles) ('unnorm_b_value').
"""

import numpy as np
import argparse
import time

import nbodykit.lab as nbk
from nbodykit.source.catalog.file import Gadget1Catalog,HDFCatalog
from nbodykit import CurrentMPIComm
from nbodykit.source.mesh.catalog import CompensateCICShotnoise

import bskit

# Routine to print script status to command line, with elapsed time
def print_status(comm,start_time,message):
    if comm.rank == 0:
        elapsed_time = time.time() - start_time
        print('%d\ts: %s' % (elapsed_time,message))


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

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('snap_type',
                    choices=['gadget','illustris','bigfile_grid'],
                    help='type of snapshot')
parser.add_argument('snap_prefix1',
                    help='complete file path to snapshot files for mesh 1, including file prefix')
parser.add_argument('nmesh',type=int,
                    help='Nmesh for FFTs')
parser.add_argument('nmesh_cic',type=int,
                    help='Nmesh to be used for CIC correction')
parser.add_argument('boxsize',type=float,
                    help='BoxSize in Mpc/h')
parser.add_argument('start_bin_index',type=int,
                    help='first bin index to use')
parser.add_argument('end_bin_index',type=int,
                    help='last bin index to use')
parser.add_argument('out_file',
                    help='file path to output file')
parser.add_argument('--k_bin_file',
                    help='complete file path to file containing k bin definitions' \
                           + '(can specify instead of kmin,kmax,dk values)')
parser.add_argument('--k_bin_info',type=float,nargs=3,
                    help='minimum and maximum k values for bispectrum triangles, ' \
                    + 'and bin width, in h/Mpc')
parser.add_argument('--high_k_bin_info',type=float,nargs=2,
                    help='number of bins with smaller dk, and dk for higher bins')
parser.add_argument('-pu','--pos_units_mpcoverh',type=float,
                    help='units of particle positions in file, in Mpc/h')
parser.add_argument('-cic','--cic_correction',action='store_true',
                    help='whether to apply extra CIC correction to input snapshot')
parser.add_argument('--snap_prefix2',
                    help='complete file path to snapshot files for mesh 2, including file prefix')
parser.add_argument('--snap_prefix3',
                    help='complete file path to snapshot files for mesh 2, including file prefix')
parser.add_argument('--triangle_type',
                    choices=['all','equilateral','squeezed','isosceles'],
                    help='type of triangles to consider: all, equilateral, squeezed, or isosceles')
parser.add_argument('--squeezed_bin_index',type=int,
                    help='k bin index of k_long bin for squeezed triangles')
parser.add_argument('--isos_mult',type=float,
                    help='multiplier for isosceles triangles, such that k_S = mult * k_L')
parser.add_argument('--isos_tol',type=float,
                    help='tolerance for relationship between k_S and k_L to count as isosceles')
parser.add_argument('--meas_type',
                    choices=['grid_info','unnorm_b_value'],default='unnorm_b_value',
                    help='type of measurement to do (grid_info or unnorm_b_value')


snap_prefix = {}

# Parse arguments
args = parser.parse_args()
snap_type = args.snap_type
snap_prefix[0] = args.snap_prefix1
snap_prefix[1] = args.snap_prefix2
snap_prefix[2] = args.snap_prefix3
Nmesh = args.nmesh
NmeshCIC = args.nmesh_cic
BoxSize = args.boxsize
start_i = args.start_bin_index
end_i = args.end_bin_index
out_file = args.out_file
do_cic = args.cic_correction

triangle_type = 'all'
if args.triangle_type is not None:
    triangle_type = args.triangle_type

for_grid_info_only = False
meas_type = args.meas_type
if meas_type == 'grid_info':
    for_grid_info_only = True

pos_fac = 1.
if args.pos_units_mpcoverh is not None:
    pos_fac = args.pos_units_mpcoverh

if (args.k_bin_file and args.k_bin_info) or (not args.k_bin_file and not args.k_bin_info):
    raise RuntimeError('Must specify either a k-bin file or kmin/kmax/dk!')

# Set some default parameters
num_low_k_bins = 0
dk_high = -1.
squeezed_bin_index = 0
isos_mult = 0
isos_tol = 0.1
    
comm = CurrentMPIComm.get()

# Set start time, and print first status message
start_time = time.time()
print_status(comm,start_time,'Starting script')
    
if args.k_bin_file is not None:
    # If k bin file has been specified, load it
    k_bin_file = args.k_bin_file
    b_k_bins = np.loadtxt(k_bin_file).T
else:
    # Otherwise, set binning scheme based on other inputs.
    kmin = args.k_bin_info[0]
    kmax = args.k_bin_info[1]
    dk = args.k_bin_info[2]

    # If widths should be different for low-k and high-k bins, 
    # take this information from the command line: the number of low-k bins
    # to use, and the high-k bin width
    if args.high_k_bin_info is not None:
        num_low_k_bins = np.int(args.high_k_bin_info[0])
        dk_high = args.high_k_bin_info[1]
        
    if triangle_type == 'squeezed':
        # If we're only measuring squeezed triangles, set bin index of k_L (squeezed side)
        print_status(comm,start_time,'Squeezed triangles requested')
        if args.squeezed_bin_index is not None:
            squeezed_bin_index = args.squeezed_bin_index
            print_status(comm,start_time,'--- k_L bin index: %d' % squeezed_bin_index)
            
    elif triangle_type == 'isosceles':
        # Or, if we're only measuring isosceles triangles, get multiplier that defines
        # them, and tolerance for side to obey isosceles condition
        print_status(comm,start_time,'Isosceles triangles requested')
        if args.isos_mult is None:
            raise ValueError('Need to set isos_mult for isosceles triangles!')
        else:
            isos_mult = args.isos_mult
            if args.isos_tol is not None:
                isos_tol = args.isos_tol
            print_status(comm,start_time,'--- isos_mult: %f' % isos_mult)
            print_status(comm,start_time,'--- isos_tol:  %f' % isos_tol)

print_status(comm,start_time,'Starting input')

# Define empty dicts for catalogs and meshes
cat = {}
mesh = {}
for i in range(3):
    mesh[i] = None

# Figure out how many input snapshots there are
num_fields = 1
if snap_prefix[1] is not None:
    num_fields += 1
if snap_prefix[2] is not None:
    if snap_prefix[1] is None:
        raise ValueError('need second mesh if third mesh is input!')
    num_fields += 1

# Get particle catalogs, or meshes, if working with bigfiles
for i in range(num_fields):
    if snap_type == 'illustris':
        cat[i] = HDFCatalog(snap_prefix[i]+'*')
        mesh[i] = cat[i].to_mesh(Nmesh=Nmesh,BoxSize=BoxSize/pos_fac,
                                     window='cic',compensated=True,
                                     position='PartType1/Coordinates')
    elif snap_type == 'gadget':
        cat[i] = Gadget1Catalog(snap_prefix[i]+'*')
        mesh[i] = cat[i].to_mesh(Nmesh=Nmesh,BoxSize=BoxSize,
                                     window='cic',compensated=True)
    elif snap_type == 'bigfile_grid':
        mesh[i] = nbk.BigFileMesh(snap_prefix[i],'Field')
        if do_cic:
            # Apply CIC window compensation to mesh
            mesh[i] = mesh[i].apply(CompensateCICShotnoiseNgrid,kind='circular',mode='complex')
            print_status(comm,start_time,'Applied CIC compensation to mesh %d' % i)

print_status(comm,start_time,'Finished input')
print_status(comm,start_time,'Starting initial paints')

if args.k_bin_file:
    # If k bin file has been specified, initialize FFTBispectrum object with the
    # corresponding binning scheme
    fftb = bskit.FFTBispectrum(mesh[0],Nmesh=Nmesh,BoxSize=np.ones(3)*BoxSize,
                               k_edges=b_k_bins,pos_units_mpcoverh=pos_fac,
                               second=mesh[1],third=mesh[2],for_grid_info_only=for_grid_info_only)
else:
    # Otherwise, initialize FFTBispectrum with binning defined by the other inputs,
    # and restricted to certain triangles if desired
    fftb = bskit.FFTBispectrum(mesh[0],Nmesh=Nmesh,BoxSize=np.ones(3)*BoxSize,
                               dk=dk,kmin=kmin,kmax=kmax,pos_units_mpcoverh=pos_fac,
                               second=mesh[1],third=mesh[2],
                              num_lowk_bins=num_low_k_bins,dk_high=dk_high,
                              triangle_type=triangle_type,squeezed_bin_index=squeezed_bin_index,
                              isos_mult=isos_mult,isos_tol=isos_tol,
                              for_grid_info_only=for_grid_info_only)
        
print_status(comm,start_time,'Finished initial paints')
num_k_bins = len(bskit.generate_bin_edge_list(fftb.attrs['kmin'],fftb.attrs['kmax'],
                                              fftb.attrs['dk'],fftb.attrs['num_lowk_bins'],
                                              fftb.attrs['dk_high']))
print_status(comm,start_time,'Number of linear k bins: %d' % num_k_bins)

if meas_type == 'grid_info':
    # If measuring grid (binning) info only
    print_status(comm,start_time,'Starting FFT binning info computation')
    fftb.measure_gridinfo_faster(imin=start_i,imax=end_i,out_file=out_file,
                                   verbose=0)
else:
    # If measuring bispectrum
    print_status(comm,start_time,'Starting FFT bispectrum')
    fftb.measure_bispectrum_faster(imin=start_i,imax=end_i,out_file=out_file,
                                   verbose=0)

print_status(comm,start_time,'Finished FFT bispectrum')
