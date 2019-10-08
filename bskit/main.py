"""Bispectrum measurement module, building on nbodykit

Module for making power spectrum and bispectrum measurements with nbodykit
using FFTs, as described e.g. in Donghui Jeong's PhD thesis. The end goal
is to integrate this functionality into nbodykit itself, but for now,
I've collected the associated routines in this module.

Main documentation and example scripts can be found in the module's
repository on github:
    https://github.com/sjforeman/bskit
Also see that repository for a citation you are requested to include in
any published work that uses this code.

Simon Foreman, Perimeter Institute, 2019
"""

import numpy as np
import logging
from scipy.interpolate import interp1d
import itertools

import nbodykit.lab as nbk
import nbodykit.algorithms.fftpower as nbkfftpower
from nbodykit import setup_logging, style, CurrentMPIComm
from nbodykit.meshtools import SlabIterator
from pmesh.pm import ParticleMesh, RealField, ComplexField
from pmesh.window import Affine


def subbox_multiindex_to_index(multiindex,nsub_per_side):
    """Convert a subbox multi-index to a single index.
    
    Given a subbox multi-index of the form (a,b,c), return a single
    index corresponding to that, via i = a*N**2 + b*N + c
    where N = nsub_per_side.
    
    Parameters
    ----------
    multiindex : array_like
        tuple of form (a,b,c) that indexes subregion of each coordinate axis
    nsub_per_side : int
        cube root of number of subboxes in full box
    
    Returns
    -------
    int
        single index corresponding to subbox referred to by `multiindex`
    """
    assert(len(multiindex)==3)
    
    return int(multiindex[0]*nsub_per_side**2 + multiindex[1]*nsub_per_side + multiindex[2])


def subbox_index_to_multiindex(i,nsub_per_side):
    """Convert a subbox single index to a multi-index
    
    Given a subbox index, return the multi-index corresponding to that,
    which is (a,b,c) where i = a*N**2 + b*N + c and N = nsub_per_side.
    
    Parameters
    ----------
    i : int
        subbox index
    nsub_per_side : int
        cube root of number of subboxes in full box
    
    Returns
    -------
    array_like
        tuple of form (a,b,c) that indexes subregion of each coordinate axis
    """
    m = np.zeros(3)
    
    div_n2 = divmod(i,nsub_per_side**2)
    div_n1 = divmod(div_n2[1],nsub_per_side)
    
    m[0] = div_n2[0]
    m[1] = div_n1[0]
    m[2] = div_n1[1]

    return m


def field_subbox_pm(box_multiindex,
                    nsub_per_side,
                    source):
    """Construct a RealField corresponding to a subregion of an input RealField.
    
    The size of the subregion is Lsub = source.BoxSize/nsub_per_side, and the
    location of the subregion is specified by the coordinate origin of the
    subregion in units of Lsub by box_multiindex.
    
    For example, for box_multiindex=(1,0,0), the location of the box within
    the full box is ((Lsub,2Lsub),(0,Lsub),(0,Lsub)).
    
    Parameters
    ----------
    box_multiindex : array_like
        tuple of form (a,b,c) that indexes subregion of each coordinate axis
    nsub_per_side : int
        cube root of number of subboxes in full box
    source : RealField
        field within full box that we want to cut the subbox from
    
    Returns
    -------
    RealField:
        field within requested subbox
    """
    assert isinstance(source,RealField)

    # Define ParticleMesh and RealField with size and Nmesh of desired subregion
    sub_pm = ParticleMesh(BoxSize=source.BoxSize/nsub_per_side,
                          Nmesh=source.Nmesh/nsub_per_side)
    sub_field = RealField(sub_pm)
    
    # Get coordinate indices of the slab of the subbox grid stored on this task
    coords = sub_field.pm.mesh_coordinates(dtype='i4')

    # Find coordinate indices of origin of subregion coordinates within full region
    start_coord = box_multiindex*source.Nmesh/nsub_per_side
    
    # Define transformation between coordinates within subregion and coordinates
    # within full-region slab stored on this task. Specifically, the subregion
    # coordinates are shifted to the specific subregion within the full-region
    # (by adding start_coord), and then shifted back such that their origin
    # coincides with the origin of the full-region slab (by subtracting source.start)
    transform = Affine(sub_field.ndim,
                translate=start_coord-source.start,
                scale=1.0,
                period=source.Nmesh)

    # Find domain decomposition of desired subregion within full region,
    # to know which tasks ask for different coordinates. We need to manually add
    # start_coord, because the transform argument passed to the decompose routine
    # only shifts based on each slab's origin (it ignores the "translate" argument
    # specified above)
    layout_src = source.pm.decompose(coords+start_coord, smoothing=0, transform=transform)

    # Find values of source grid in subregion, as 1d list matched to 1d coords list
    vals = source.readout(coords, resampler='nnb', transform=transform, layout=layout_src)
    
    # Make new RealField corresponding to subregion. No "layout" argument is needed
    # because the coords list was fetched directly from the subbox slab on this task
    # i.e. each task will only paint its own values to the field
    return sub_field.pm.paint(coords, mass=vals, resampler='nnb', 
                               transform=sub_field.pm.affine_grid )


def field_subbox_slice(field,bounds):
    """Construct a mesh corresponding to a subregion of an input RealField.
    
    The bounds argument is a list of (min,max) coordinate values that define
    the subregion.
    
    Parameters
    ----------
    field : RealField
        field within full box that we want to cut the subbox from
    bounds : array_like
        array with shape(ndim,2) and ndim = 2 or 3, specifying coordinate bounds of subregion
        
    Returns
    -------
    ArrayMesh
        mesh corresponding to field within subregion
    """
    
    comm = CurrentMPIComm.get()
    
    # Check type and dimensions of input array
    if not isinstance(field, RealField):
        raise TypeError("input field must be RealField")
    ndim = len(field.x)
    if ndim not in [2,3]:
        raise NotImplementedError("field_subbox_slice can only be used on 3D or 2D arrays")
    
    # For each coordinate axis, construct array of indices corresponding to
    # coordinates that fall within the desired subbox bounds
    coords = {}
    ci = {}
    for i in range(ndim):
        coords[i] = field.x[i].copy()
        coords[i] += field.BoxSize[i]*(coords[i]<0)
        ci[i] = np.where((coords[i]>=bounds[i][0]) & (coords[i]<bounds[i][1]))[i]
        
    # Get side lengths of subbox
    box_size = [bounds[i][1]-bounds[i][0] for i in range(ndim)]
    
    # Return ArrayMesh corresponding to field within subbox
    if ndim == 2:
        return nbk.ArrayMesh(field[np.ix_(ci[0],ci[1])],BoxSize=box_size)
    elif ndim == 3:
        return nbk.ArrayMesh(field[np.ix_(ci[0],ci[1],ci[2])],BoxSize=box_size)
    
    
# def SubboxBoundIterator(Lbox,N_per_axis):
#     '''
#     Iterator for looping through N_per_axis^3 subregions of a cubical box 
#     with side length Lbox. The syntax for using is:
    
#     for b in SubboxBoundIterator(1000,2):
#         do_something(b)
        
#     The format of b is ((xmin,xmax),(ymin,ymax),(zmin,zmax))
#     '''
#     dL = Lbox/N_per_axis
    
#     for i in range(N_per_axis):
#         for j in range(N_per_axis):
#             for k in range(N_per_axis):
#                 yield ( (i*dL,(i+1)*dL) , (j*dL,(j+1)*dL) , (k*dL,(k+1)*dL) )


# def number_field(box_size,
#                  n_mesh,
#                  kmin,
#                  kmax):
#     '''
#     For a Fourier grid specified by box_size and n_mesh, return
#     the inverse FFT of a grid with ones where the k vector norm
#     is between kmin and kmax, and zeros elsewhere
#     '''
#     return k_field(box_size,n_mesh,kmin,kmax,0.)


def k_field(box_size,
            n_mesh,
            kmin,
            kmax,
            p):
    """Construct filtered grid of |k|^p values.
    
    For a specified by box_size and n_mesh, construct a Fourier grid with
    values equal to |\vec{k}|^p where the k vector norm is between kmin and kmax,
    and zeros elsewhere. The routine returns the inverse FT of this grid,
    which can then be used to compute the average |\vec{k}| values within the
    the k bin or a bispectrum triangle bin.
    
    Parameters
    ----------
    box_size : array_like
        array with shape(3) specifying the x,y,z box side lengths
    n_mesh : array_like
        array with shape(3) specifying the x,y,z mesh resolutions
    kmin : float
        minimum k value for filtering
    kmax : float
        maximum k value for filtering
    p : float
        power to raise |\vec{k}| values to on grid
    
    Returns
    -------
    RealField
        filtered real-space field as described above
    """

    # Check that we've gotten list of box size and Nmesh in x,y,z
    if len(box_size) != 3: raise ValueError('box_size must be 3-element list!')
    if len(n_mesh) != 3: raise ValueError('n_mesh must be 3-element list!')

    # Function to apply to mesh, that sets modes to unity within bin and zero elsewhere
    def nbin_mask(k,v):
        # Compute |\vec{k}| values on the mesh
        kk = sum(ki ** 2. for ki in k)**0.5
        # Return array of ones, masked according to bin boundaries
        return kk**p * ((kk<=kmax) & (kk>=kmin))

    # Make a mesh with the desired box size and resolution.
    pm = ParticleMesh(BoxSize=np.ndarray.tolist(box_size),Nmesh=np.ndarray.tolist(n_mesh))
    cf = ComplexField(pm)
    mesh = nbk.FieldMesh(cf)
    mesh = mesh.apply(nbin_mask,mode='complex',kind='wavenumber')
    painted_real_field = mesh.paint(mode='real')

    return painted_real_field


def number_field(box_size,
                n_mesh,
                kmin,
                kmax):
    """Construct filtered grid of |k|^p values.
    
    For a specified by box_size and n_mesh, construct a Fourier grid with
    values equal to 1 where the k vector norm is between kmin and kmax,
    and zeros elsewhere. The routine returns the inverse FT of this grid,
    which can then be used to compute the number of modes within the
    the k bin or a bispectrum triangle bin.
    
    Parameters
    ----------
    box_size : array_like
        array with shape(3) specifying the x,y,z box side lengths
    n_mesh : array_like
        array with shape(3) specifying the x,y,z mesh resolutions
    kmin : float
        minimum k value for filtering
    kmax : float
        maximum k value for filtering
    
    Returns
    -------
    RealField
        filtered real-space field as described above
    """

    # Check that we've gotten list of box size and Nmesh in x,y,z
    if len(box_size) != 3: raise ValueError('box_size must be 3-element list!')
    if len(n_mesh) != 3: raise ValueError('n_mesh must be 3-element list!')

    # Function to apply to mesh, that sets modes to unity within bin and zero elsewhere
    def nbin_mask(k,v):
        # Compute |\vec{k}| values on the mesh
        kk = sum(ki ** 2. for ki in k)**0.5
        # Make array of ones with same shape as input mesh
        vc = np.ones(v.shape)
        # Return array of ones, masked according to bin boundaries
        return vc * ((kk<=kmax) & (kk>=kmin))

    # Make a mesh with the desired box size and resolution.
    pm = ParticleMesh(BoxSize=np.ndarray.tolist(box_size),Nmesh=np.ndarray.tolist(n_mesh))
    cf = ComplexField(pm)
    mesh = nbk.FieldMesh(cf)
    mesh = mesh.apply(nbin_mask,mode='complex',kind='wavenumber')
    painted_real_field = mesh.paint(mode='real')

    return painted_real_field



def pk_FFT(mesh,kmin,kmax):
    """Measure binned power spectrum using an FFT-based estimator.
    
    For a k bin specified by kmin and kmax, the routine computes the number of modes
    within the bin, the mean |\vec{k}| within the bin, and the mean P(k) within the bin,
    using FFT-based methods.
    
    Parameters
    ----------
    mesh : MeshSource
        mesh to measure power spectrum of
    kmin : float
        lower k bin edge
    kmax : float
        upper k bin edge
        
    Returns
    -------
    float
        mean P(k) within bin
    float
        number of modes within bin
    float
        mean k within bin
    """
    comm = CurrentMPIComm.get()

    # Get real-space field that's the inverse FFT of a grid with ones for
    # k vectors within the desired k-bin, and zeros elsewhere
    painted_real_field = number_field(mesh.attrs['BoxSize'],mesh.attrs['Nmesh'],
                                     kmin,kmax)

    # From masked mesh, compute number of discrete modes within bin.
    # The sum below will be computed separately on each process, so we
    # need to gather all the results together
    Nbin_local = np.sum(painted_real_field**2.) / np.prod(mesh.attrs['Nmesh'])
    Nbin_sum = np.sum(comm.gather(Nbin_local))
    Nbin = comm.bcast(Nbin_sum)

    # Get real-space field that's the inverse FFT of a grid with |\vec{k}|^0.5 for
    # k vectors within the desired k-bin, and zeros elsewhere. Squaring this grid and
    # dividing by Nmesh^3 * Nbin will give the average value of |\vec{k}| over all
    # discrete modes in the bin
    painted_real_field = k_field(mesh.attrs['BoxSize'],mesh.attrs['Nmesh'],kmin,kmax,0.5)
    k_mean_local = np.sum(painted_real_field**2.) / np.prod(mesh.attrs['Nmesh']) / Nbin
    k_mean_sum = np.sum(comm.gather(k_mean_local))
    k_mean = comm.bcast(k_mean_sum)

    # Function to apply to mesh, to mask modes outside bin
    def mask(k,v):
        # Compute |\vec{k}| values on the mesh
        kk = sum(ki ** 2. for ki in k)**0.5
        # Mask out values outside desired bin
        return v * ((kk<=kmax) & (kk>=kmin))

    # Make another local copy of the mesh, with this mask applied.
    # Currently, I implement an "apply" function myself, instead of
    # using mesh.apply, since this doesn't appear to work when you have
    # a CatalogMesh with compensated=True.
#     mesh_copy = mesh.apply(mask,mode='complex',kind='wavenumber')
    mesh_copy = mesh.view()
    mesh_copy._actions.append(('complex', mask, 'wavenumber'))
    painted_real_field = mesh_copy.paint(mode='real')

    # Compute binned power spectrum, using previously-calculated value of Nbin
    P_local = np.sum(painted_real_field**2.) \
                * mesh_copy.attrs['BoxSize'].prod() \
                / mesh_copy.attrs['Nmesh'].prod() \
                / Nbin
    P_sum = np.sum(comm.gather(P_local))
    P = comm.bcast(P_sum)
    
    return P, Nbin, k_mean


def geometric_k_mean(kmin,
                     kmax):
    """Analytically approximate mean |\vec{k}| between kmin and kmax
    
    We find the geometric mean, given by the following formula:
        k_mean = [ \int_{k_min}^{k_max} d^3\vec{k} k ] / [ \int_{k_min}^{k_max} d^3\vec{k} ]
               = (3/4) [ k_max^4 - k_min^4 ] / [ k_max^3 - k_min^3 ]
               
    Parameters
    ----------
    kmin : float
        lower k bin edge
    kmax : float
        upper k bin edge
    
    Returns
    -------
    float
        k_mean
    """
    return 0.75 * (kmax**4 - kmin**4) / (kmax**3 - kmin**3)


def compute_Nbin(box_size,
                 n_mesh,
                 bins,
                 verbose,
                 return_number_fields = False):
    """Compute number of triangles in a bin in (k1,k2,k3)
    
    This uses the FFT-based algorithm.
    
    Parameters
    ----------
    box_size : array_like
        array with shape(3) specifying the x,y,z box side lengths
    n_mesh : array_like
        array with shape(3) specifying the x,y,z mesh resolutions
    bins : array_like
        array with shape (3,2) specifying kmin and kmax for the k1,k2,k3 bins
    verbose : int
        verbosity level of output: 0 for nothing, >0 for more than nothing
    return_number_fields : bool, optional
        whether to also return the number_fields used in the computation
    
    Returns
    -------
    float
        number of triangles in bin
    dict
        dict of RealField objects for k1,k2,k3 (with keys 0,1,2)
    """
    
    # Tolerance for testing whether k bins are the same
    tol = 0.01
    
    comm = CurrentMPIComm.get()
    
    # Initialize dict of up to three number_fields
    number_fields = {0 : None, 1 : None, 2 : None}

    if verbose > 0 and comm.rank == 0: print('Starting Nbin computations')
        
    # Loop over k1,k2,k3 bins that define triangle bin
    for i in range(3):
        if verbose > 0 and comm.rank == 0: print('\tConstructing masked field for bin %d' % i)

        # For k2 (i=1), check whether k1 (i=0) bin is the same. For k3 (i=2), check
        # whether k0 or k1 are the same. If any of these are true, we can re-use
        # the corresponding previously-computed number_field
        for j in range(i):
            if (np.abs(bins[i][0]-bins[j][0])/bins[i][0] < tol) \
             and (np.abs(bins[i][1]-bins[j][1])/bins[i][1] < tol):
                number_fields[i] = number_fields[j]
        
        # If the number field for this i has not been computed already, do it
        if number_fields[i] is None:
            number_fields[i] = number_field(box_size,n_mesh,bins[i][0],bins[i][1])
            
    # From masked mesh, compute number of discrete modes within bin
    if verbose > 0 and comm.rank == 0: print('\tSumming over product of masked fields')
    Nbin_local = np.sum(number_fields[0]*number_fields[1]*number_fields[2]) \
                 / np.prod(n_mesh)
    Nbin_sum = np.sum(comm.gather(Nbin_local))
    Nbin = comm.bcast(Nbin_sum)
    
    # Return Nbin, and if requested, also return number_fields (for possible use
    # by the calling routine)
    if return_number_fields:
        return Nbin,number_fields
    else:
        return Nbin


def compute_k_means_on_grid(box_size,
                            n_mesh,
                            bins,
                            verbose,
                            Nbin,
                            number_fields):
    """Compute mean k_i (vector norm) values in a bin in (k1,k2,k3)
    
    This uses the FFT-based algorithm, and takes as inputs the number of
    triangles in the bin and the dict of RealField objects used for computing
    the number of triangles.
    
    Parameters
    ----------
    box_size : array_like
        array with shape(3) specifying the x,y,z box side lengths
    n_mesh : array_like
        array with shape(3) specifying the x,y,z mesh resolutions
    bins : array_like
        array with shape (3,2) specifying kmin and kmax for the k1,k2,k3 bins
    verbose : int
        verbosity level of output: 0 for nothing, >0 for more than nothing
    Nbin : float
        number of triangles in bin
    number_fields : dict
        dict of number_fields obtained from compute_Nbin
    
    Returns
    -------
    dict
        dict of mean k values
    """
    
    # Tolerance for testing whether k bins are the same
    tol = 0.01
    
    comm = CurrentMPIComm.get()
    
    k_mean = {}
    
    if verbose > 0 and comm.rank == 0: print('Starting k_mean computations')
    
    # Loop over k1,k2,k3
    for i in range(3):
        
        # Compute k_field if i=0 or if the bin for this i is different from
        # that for i=1. (If bin is the same, no need to recompute k_field.)
        if i == 0 or (np.abs(bins[i][0]-bins[i-1][0])/bins[i][0] > tol) \
         or (np.abs(bins[i][1]-bins[i-1][1])/bins[i][1] > tol):
            kk_field = k_field(box_size,n_mesh,bins[i][0],bins[i][1],1.)
    
        # Compute mean k1, k2, or k3 in triangle bin
        if i == 0:
            k_mean_local = np.sum(kk_field*number_fields[1]*number_fields[2]) \
                            / np.prod(n_mesh) / Nbin
        elif i == 1:
            k_mean_local = np.sum(number_fields[0]*kk_field*number_fields[2]) \
                            / np.prod(n_mesh) / Nbin
        elif i == 2:
            k_mean_local = np.sum(number_fields[0]*number_fields[1]*kk_field) \
                            / np.prod(n_mesh) / Nbin

        k_mean_sum = np.sum(comm.gather(k_mean_local))
        k_mean[i] = comm.bcast(k_mean_sum)
        
    return k_mean


def compute_bk_FFT_value(mesh,
                         bins,
                         Nbin = 1,
                         verbose = 0,
                         second_mesh = None,
                         third_mesh = None):
    """Compute bispectrum value within specified bin in (k1,k2,k3)
    
    This uses the FFT-based algorithm, and takes as input the number of
    triangles in the bin. A cross bispectrum between 2 or 3 fields can
    be computed as well. If 2 fields A and B are specified, <AAB> is
    computed.
    
    Parameters
    ----------
    mesh : MeshSource
        mesh with grid to measure bispectrum of
    bins : array_like
        array with shape (3,2) specifying kmin and kmax for the k1,k2,k3 bins
    Nbin : float, optional
        number of triangles in bin (set to 1 if not specified)
    verbose : int, optional
        verbosity level of output: 0 (default) for nothing, >0 for more than nothing
    second_mesh : MeshSource, optional
        second mesh for computing cross bispectrum
    third_mesh : MeshSource, optional
        third mesh for computing cross bispectrum
    
    Returns
    -------
    float
        bispectrum value in bin
    """
    
    # Tolerance for testing whether k bins are the same
    tol = 0.01
    
    comm = CurrentMPIComm.get()
    
    delta_fields = {}

    if verbose > 0: print('Starting bispectrum computations')
        
    # Loop over k1,k2,k3 bins
    for i in range(3):
        # Function to apply to mesh, to mask modes outside bin
        def mask(k,v):
            # Compute |\vec{k}| values on the mesh
            kk = sum(ki ** 2. for ki in k)**0.5
            # Mask out values outside desired bin
            return v * ((kk<=bins[i][1]) & (kk>=bins[i][0]))
        
        if verbose > 0 and comm.rank == 0: print('\tConstructing masked field for bin %d' % i)

        # Make another local copy of the mesh.
        # If one mesh, take <AAA>. If two, take <AAB>. If three, take <ABC>.
        # (A,B,C denote up to three input meshes.)
        if i == 0:
            mesh_copy = mesh.view()
        elif i == 1:
            if third_mesh is not None:
                mesh_copy = second_mesh.view()
            else:
                mesh_copy = mesh.view()
        elif i == 2:
            if third_mesh is not None:
                mesh_copy = third_mesh.view()
            elif second_mesh is not None:
                mesh_copy = second_mesh.view()
            else:
                mesh_copy = mesh.view()
                
        # Apply mask and paint, or re-use previously computed field, if bins
        # match and the same mesh is used
        if i == 0:
            mesh_copy._actions.append(('complex', mask, 'wavenumber'))
            delta_fields[i] = mesh_copy.paint(mode='real')
        elif i == 1:
            if third_mesh is None and ((np.abs(bins[i][0]-bins[i-1][0])/bins[i][0] < tol) \
             or (np.abs(bins[i][1]-bins[i-1][1])/bins[i][1] < tol)):
                delta_fields[i] = delta_fields[i-1]
            else:
                mesh_copy._actions.append(('complex', mask, 'wavenumber'))
                delta_fields[i] = mesh_copy.paint(mode='real')
        elif i == 2:
            if second_mesh is None and ((np.abs(bins[i][0]-bins[i-1][0])/bins[i][0] < tol) \
             or (np.abs(bins[i][1]-bins[i-1][1])/bins[i][1] < tol)):
                delta_fields[i] = delta_fields[i-1]
            else:
                mesh_copy._actions.append(('complex', mask, 'wavenumber'))
                delta_fields[i] = mesh_copy.paint(mode='real')
        
    # Compute binned bispectrum, using value of Nbin passed to routine.
    # Return B
    if verbose > 0 and comm.rank == 0: print('\tSumming over product of masked fields')
    B_local = np.sum(delta_fields[0]*delta_fields[1]*delta_fields[2]) \
                * mesh_copy.attrs['BoxSize'].prod()**2 \
                / np.prod(mesh_copy.attrs['Nmesh']) \
                / Nbin
    B_sum = np.sum(comm.gather(B_local))
    B = comm.bcast(B_sum)
    
    return B


def bk_FFT_full(mesh,
           bin0,bin1,bin2,
           verbose = 0,
           second_mesh = None,
           third_mesh = None,
           approximate_k_means = False):
    """Compute bispectrum value, number of triangles, and mean k_i values within specified bin in (k1,k2,k3)
    
    This uses the FFT-based algorithm. A cross bispectrum between 2 or 3 fields can
    be computed as well. If 2 fields A and B are specified, <AAB> is computed.
    
    This is essentially a wrapper for compute_Nbin, compute_k_means_on_grid, and
    compute_bk_FFT_value.
    
    Parameters
    ----------
    mesh : MeshSource
        mesh with grid to measure bispectrum of
    bin0, bin1, bin2 : array_like
        arrays with shape(2) specifying kmin and kmax for the k1,k2,k3 bins
    verbose : int, optional
        verbosity level of output: 0 (default) for nothing, >0 for more than nothing
    second_mesh : MeshSource, optional
        second mesh for computing cross bispectrum
    third_mesh : MeshSource, optional
        third mesh for computing cross bispectrum
    approximate_k_means : bool, optional
        if True, return approximate k_mean values instead of explicit grid computations
        
    Returns
    -------
    float
        bispectrum value in bin
    float
        number of triangles in bin
    float, float, float
        mean |\vec{k}_i| values in triangle bin
    """

    comm = CurrentMPIComm.get()
    
    if third_mesh != None and second_mesh == None:
        raise ValueError('Must specify second_mesh if specifying third_mesh!')
    
    if second_mesh != None and not np.array_equal(mesh.attrs['BoxSize'],second_mesh.attrs['BoxSize']):
        raise ValueError('BoxSize mismatch between first and second mesh')
    if third_mesh != None and not np.array_equal(mesh.attrs['BoxSize'],third_mesh.attrs['BoxSize']):
        raise ValueError('BoxSize mismatch between first and third mesh')
    if second_mesh != None and not np.array_equal(mesh.attrs['Nmesh'],second_mesh.attrs['Nmesh']):
        raise ValueError('Nmesh mismatch between first and second mesh')
    if third_mesh != None and not np.array_equal(mesh.attrs['Nmesh'],third_mesh.attrs['Nmesh']):
        raise ValueError('Nmesh mismatch between first and third mesh')
    
    bins = {0:bin0, 1:bin1, 2:bin2}
    
    # Get number of triangles in bin, and number fields (for re-use in compute_k_means_on_grid)
    Nbin,number_fields = compute_Nbin(mesh.attrs['BoxSize'],mesh.attrs['Nmesh'],bins,verbose,
                                      return_number_fields=True)

    # Getting k_mean for each k argument
    if not approximate_k_means:
        # Compute explicitly
        k_mean = compute_k_means_on_grid(mesh.attrs['BoxSize'],mesh.attrs['Nmesh'],
                                         bins,verbose,Nbin,number_fields)
    else:
        # Or, if desired, approximate k_mean analytically (manually de-assigning
        # number_fields dict to free up memory)
        number_fields = {}
        k_mean = {}
        for i in range(3):
            k_mean[i] = geometric_k_mean(bins[i][0],bins[i][1])

    # Compute mean B over bin
    B = compute_bk_FFT_value(mesh,bins,Nbin,verbose,second_mesh,third_mesh)

    return B,Nbin,k_mean[0],k_mean[1],k_mean[2]


def bk_FFT_grid_info(mesh,
                     bin0,bin1,bin2,
                     verbose = 0):
    """Compute number of triangles and mean k_i values within specified bin in (k1,k2,k3)
    
    This uses the FFT-based algorithm
    
    This is essentially a wrapper for compute_Nbin and compute_k_means_on_grid.
    
    Parameters
    ----------
    mesh : MeshSource
        mesh with grid to compute information for
    bin0, bin1, bin2 : array_like
        arrays with shape(2) specifying kmin and kmax for the k1,k2,k3 bins
    verbose : int, optional
        verbosity level of output: 0 (default) for nothing, >0 for more than nothing
        
    Returns
    -------
    float
        number of triangles in bin
    float, float, float
        mean |\vec{k}_i| values in triangle bin
    """

    comm = CurrentMPIComm.get()
    
    bins = {0:bin0, 1:bin1, 2:bin2}
    
    # Get number of triangles in bin
    Nbin,number_fields = compute_Nbin(mesh.attrs['BoxSize'],mesh.attrs['Nmesh'],bins,verbose,
                                      return_number_fields=True)

    # Get k_mean for each k argument
    k_mean = compute_k_means_on_grid(mesh.attrs['BoxSize'],mesh.attrs['Nmesh'],
                                     bins,verbose,Nbin,number_fields)

    return Nbin,k_mean[0],k_mean[1],k_mean[2]


def bk_FFT_unnormalized_value(mesh,
                              bin0,bin1,bin2,
                              verbose = 0,
                              second_mesh = None,
                              third_mesh = None):
    """Compute unnormalized bispectrum value within specified bin in (k1,k2,k3).
    
    Unnormalized means that the value is not divided by the number of triangles in
    the bin - effectively, it is a sum of \delta(k_1)\delta(k_2)\delta(k_3)
    for all triplets of (k_1,k_2,k_3) vectors in the bin.
    
    Parameters
    ----------
    mesh : MeshSource
        mesh with grid to compute information for
    bin0, bin1, bin2 : array_like
        arrays with shape(2) specifying kmin and kmax for the k1,k2,k3 bins
    verbose : int, optional
        verbosity level of output: 0 (default) for nothing, >0 for more than nothing
    second_mesh : MeshSource, optional
        second mesh for computing cross bispectrum
    third_mesh : MeshSource, optional
        third mesh for computing cross bispectrum
        
    Returns
    -------
    float
        unnormalized bispectrum value in bin
    """

    comm = CurrentMPIComm.get()
    
    bins = {0:bin0, 1:bin1, 2:bin2}
    
    B = compute_bk_FFT_value(mesh,bins,Nbin=1,verbose=verbose,
                             second_mesh=second_mesh,third_mesh=third_mesh)

    return B


def bk_binned(mesh,
              bin0,bin1,bin2,
              verbose = 0):
    """Compute bispectrum value, number of triangles, and mean k_i values within specified bin in (k1,k2,k3) using an inefficient 'counting triangles' estimator.

    NOTE: This hasn't been tested in MPI, so it's very likely that it will give incorrect
    results when run with more than one process! For checking things within an iPython
    notebook (running in a single process), though, it gives reliable results.
    
    Parameters
    ----------
    mesh : MeshSource
        mesh with grid to measure bispectrum of
    bin0, bin1, bin2 : array_like
        arrays with shape(2) specifying kmin and kmax for the k1,k2,k3 bins
    verbose : int, optional
        verbosity level of output: 0 (default) for nothing, >0 for more than nothing
        
    Returns
    -------
    float
        bispectrum value in bin
    float
        number of triangles in bin
    float, float, float
        mean |\vec{k}_i| values in triangle bin
    """
    from nbodykit.meshtools import SlabIterator

    # Define lists of bin and dimension indices
    bin_i_list = [0,1,2]
    dims = [0,1,2]

    # Define list of bin boundaries and squared bin boundaries
    bins = {0:bin0, 1:bin1, 2:bin2}
    binsq = {}
    for i in bin_i_list:
        binsq[i] = [x**2 for x in bins[i]]

    # From input mesh, get complex field and k coordinates (called x here)
    y3d = mesh.paint(mode='complex')
    x3d = y3d.x

    # From input mesh, get k_f = 2\pi/L_box
    kf = 2*np.pi/mesh.attrs['BoxSize'].mean()

    # Assume input mesh is Hermitian symmetric, and set symmetry axis appropriately
    symmetry_axis = -1

    # Define empty dicts/arrays to hold k values, k vector components, and
    # \delta values for points falling within each k-bin. The first index of each
    # (which we use b for) denotes the bin index, and the second index of ki_bin_arr
    # denotes spatial dimension (x,y,z)
    k_bin_arr = {0:np.empty(0), 1:np.empty(0), 2:np.empty(0)}
    ki_bin_arr = {}
    for b in bin_i_list:
        ki_bin_arr[b] = {0:np.empty(0), 1:np.empty(0), 2:np.empty(0)}
    y3d_bin_arr = {0:np.empty(0), 1:np.empty(0), 2:np.empty(0)}

    # Loop over y-z planes of the coordinate mesh, grabbing points within each k-bin
    # and adding them to the dicts/array defined above
    if verbose > 0: print('Starting loop over slabs')
    for slab in SlabIterator(x3d, axis=0, symmetry_axis=symmetry_axis):

        # Get \delta(k) values on slab
        y3d_copy = y3d.copy()
        y3d_slab = y3d_copy[slab.index]

        # Get k^2 values on slab
        xslab = slab.norm2()
        if len(xslab.flat) == 0: continue

        # Get kx,ky,kz values on slab, as arrays with same dimensions as
        # xslab and y3d_slab
        ki = {}
        ki[0] = slab.coords(0)*np.ones_like(slab.coords(1))*np.ones_like(slab.coords(2))
        ki[1] = np.ones_like(slab.coords(0))*slab.coords(1)*np.ones_like(slab.coords(2))
        ki[2] = np.ones_like(slab.coords(0))*np.ones_like(slab.coords(1))*slab.coords(2)

        # Loop over bins
        for b in bin_i_list:

            # Get a (flattened) array of bin indices that each point in slab falls into.
            # Since we're only testing one bin at a time, 0=below bin, 1=in bin, 2=above bin
            dig_k = np.digitize(xslab.flat, binsq[b])

            # Make slab mask corresponding to points in bin
            bin_mask = dig_k==1

            # Collect k,kx,ky,kz,y3d values for points within bin, as well as
            # Hermitian weights for these points
            k_here = xslab.flat[bin_mask]**0.5
            ki_here = {}
            for i in dims:
                ki_here[i] = ki[i].flat[bin_mask]
            herm_weights_here = slab.hermitian_weights.flat[bin_mask]
            y3d_here = y3d_slab.flat[bin_mask]

            # If a point has Hermitian weight 2, it means
            # that the value of \delta(-\vec{k}) is not explicitly stored on the grid,
            # so we append it to our lists, using the fact that
            # \delta(-\vec{k}) = \delta(\vec{k})*
            for i,w in enumerate(herm_weights_here):
                if np.abs(w-2)<1e-4:
                    k_here = np.append(k_here,k_here[i])
                    for j in dims:
                        ki_here[j] = np.append(ki_here[j],-ki_here[j][i])
                    y3d_here = np.append(y3d_here,np.conj(y3d_here[i]))

            # Append points we found in this slab to the master lists for the given k-bin
            k_bin_arr[b] = np.append(k_bin_arr[b],k_here)
            for j in dims:
                ki_bin_arr[b][j] = np.append(ki_bin_arr[b][j],ki_here[j])
            y3d_bin_arr[b] = np.append(y3d_bin_arr[b],y3d_here)

            # Print the lists of stuff we obtained for this slab, for debugging
            if verbose > 2:
                print('new iteration, bin %d' % b)
                print(herm_weights_here)
                print(k_here)
                print(ki_here[0])
                print(ki_here[1])
                print(ki_here[2])
                print(y3d_here)
                print('')

    # Initialize bispectrum variable, k_mean variables, and bin counter
    B = 0.
    Nbin = 0
    k_mean = np.array([0.,0.,0.])

    # Step through combinations of points in the 3 bins using a nested loop.
    # For each combination, check whether the k vectors form a closed triangle.
    # If so, increase the bin counter and add the product of \delta values to B
    if verbose > 0: print('Starting loop over points in bins')
    for i0,k0 in enumerate(k_bin_arr[0]):
        for i1,k1 in enumerate(k_bin_arr[1]):
            for i2,k2 in enumerate(k_bin_arr[2]):

                # Square of k_1+k_2+k_3
                k_tot_sq = np.sum(np.fromiter(((ki_bin_arr[0][j][i0] \
                                   + ki_bin_arr[1][j][i1] \
                                   + ki_bin_arr[2][j][i2])**2. for j in dims), float))

                # Test if it's zero, and if so, add to B, k_mean array, and Nbin
                if np.abs(k_tot_sq)<1e-4*kf:
                    if verbose > 1:
                        print('Found a triangle! i0,i1,i2 = %d %d %d' % (i0,i1,i2))
                        print([ki_bin_arr[0][j][i0] for j in dims])
                        print([ki_bin_arr[1][j][i1] for j in dims])
                        print([ki_bin_arr[2][j][i2] for j in dims])
                        print(y3d_bin_arr[0][i0]*y3d_bin_arr[1][i1]*y3d_bin_arr[2][i2])

                    Nbin += 1
                    for b,i in zip(bin_i_list,[i0,i1,i2]):
                        k_mean[b] += k_bin_arr[b][i]

                    B += y3d_bin_arr[0][i0]*y3d_bin_arr[1][i1]*y3d_bin_arr[2][i2]

    # Normalize B by L^6
    B *= mesh.attrs['BoxSize'].prod()**2.

    # If Nbin=0, don't divide by it!
    if Nbin > 0:
        B /= Nbin
        for j in bin_i_list:
            k_mean[j] /= Nbin

    # Return results
    if verbose > 0: print('B, Nbin = %g %g' % (B,Nbin))
    return B,Nbin,k_mean[0],k_mean[1],k_mean[2]



def generate_bin_edge_list(kmin = -1.,
                           kmax = -1.,
                           dk = -1.,
                           num_lowk_bins = 0,
                           dk_high = -1.):
    """Generate list of k bin edges, based on input kmin, kmax,
    dk.
    
    If desired, low-k and high-k bins can have different
    widths - this is specified through num_lowk_bins and dk_high.
    For example, (kmin=0.05, kmax=1., dk=0.05, num_lowk_bins=10, dk_high=0.1)
    will result in the 10 lowest-k bins having width 0.05, and the remaining
    higher-k bins having width 0.1.
    
    Parameters
    ----------
    kmin : float
        lower edge of lowest bin
    kmax : float
        upper edge of highest bin
    dk : float
        bin width
    num_lowk_bins : float, optional
        number of low-k bins with width `dk`
    dk_high : float, optional
        width for higher-k bins
    
    Returns
    -------
    array_like
        shape(Nbins,2) array of lower and upper bin edges
    """
    
    # Check that kmin, kmax, and dk are greater than zero
    if kmin <= 0.: raise ValueError('kmin must be > 0!')
    if kmax <= 0.: raise ValueError('kmax must be > 0!')
    if dk <= 0.: raise ValueError('dk must be > 0!')
    
    # Check for different dk at high k
    if num_lowk_bins > 0 and dk_high < 0.:
        raise ValueError('Must specify dk_high if num_lowk_bins > 0!')
    
    # Make evenly-spaced arrays of k_lower, k_upper, and dk values,
    # with k_lower from kmin to kmax-dk
    kL_arr = np.arange(kmin,kmax-dk,dk)
    dk_arr = np.ones_like(kL_arr)*dk
    kU_arr = kL_arr+dk_arr
    
    # If we want to use a different dk at higher k, make separate
    # lists of the higher-k bin edges, and then appending them to the
    # low-k bin edge lists
    if num_lowk_bins > 0 and num_lowk_bins < len(kL_arr):
        kLhigh_arr = np.arange(kU_arr[num_lowk_bins-1],kmax-dk_high,dk_high)
        dkhigh_arr = np.ones_like(kLhigh_arr)*dk_high
        kUhigh_arr = kLhigh_arr+dkhigh_arr
        
        kL_arr = np.hstack((kL_arr[:num_lowk_bins],kLhigh_arr))
        dk_arr = np.hstack((dk_arr[:num_lowk_bins],dkhigh_arr))
        kU_arr = np.hstack((kU_arr[:num_lowk_bins],kUhigh_arr))
    
    # Join the upper and lower bin edge lists into one list of (kL,kU)
    kLU_arr = np.vstack((kL_arr[:],kU_arr[:])).T
    
    return kLU_arr


def generate_equilateral_triangle_bin_list(kmin = -1.,
                                           kmax = -1.,
                                           dk = -1.,
                                           num_lowk_bins = 0,
                                           dk_high = -1.,
                                           return_indices = False):
    """Generate list of (k1,k2,k3) bin boundaries for equilateral bispectrum triangles.
    
    The output array has shape(N_tri,6) where N_tri is the number of
    triangle bins, and the second index is over
        [k1_min,k1_max,k2_min,k2_max,k3_min,k3_max] .
        
    Parameters
    ----------
    kmin : float
        lower edge of lowest bin
    kmax : float
        upper edge of highest bin
    dk : float
        bin width
    num_lowk_bins : float, optional
        number of low-k bins with width `dk`
    dk_high : float, optional
        width for higher-k bins
    return_indices : bool, optional
        instead of bin edges, return shape(N_tri,3) array of k-bin indices
    
    Returns
    -------
    array_like
        either shape(N_tri,6) array of bin edges or shape(N_tri,3) array of k-bin indices
    """
    
    # Get list of k bin edges
    kLU_arr = generate_bin_edge_list(kmin,kmax,dk,num_lowk_bins,dk_high)
    
    if return_indices:
        # Generate list of k1,k2,k3 bin indices, which are all the same
        index_arr = range(len(kLU_arr))
        tri_bin_indices = np.vstack((index_arr,index_arr,index_arr)).T
        return tri_bin_indices
    else:
        # Generate list of k1,k2,k3 bin edges, which are all the same
        tri_bin_edges = np.hstack((kLU_arr,kLU_arr,kLU_arr))
        return tri_bin_edges


def generate_squeezed_triangle_bin_list(kmin = -1.,
                                        kmax = -1.,
                                        dk = -1.,
                                        squeezed_bin_index = 0,
                                        num_lowk_bins = 0,
                                        dk_high = -1.,
                                        return_indices = False):
    """Generate list of (k1,k2,k3) bin boundaries for squeezed bispectrum triangles.
    
    The output array has shape(N_tri,6) where N_tri is the number of
    triangle bins, and the second index is over
        [k1_min,k1_max,k2_min,k2_max,k3_min,k3_max] .
        
    The squeezed triangles are defined isosceles, with the two 'short' k bins identical
    and the 'long' bin specified by its k bin index through squeezed_bin_index.
    
    Parameters
    ----------
    kmin : float
        lower edge of lowest bin
    kmax : float
        upper edge of highest bin
    dk : float
        bin width
    squeezed_bin_index : int
        k bin index of k_long bin
    num_lowk_bins : float, optional
        number of low-k bins with width `dk`
    dk_high : float, optional
        width for higher-k bins
    return_indices : bool, optional
        instead of bin edges, return shape(N_tri,3) array of k-bin indices
    
    Returns
    -------
    array_like
        either shape(N_tri,6) array of bin edges or shape(N_tri,3) array of k-bin indices
    """
    
    # Get list of k bin edges
    kLU_arr = generate_bin_edge_list(kmin,kmax,dk,num_lowk_bins,dk_high)
    
    if return_indices:
        # Return array of triplets of k bin indices
        short_index_arr = range(squeezed_bin_index+1,len(kLU_arr))
        squeezed_index_arr = [squeezed_bin_index for i in range(len(short_index_arr))]
        tri_bin_indices = np.vstack((squeezed_index_arr,short_index_arr,short_index_arr)).T
        return tri_bin_indices
    else:
        # Construct arrays of bin edges for squeezed side and non-squeezed sides
        squeezed_bin = kLU_arr[squeezed_bin_index]
        kLU_arr = kLU_arr[squeezed_bin_index+1:]
        squeezed_bin_arr = np.array([squeezed_bin for i in range(len(kLU_arr))])

        # Generate list of k1,k2,k3 bin edges
        tri_bin_edges = np.hstack((squeezed_bin_arr,kLU_arr,kLU_arr))

        return tri_bin_edges


def generate_isosceles_triangle_bin_list(kmin = -1.,
                                        kmax = -1.,
                                        dk = -1.,
                                        isos_mult = 0,
                                        isos_tol = 0.1,
                                        num_lowk_bins = 0,
                                        dk_high = -1.,
                                        return_indices = False):
    """Generate list of (k1,k2,k3) bin boundaries for isosceles bispectrum triangles.
    
    The output array has shape(N_tri,6) where N_tri is the number of
    triangle bins, and the second index is over
        [k1_min,k1_max,k2_min,k2_max,k3_min,k3_max] .
        
    Isosceles triangles are defined by the multiplier isos_mult, with k2 and k3
    equal to isos_mult * k1. The isosceles condition
    is not enforced for each triangle with a bin; rather, for each k_short bin, we
    find the corresponding k_long bin such that the mean k within the k_long bin
    is as close as possible to satisfying the condition. If there is no bin for
    which this is true to within isos_tol tolerance, that k_short bin is skipped.
    
    Parameters
    ----------
    kmin : float
        lower edge of lowest bin
    kmax : float
        upper edge of highest bin
    dk : float
        bin width
    isos_mult : float
        multiplier defining triangles
    isos_tol : float
        fractioanl tolerance for the isosceles condition to be satisfied
    num_lowk_bins : float, optional
        number of low-k bins with width `dk`
    dk_high : float, optional
        width for higher-k bins
    return_indices : bool, optional
        instead of bin edges, return shape(N_tri,3) array of k-bin indices
    
    Returns
    -------
    array_like
        either shape(N_tri,6) array of bin edges or shape(N_tri,3) array of k-bin indices
    """
    
    # Check that isos_mult is sensible
    if isos_mult < 1.:
        raise ValueError('isos_mult must be greater than 1! Come on man...')
    
    # Get list of k bin edges
    kLU_arr = generate_bin_edge_list(kmin,kmax,dk,num_lowk_bins,dk_high)
    
    # Compute approximate k_mean for each bin
    kmean_arr = np.array([geometric_k_mean(x[0],x[1]) for x in kLU_arr]) 
    
    # Compute desired k_mean of each corresponding kL bin
    target_klong_arr = kmean_arr/isos_mult

    # Routine that finds the index of the closest array element to value, and
    # returns that index along with the fractional deviation of that element
    # from value
    def find_nearest(array,value):
        arr = np.asarray(array)
        dev = (np.abs(array - value)).min()/value
        closest_ind = (np.abs(array - value)).argmin()
        return closest_ind,dev
    
    # Make array of indices of best-match kL bin for each kS bin, along with
    # fractional deviation of desired and actual kL_mean
    nearest_index_arr = np.array([find_nearest(kmean_arr,x) for x in target_klong_arr])

    # Mask out all kL proposals whose fractional deviation exceeds isos_tol
    mask = [x[1]<isos_tol for x in nearest_index_arr]
    
    if return_indices:
        kS_index_arr = np.arange(len(kLU_arr))[mask]
        kL_index_arr = np.arange(len(kLU_arr))[nearest_index_arr[mask][:,0].astype(int)]
        tri_bin_indices = np.vstack((kL_index_arr,kS_index_arr,kS_index_arr)).T
        return tri_bin_indices
    else:
        # Make arrays of kS and kL bins that are accepted according to our tolerance
        kS_arr = kLU_arr[mask]
        kL_arr = kLU_arr[nearest_index_arr[mask][:,0].astype(int)]

        # Make array of (kL,kS,kS) bin edges
        tri_bin_edges = np.hstack((kL_arr,kS_arr,kS_arr))

        return tri_bin_edges


def generate_triangle_bin_list(kmin = -1.,
                               kmax = -1.,
                               dk = -1.,
                               mu_min = None,
                               mu_max = None,
                               dmu = None,
                               num_fields = 1,
                               num_lowk_bins = 0,
                               dk_high = -1.,
                               return_indices = False):
    """Generate list of (k1,k2,k3) bin boundaries for general bispectrum triangles.
    
    The output array has shape(N_tri,6) where N_tri is the number of
    triangle bins, and the second index is over
        [k1_min,k1_max,k2_min,k2_max,k3_min,k3_max] ,
    subject to restrictions
        k1_min >= k2_min >= k3_min
    and
        k2_max + k3_max >= k1_min .
    If num_fields > 1, these conditions are altered to ensure that all non-redundant
    triangle bins are included in the list.
        
    In the future, binning can be specified as bins in mu (cosine of the angle between
    the two equal-length sides), but this feature has not been implemented yet.
    
    Parameters
    ----------
    kmin : float
        lower edge of lowest bin
    kmax : float
        upper edge of highest bin
    dk : float
        bin width
    mu_min : float, optional
        lower edge of lowest mu bin (NOT IMPLEMENTED)
    mu_max : float, optional
        upper edge of highest mu bin (NOT IMPLEMENTED) 
    num_fields : int, optional
        number of fields being correlated in bispectrum
    num_lowk_bins : float, optional
        number of low-k bins with width `dk`
    dk_high : float, optional
        width for higher-k bins
    return_indices : bool, optional
        instead of bin edges, return shape(N_tri,3) array of k-bin indices
    
    Returns
    -------
    array_like
        either shape(N_tri,6) array of bin edges or shape(N_tri,3) array of k-bin indices
    """
    
    # Check that num_fields is valid
    if num_fields not in [1,2,3]: raise ValueError('num_fields must be 1, 2, or 3!')
    
    # Check for mu information
    if mu_min is not None or mu_max is not None or dmu is not None:
        raise NotImplementedError('Mu binning not implemented yet!')
        
    # Get list of k bin edges
    kLU_arr = generate_bin_edge_list(kmin,kmax,dk,num_lowk_bins,dk_high)
    
    def filter_triangles(arr):
        # Remove bin combinations for which kA_max+kB_max<kC_min,
        # for (A,B,C) = any permutation of (1,2,3).
        # The input array is of triplets of (k_lower,k_upper)
        return np.array(list(filter(lambda x: x[0][1]+x[1][1]>=x[2][0]
                                                         and x[0][1]+x[2][1]>=x[1][0]
                                                         and x[1][1]+x[2][1]>=x[0][0], arr)))
    
    if num_fields == 1:
        # Make array of all combinations of 3 of k values, with replacement
        tri_arr = np.array(list(itertools.combinations_with_replacement(kLU_arr,3)))
        # Remove invalid triangles
        tri_arr_filt = filter_triangles(tri_arr)
        # Currently, the columns are such that 0<=1<=2. The output will be more intuitive if we
        # swap them around
        temp = tri_arr_filt[:,2].copy()
        tri_arr_filt[:,2] = tri_arr_filt[:,0].copy()
        tri_arr_filt[:,0] = temp
        # Reshape array into list of 6-element items
        tri_arr_filt = tri_arr_filt.reshape(tri_arr_filt.shape[0],6)
        
    elif num_fields == 2:
        # Make array of all combinations of 2 k values, with replacement
        tri_arr_2 = np.array(list(itertools.combinations_with_replacement(kLU_arr,2)))
        tri_arr = [[a[0],a[1],b] for a in tri_arr_2 for b in kLU_arr]
        # Remove invalid triangles
        tri_arr_filt = filter_triangles(tri_arr)
        # Currently, the columns are such that 0<=1. The output will be more intuitive if we
        # swap them around
        temp = tri_arr_filt[:,1].copy()
        tri_arr_filt[:,1] = tri_arr_filt[:,0].copy()
        tri_arr_filt[:,0] = temp
        # Reshape array into list of 6-element items
        tri_arr_filt = tri_arr_filt.reshape(tri_arr_filt.shape[0],6)
        
    elif num_fields == 3:
        # Make array of all permutations of 3 k values with repetition
        tri_arr = list(itertools.product(kLU_arr,repeat=3))
        # Remove invalid triangles
        tri_arr_filt = filter_triangles(tri_arr)
        # Reshape array into list of 6-element items
        tri_arr_filt = tri_arr_filt.reshape(tri_arr_filt.shape[0],6)
        
    # Finally, we sort the list by the first k_lower, then second, then third, so that
    # the first column increases most quickly. This ensures that all triangles with a lower
    # maximum side length appear in the list before the maximum side length is incremented.
    tri_bin_edges = np.sort(tri_arr_filt.view('f8,f8,f8,f8,f8,f8'), 
                           order=['f0','f2','f4'], axis=0).view(np.float)
    
    if not return_indices:
        # Return array of triangle bin edges
        return tri_bin_edges
    else:
        # Return array of k bin indices corresponding to k bins for each side
        return np.array([ [ np.where(kLU_arr[:,0]==tri_bin_edges[t][i])[0][0] 
           for i in [0,2,4] ] for t in range(len(tri_bin_edges)) ])
        

class FFTBispectrum(nbkfftpower.FFTBase):
    """A class to organize bispectrum measurements
    
    This class contains methods and containers for bispectrum measurements from
    a single grid, or in the case of cross bispectra, and given set of 2 or 3
    grids. Its structure is based on nbodykit's FFTPower.
    
    Upon initialization, the class converts and paints the input mesh(es)
    in Fourier space, and pre-computes the list of triangle grids over which
    the measurements will occur. Measurements are made with measure_bispectrum
    (slower algorithm) or measure_bispectrum_faster (measure unnormalized B values
    only, using faster algorithm) and measure_gridinfo_faster (measure triangle
    counts and k_mean values within triangle bins). Final results can be saved
    to file using save_bispectrum, or can also be saved to file at intervals
    during execution (recommended).
    
    The results from the three measurement routines are stored in a dict b,
    with keys 'index' (triangle bin indices), 'k_edge' (k1,k2,k3 bin 
    boundaries for each triangle bin), 'k_mean' (k_mean
    values in bins), 'B' (bispectrum values in bins), and 'N_tri' (triangle counts
    in bins).
    
    Parameters
    ----------
    source : CatalogSource, MeshSource
        the source for the first field; if a CatalogSource is provided, it
        is automatically converted to MeshSource using the default painting
        parameters (via :func:`~nbodykit.base.catalogmesh.CatalogMesh.to_mesh`)
    Nmesh : float, optional
        the number of cells per side in the particle mesh used to paint the source
    BoxSize : float, optional
        the size of the box
    dk : float, optional
        the linear spacing of k bins to use; if not provided, the
        fundamental mode of the box is used
    kmin : float, optional
        the lower edge of the first k bin to use. only optional if k_edges is specified
    kmax : float, optional
        the upper edge of the last k bin to use. only optional if k_edges is specified
    num_lowk_bins : float, optional
        number of low-k bins with width `dk`
    dk_high : float, optional
        width for higher-k bins
    dmu : float, optional
        NOT IMPLEMENTED
    mu_min : float, optional
        NOT IMPLEMENTED
    mu_max : float, optional
        NOT IMPLEMENTED
    pos_units_mpcoverh : float, optional
        conversion factor between position units of input grid and Mpc/h
    k_edges : float, optional
        array with shape (N_tri,6) of triangle bin edges, where second index is over
            [ k1_low k1_high k2_low k2_high k3_low k3_high ]
    second : CatalogSource, MeshSource
        second mesh to cross-correlate with
    third : CatalogSource, MeshSource
        third mesh to cross-correlate with
    triangle_type : {'all', 'equilateral', 'isosceles', 'squeezed'}
        type of triangle to measure bispectrum on
    isos_mult : float, optional
        multiplier defining triangles
    isos_tol : float, optional
        fractioanl tolerance for the isosceles condition to be satisfied
    squeezed_bin_index : int, optional
        k bin index of k_long bin
    for_grid_info_only : bool, optional
        if true, don't paint input mesh(es) - saves time if only grid info is wanted
    
    Methods
    -------
    set_k_edges(k_edges)
        set k bin edge array to an input array
    measure_bispectrum(imin=None, imax=None, kmeas_min=None, kmeas_max=None,
                       out_file=None, verbose=0)
        measure bispectrum (along with triangle counts and mean k values) 
        on all triangles bins with indices between imin and imax,
        intermittently writing results to out_file
    measure_bispectrum_faster(imin=None, imax=None, kmeas_min=None, kmeas_max=None,
                       out_file=None, verbose=0)
        measure unnormalized bispectrum  
        on all triangle bins with indices between imin and imax,
        intermittently writing results to out_file
    measure_gridinfo_faster(imin=None, imax=None, kmeas_min=None, kmeas_max=None,
                       out_file=None, verbose=0)
        measure  triangle counts and mean k values) 
        on all triangle bins with indices between imin and imax,
        intermittently writing results to out_file
    save_bispectrum(out_file)
        write bispectrum results to out_file
    """
    logger = logging.getLogger('FFTBispectrum')

    def __init__(self,
                 source,
                 Nmesh = None,
                 BoxSize = None, 
                 dk = None,
                 kmin = None,
                 kmax = None,
                 num_lowk_bins = 0,
                 dk_high = -1.,
                 dmu = None,
                 mu_min = None,
                 mu_max = None,
                 pos_units_mpcoverh = 1.,
                 k_edges = None,
                 second = None,
                 third = None,
                 triangle_type='all',
                 isos_mult = 0,
                 isos_tol = 0.1,
                 squeezed_bin_index = 0,
                 for_grid_info_only = False
                 ):
        
        # Via FFTBase initializer, cast source to MeshSource and set some attributes
        nbkfftpower.FFTBase.__init__(self, source, None, Nmesh, BoxSize)
        self.mesh = self.first
        self.second = None
        self.third = None
        
        # If second or third source arguments exist, cast them to meshes.
        # Also determine number of fields here
        self.num_fields = 1
        if second is not None:
            self.second = nbkfftpower._cast_source(second,Nmesh=Nmesh,BoxSize=BoxSize)
            self.num_fields += 1
        if third is not None:
            if second is None:
                raise ValueError('Need second source defined if third source is defined!')
            self.third = nbkfftpower._cast_source(third,Nmesh=Nmesh,BoxSize=BoxSize)
            self.num_fields += 1

        # Set dk to k_f if not specified
        if dk is None:
            dk = 2 * np.pi / self.attrs['BoxSize'].min()
           
        # Check for input kmin/kmax or k_edges
        if (kmin is None or kmax is None) and k_edges is None:
            raise ValueError('Must specify either {kmin,kmax} values or k_edges array!')

        # Save metadata
        self.attrs['dk'] = dk
        self.attrs['kmin'] = kmin
        self.attrs['kmax'] = kmax
        self.attrs['dmu'] = dmu
        self.attrs['mu_min'] = mu_min
        self.attrs['mu_max'] = mu_max
        self.attrs['pos_units_mpcoverh'] = pos_units_mpcoverh
        self.attrs['triangle_type'] = triangle_type
        self.attrs['num_lowk_bins'] = num_lowk_bins
        self.attrs['dk_high'] = dk_high
        self.attrs['isos_mult'] = isos_mult
        self.attrs['isos_tol'] = isos_tol
        self.attrs['squeezed_bin_index'] = squeezed_bin_index
        
        # Generate list of triangle bin edges, or take directly from input array
        if k_edges is not None:
            self.k_edges = k_edges
        else:
            if triangle_type == 'all':
                self.k_edges = generate_triangle_bin_list(kmin=kmin,kmax=kmax,dk=dk,
                                                          dmu=dmu,mu_min=mu_min,mu_max=mu_max,
                                                          num_fields=self.num_fields,
                                                   num_lowk_bins=num_lowk_bins,dk_high=dk_high)
                self.k_indices = generate_triangle_bin_list(kmin=kmin,kmax=kmax,dk=dk,
                                                          dmu=dmu,mu_min=mu_min,mu_max=mu_max,
                                                          num_fields=self.num_fields,
                                                   num_lowk_bins=num_lowk_bins,dk_high=dk_high,
                                                   return_indices=True)
            elif triangle_type == 'equilateral':
                self.k_edges = generate_equilateral_triangle_bin_list(kmin=kmin,kmax=kmax,dk=dk,
                                                        num_lowk_bins=num_lowk_bins,dk_high=dk_high)
                self.k_indices = generate_equilateral_triangle_bin_list(kmin=kmin,kmax=kmax,dk=dk,
                                                        num_lowk_bins=num_lowk_bins,dk_high=dk_high,
                                                        return_indices=True) 
            elif triangle_type == 'squeezed':
                self.k_edges = generate_squeezed_triangle_bin_list(kmin,kmax,dk,
                                                             squeezed_bin_index=squeezed_bin_index,
                                                        num_lowk_bins=num_lowk_bins,dk_high=dk_high)
                self.k_indices = generate_squeezed_triangle_bin_list(kmin,kmax,dk,
                                                             squeezed_bin_index=squeezed_bin_index,
                                                        num_lowk_bins=num_lowk_bins,dk_high=dk_high,
                                                        return_indices=True)
            elif triangle_type == 'isosceles':
                self.k_edges = generate_isosceles_triangle_bin_list(kmin,kmax,dk,
                                                        isos_mult=isos_mult,isos_tol=isos_tol,
                                                        num_lowk_bins=num_lowk_bins,dk_high=dk_high)
                self.k_indices = generate_isosceles_triangle_bin_list(kmin,kmax,dk,
                                                        isos_mult=isos_mult,isos_tol=isos_tol,
                                                        num_lowk_bins=num_lowk_bins,dk_high=dk_high,
                                                        return_indices=True)
                
        # Initialize variable for dicts of measured quantities
        self.b = None
        
        # Paint mesh(es) explicitly, to prevent explicit re-painting with each
        # bispectrum measurement
        if not for_grid_info_only:
            self._paint_meshes()
        else:
            self.attrs['painted'] = False
            
        

    def __getstate__(self):
        state = dict(b=self.b,k_edges=self.k_edges,attrs=self.attrs)
        return state

    def __setstate__(self,
                     state):
        self.attrs = state['attrs']
        self.k_edges = state['k_edges']
        self.b = state['b']

    def _paint_meshes(self):
        """Paint mesh(es), to prevent explicit re-painting with each bispectrum measurement
        """
        
        cfield = self.mesh.paint(mode='complex')
        self.mesh = nbk.FieldMesh(cfield)
        if self.num_fields > 1:
            cfield = self.second.paint(mode='complex')
            self.second = nbk.FieldMesh(cfield)
        if self.num_fields > 2:
            cfield = self.third.paint(mode='complex')
            self.third = nbk.FieldMesh(cfield)
            
        self.attrs['painted'] = True
        
        
    def set_k_edges(self,
                    k_edges):
        """Set array of k bin edges to input array
        
        Parameters
        ----------
        k_edges : array_like
            array with shape (N_tri,6) of triangle bin edges, where second index is over
                [ k1_low k1_high k2_low k2_high k3_low k3_high ]
        """
        self.k_edges = k_edges

        
    def measure_bispectrum(self,
                           imin = None,
                           imax = None,
                           kmeas_min = None,
                           kmeas_max = None,
                           out_file = None,
                           verbose = 0,
                           meas_type = 'full'):
        """Measure bispectrum over range of triangle bins.
        
        The range of triangle bins to use is defined by their minimum and maximum indices
        in the stored `k_edges` array. In the future, the user will be able to specify min
        and max k values, but that feature has not been implemented yet.
        
        The meas_type flag controls what information is computed and output:
            - 'full': bispectrum values, triangle counts in bins, and k_means in bins
            - 'grid_info': triangle counts in bins and k_means in bins
            - 'unnorm_b_value': unnormalized bispectrum values (not divided by N_triangles in bin)
        
        This routine uses the slower, less memory-intensive version of the FFT-based
        bispectrum algorithm, which is slower because it executes many redundant computations
        in determining the k_means and triangle counts. The measure_bispectrum_faster
        and measure_gridinfo_faster routines are recommended if measuring a large number
        of triangle bins, or measuring from multiple snapshots with the same grid size
        and binning scheme.
        
        Parameters
        ----------
        imin, imax : int
            minimum and maximum triangle bin indices to measure bispectrum between
        kmeas_min, kmeas_max : float, optional
            NOT IMPLEMENTED
        out_file : str, optional
            filename of output file to save intermittent results to
        verbose : int, optional
            verbosity level of output: 0 (default) for nothing, >0 for more than nothing
        meas_type : {'full', 'grid_info', 'unnorm_b_value'}
            measurement type, as described above
        
        Returns
        -------
        dict
            b dict (see doctring for FFTBispectrum class) containing information for
            bins looped over in this function call
        """

        # Make list of triangle bins to measure over
        if (kmeas_min is not None) and (kmeas_max is not None):
            raise NotImplementedError('kmeas_min/kmeas_max not implemented yet!')
        else:
            tri_list = self.k_edges[imin:min(imax,len(self.k_edges))]
            index_list = range(imin,min(imax,len(self.k_edges)))
            
        comm = CurrentMPIComm.get()
        
        unit_fac = self.attrs['pos_units_mpcoverh']
                
        # Main loop over triangle list
        Bi = 0.
        Ntri = 0.
        k1m = 0.
        k2m = 0.
        k3m = 0.
        for i,tri in enumerate(tri_list):
            if comm.rank == 0 and verbose > 0:
                print('Starting triangle %d' % i)
                
            # Measure quantities on given triangle bin, converting units of k bin edges
            # if necessary
            if meas_type == 'full':
                Bi,Ntri,k1m,k2m,k3m \
                    = bk_FFT_full(self.mesh,tri_list[i][0:2]*unit_fac,
                                       tri_list[i][2:4]*unit_fac,
                                       tri_list[i][4:6]*unit_fac,
                                       second_mesh = self.second,third_mesh = self.third)
            elif meas_type == 'grid_info':
                Ntri,k1m,k2m,k3m \
                    = bk_FFT_grid_info(self.mesh,tri_list[i][0:2]*unit_fac,
                                       tri_list[i][2:4]*unit_fac,
                                       tri_list[i][4:6]*unit_fac)
            elif meas_type == 'unnorm_b_value':
                Bi = bk_FFT_unnormalized_value(self.mesh,
                                               tri_list[i][0:2]*unit_fac,
                                               tri_list[i][2:4]*unit_fac,
                                               tri_list[i][4:6]*unit_fac,
                                               second_mesh = self.second,third_mesh = self.third)
                
            if comm.rank == 0 and verbose > 0:
                print('Finished triangle %d' % i)
                
            # Convert mean k values of each side to h Mpc^-1, and also convert units
            # of B value to h^-1 Mpc^6
            kmean = [k1m,k2m,k3m]
            for j in range(len(kmean)):
                kmean[j] /= unit_fac
            Bi *= unit_fac**6.
            
            # Add measurement info to stored arrays of previous measurements
            if self.b is None:
                self.b = { 'index': np.array([index_list[i]]),
                              'k_edge': np.array([tri_list[i]]), 
                              'k_mean': np.array([kmean]),
                              'B': np.array([Bi]),
                              'N_tri': np.array([Ntri]) }
            else:
                self.b['index'] = np.append(self.b['index'],[index_list[i]],axis=0)
                self.b['k_edge'] = np.append(self.b['k_edge'],[tri_list[i]],axis=0)
                self.b['k_mean'] = np.append(self.b['k_mean'],[kmean],axis=0)
                self.b['B'] = np.append(self.b['B'],[Bi],axis=0)
                self.b['N_tri'] = np.append(self.b['N_tri'],[Ntri],axis=0)
            
            # If out_file is specified, write measurement to disk, by appending
            # to file
            if comm.rank == 0:
                if verbose > 0:
                    print('FFT Bk: %d %g %g %g %g %g' % (i,self.b['k_mean'][-1][0],
                                                         self.b['k_mean'][-1][1],
                                                         self.b['k_mean'][-1][2],
                                                         self.b['B'][-1], 
                                                         self.b['N_tri'][-1]) )
                if out_file is not None:
                    with open(out_file, "a") as f:
                        if meas_type == 'full':
                            f.write('%d %e %e %e %e %e %e %e %e %e %e %e\n' % \
                                    (index_list[i],kmean[0],kmean[1],kmean[2], \
                                    tri_list[i][0],tri_list[i][1],tri_list[i][2], \
                                    tri_list[i][3],tri_list[i][4],tri_list[i][5], \
                                    Bi,Ntri))
                        elif meas_type == 'grid_info':
                            f.write('%d %e %e %e %e %e %e %e %e %e %e\n' % \
                                    (index_list[i],kmean[0],kmean[1],kmean[2], \
                                    tri_list[i][0],tri_list[i][1],tri_list[i][2], \
                                    tri_list[i][3],tri_list[i][4],tri_list[i][5], \
                                    Ntri))
                        elif meas_type == 'unnorm_b_value':
                            f.write('%d %e\n' % (index_list[i],Bi))

        # Return dict of measurements made in this function call
        b_return = { 'index': self.b['index'][-len(tri_list):],
                      'k_edge': self.b['k_edge'][-len(tri_list):], 
                      'k_mean': self.b['k_mean'][-len(tri_list):],
                      'B': self.b['B'][-len(tri_list):],
                      'N_tri': self.b['N_tri'][-len(tri_list):] }
        return b_return
    
    
    
    def measure_bispectrum_faster(self,
                           imin = None,
                           imax = None,
                           kmeas_min = None,
                           kmeas_max = None,
                           out_file = None,
                           verbose = 0):
        """Measure unnormalized bispectrum over range of triangle bins.
        
        The range of triangle bins to use is defined by their minimum and maximum indices
        in the stored `k_edges` array. In the future, the user will be able to specify min
        and max k values, but that feature has not been implemented yet.
        
        This routine only measures the unnormalized bispectrum values in each bin, and
        only stores the bin indices, bin bounds and unnormalized B values to the `b` dict
        afterwards. To obtain k_means and triangle counts (necessary for properly
        normalizing the B values from this routine), one should have run the
        measure_gridinfo_faster routine first and stored its output.
        
        Parameters
        ----------
        imin, imax : int
            minimum and maximum triangle bin indices to measure bispectrum between
        kmeas_min, kmeas_max : float, optional
            NOT IMPLEMENTED
        out_file : str, optional
            filename of output file to save intermittent results to
        verbose : int, optional
            verbosity level of output: 0 (default) for nothing, >0 for more than nothing
        meas_type : {'full', 'grid_info', 'unnorm_b_value'}
            measurement type, as described above
        
        Returns
        -------
        dict
            b dict (see doctring for FFTBispectrum class) containing information for
            bins looped over in this function call
        """

        num_tri_to_write = 100
        
        comm = CurrentMPIComm.get()
        
        # Generate lists of k bin edges and indices
        k_bin_edges = generate_bin_edge_list(self.attrs['kmin'],
                           self.attrs['kmax'],
                           self.attrs['dk'],
                           self.attrs['num_lowk_bins'],
                           self.attrs['dk_high'])
        k_bin_indices = np.arange(len(k_bin_edges))
        
        # Also fetch arrays of triangle bin edges and indices
        triangle_bin_edge_list = self.k_edges[imin:min(imax,len(self.k_edges))]
        triple_index_list = self.k_indices[imin:min(imax,len(self.k_indices))]
        triangle_index_list = range(imin,min(imax,len(self.k_indices)))
        
        # Make empty dict to hold delta fields
        delta_fields = {}
        if verbose > 0 and comm.rank == 0: print('Starting delta_field computations')
            
        # Generate real-space delta fields for each k bin, and add to dict, indexed
        # by k bin index
        for i in k_bin_indices:
            
            def mask(k,v):
                # Compute |\vec{k}| values on the mesh
                kk = sum(ki ** 2. for ki in k)**0.5
                # Mask out values outside desired bin
                return v * ((kk<=k_bin_edges[i][1]) & (kk>=k_bin_edges[i][0]))
        
            if verbose > 0 and comm.rank == 0: print('\tConstructing masked field for bin %d' % i)
                
            if self.second is not None:
                raise NotImplementedError('Only auto bispectrum implemented so far!')
        
            mesh_copy = self.mesh.view()
            mesh_copy._actions.append(('complex', mask, 'wavenumber'))
            delta_fields[i] = mesh_copy.paint(mode='real')

        if verbose > 0 and comm.rank == 0: print('\tStarting sum over triangles')
            
        # Now, loop over triangle bins
        B = 0.
        Ntri = 0.
        k1m = 0.
        k2m = 0.
        k3m = 0.
        for t,triple in enumerate(triple_index_list):
            
            # Compute unnormalized B value in triangle bin, using pre-computed delta fields
            # retrieved from dict
            B_local = np.sum(delta_fields[triple[0]]*delta_fields[triple[1]]*delta_fields[triple[2]]) \
                        * self.mesh.attrs['BoxSize'].prod()**2 \
                        / np.prod(self.mesh.attrs['Nmesh']) # / Nbin
            B_sum = np.sum(comm.gather(B_local))
            B = comm.bcast(B_sum)
        
            # Convert units of B value to h^-1 Mpc^6
            B *= self.attrs['pos_units_mpcoverh']**6.
            
            # Add measurement info to stored arrays of previous measurements
            if self.b is None:
                self.b = { 'index': np.array([triangle_index_list[t]]),
                              'k_edge': np.array([triangle_bin_edge_list[t]]), 
                              'k_mean': np.array([0.]),
                              'B': np.array([B]),
                              'N_tri': np.array([Ntri]) }
            else:
                self.b['index'] = np.append(self.b['index'],[triangle_index_list[t]],axis=0)
                self.b['k_edge'] = np.append(self.b['k_edge'],[triangle_bin_edge_list[t]],axis=0)
                self.b['k_mean'] = np.append(self.b['k_mean'],[0.],axis=0)
                self.b['B'] = np.append(self.b['B'],[B],axis=0)
                self.b['N_tri'] = np.append(self.b['N_tri'],[Ntri],axis=0)
            
            # If out_file is specified, write measurement to disk by appending to file.
            # Perform this writing operation for every num_tri_to_write triangles
            # at a time, to lessen I/O
            if comm.rank == 0:
                if verbose > 0:
                    print('FFT Bk: %d %d %g %g %g %g %g' % (t,triangle_index_list[t],
                                                            self.b['k_mean'][-1][0],
                                                            self.b['k_mean'][-1][1],
                                                            self.b['k_mean'][-1][2],
                                                            self.b['B'][-1],
                                                             self.b['N_tri'][-1]) )
                if out_file is not None:
                    if (t > 0) and ((t+1) % num_tri_to_write == 0):
                        with open(out_file, "a") as f:
                            for tt in range(num_tri_to_write):
                                f.write('%d %e %e %e %e %e %e %e\n' 
                                        % (triangle_index_list[t-(num_tri_to_write-1)+tt],triangle_bin_edge_list[t-99+tt][0],
                                           triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][1],triangle_bin_edge_list[t-99+tt][2],
                                           triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][3],triangle_bin_edge_list[t-99+tt][4],
                                           triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][5],self.b['B'][-num_tri_to_write+tt]))

        # After loop has finished, write any remaining triangles that have been computed
        # but have not yet been written out
        if comm.rank == 0 and out_file is not None and (t > 0) and ((t+1) % num_tri_to_write != 0):
            with open(out_file, "a") as f:
                t_left = (t+1) % num_tri_to_write
                for tt in range(t_left):
                    f.write('%d %e %e %e %e %e %e %e\n' 
                                % (triangle_index_list[t-t_left+1+tt],triangle_bin_edge_list[t-t_left+1+tt][0],
                                   triangle_bin_edge_list[t-t_left+1+tt][1],triangle_bin_edge_list[t-t_left+1+tt][2],
                                   triangle_bin_edge_list[t-t_left+1+tt][3],triangle_bin_edge_list[t-t_left+1+tt][4],
                                   triangle_bin_edge_list[t-t_left+1+tt][5],self.b['B'][-t_left+tt]))


        # Return dict of measurements made in this function call
        b_return = { 'index': self.b['index'][-len(triangle_index_list):],
                      'k_edge': self.b['k_edge'][-len(triangle_index_list):], 
                      'k_mean': self.b['k_mean'][-len(triangle_index_list):],
                      'B': self.b['B'][-len(triangle_index_list):],
                      'N_tri': self.b['N_tri'][-len(triangle_index_list):] }
        return b_return
    
    
    def measure_gridinfo_faster(self,
                           imin = None,
                           imax = None,
                           kmeas_min = None,
                           kmeas_max = None,
                           out_file = None,
                           verbose = 0):
        """Compute bispectrum binning info over range of triangle bins.
        
        The range of triangle bins to use is defined by their minimum and maximum indices
        in the stored `k_edges` array. In the future, the user will be able to specify min
        and max k values, but that feature has not been implemented yet.
        
        This routine only measures the triangle counts in bins and k_means in bins, to be
        saved and later combined with the unnormalized bispectrum values from
        measure_bispectrum_faster into a full set of measurements.
        
        Parameters
        ----------
        imin, imax : int
            minimum and maximum triangle bin indices to measure bispectrum between
        kmeas_min, kmeas_max : float, optional
            NOT IMPLEMENTED
        out_file : str, optional
            filename of output file to save intermittent results to
        verbose : int, optional
            verbosity level of output: 0 (default) for nothing, >0 for more than nothing
        meas_type : {'full', 'grid_info', 'unnorm_b_value'}
            measurement type, as described above
        
        Returns
        -------
        dict
            b dict (see doctring for FFTBispectrum class) containing information for
            bins looped over in this function call
        """
        
        num_tri_to_write = 100
        
        comm = CurrentMPIComm.get()
        
        n_mesh = self.mesh.attrs['Nmesh']
        box_size = self.mesh.attrs['BoxSize']
        
        # Generate lists of k bin edges and indices
        k_bin_edges = generate_bin_edge_list(self.attrs['kmin'],
                           self.attrs['kmax'],
                           self.attrs['dk'],
                           self.attrs['num_lowk_bins'],
                           self.attrs['dk_high'])
        k_bin_indices = np.arange(len(k_bin_edges))
        
        # Also fetch arrays of triangle bin edges and indices
        triangle_bin_edge_list = self.k_edges[imin:min(imax,len(self.k_edges))]
        triple_index_list = self.k_indices[imin:min(imax,len(self.k_indices))]
        triangle_index_list = range(imin,min(imax,len(self.k_indices)))
        
        # Make empty dicts to hold number_fields and k_fields
        number_fields = {}
        kk_fields = {}
        if verbose > 0 and comm.rank == 0: 
            print('Starting number_field and/or k_field computations')
            
        # Generate real-space number and k fields for each k bin, and add to dict, indexed
        # by k bin index
        for i in k_bin_indices:
            if verbose > 0 and comm.rank == 0:
                print('\tConstructing number_fields and/or k_fields for bin %d' % i)
            
            number_fields[i] = number_field(box_size,n_mesh,
                                            k_bin_edges[i][0],k_bin_edges[i][1])
            kk_fields[i] = k_field(box_size,n_mesh,
                                    k_bin_edges[i][0],k_bin_edges[i][1],1.)


        if verbose > 0 and comm.rank == 0: print('\tStarting sum over triangles')
            
        # Now, loop over triangle bins
        B = 0.
        Ntri = 0.
        k_mean_local = np.zeros(3)
        k_mean_sum = np.zeros(3)
        k_mean = np.zeros(3)
        for t,triple in enumerate(triple_index_list):
            
            # Compute number of triangles in bin
            Nbin_local = np.sum(number_fields[triple[0]] \
                                *number_fields[triple[1]] \
                                *number_fields[triple[2]]) / np.prod(n_mesh)
            Nbin_sum = np.sum(comm.gather(Nbin_local))
            Nbin = comm.bcast(Nbin_sum)

            # Compute mean values of k1,k2,k3 over triangle bin
            k_mean_local[0] = np.sum(kk_fields[triple[0]] \
                               *number_fields[triple[1]] \
                               *number_fields[triple[2]]) / np.prod(n_mesh) / Nbin
            k_mean_local[1] = np.sum(number_fields[triple[0]] \
                               *kk_fields[triple[1]] \
                               *number_fields[triple[2]]) / np.prod(n_mesh) / Nbin
            k_mean_local[2] = np.sum(number_fields[triple[0]] \
                               *number_fields[triple[1]] \
                               *kk_fields[triple[2]]) / np.prod(n_mesh) / Nbin
            for j in range(3):
                k_mean_sum[j] = np.sum(comm.gather(k_mean_local[j]))
                k_mean[j] = comm.bcast(k_mean_sum[j])
                
            # Convert units of k_means
            k_mean /= self.attrs['pos_units_mpcoverh']
            
            # Add measurement info to stored arrays of previous measurements
            if self.b is None:
                self.b = { 'index': np.array([triangle_index_list[t]]),
                              'k_edge': np.array([triangle_bin_edge_list[t]]), 
                              'k_mean': np.array([k_mean]),
                              'B': np.array([0.]),
                              'N_tri': np.array([Nbin]) }
            else:
                self.b['index'] = np.append(self.b['index'],[triangle_index_list[t]],axis=0)
                self.b['k_edge'] = np.append(self.b['k_edge'],[triangle_bin_edge_list[t]],axis=0)
                self.b['k_mean'] = np.append(self.b['k_mean'],[k_mean],axis=0)
                self.b['B'] = np.append(self.b['B'],[B],axis=0)
                self.b['N_tri'] = np.append(self.b['N_tri'],[Nbin],axis=0)
            
            if comm.rank == 0 and verbose > 1:
                print('Computed Bk for triangle index %d' % t)
            
            # If out_file is specified, write measurement to disk by appending to file.
            # Perform this writing operation for every num_tri_to_write triangles
            # at a time, to lessen I/O
            if comm.rank == 0:
                if verbose > 0 and ((t+1) % 100 == 0):
                    print('FFT Bk: %d %d %g %g %g %g %g' % (t,triangle_index_list[t],
                                                            self.b['k_mean'][-1][0],
                                                            self.b['k_mean'][-1][1],
                                                            self.b['k_mean'][-1][2],
                                                            self.b['B'][-1],
                                                             self.b['N_tri'][-1]) )
                if out_file is not None:
                    if (t > 0) and ((t+1) % num_tri_to_write == 0):
                        with open(out_file, "a") as f:
                            for tt in range(num_tri_to_write):
                                f.write('%d %e %e %e %e %e %e %e %e %e %e\n' 
                                 % (triangle_index_list[t-(num_tri_to_write-1)+tt],
                                    self.b['k_mean'][-num_tri_to_write+tt][0],
                                    self.b['k_mean'][-num_tri_to_write+tt][1],
                                    self.b['k_mean'][-num_tri_to_write+tt][2],
                                    triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][0],
                                triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][1],triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][2],
                                triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][3],triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][4],
                                triangle_bin_edge_list[t-(num_tri_to_write-1)+tt][5],self.b['N_tri'][-num_tri_to_write+tt]))
                        if verbose > 0:
                            print('Wrote %d B values, starting with triangle index %d' % (num_tri_to_write,triangle_index_list[t-(num_tri_to_write-1)+tt]))

        # After loop has finished, write any remaining triangles that have been computed
        # but have not yet been written out
        if comm.rank == 0 and out_file is not None and (t > 0) and ((t+1) % num_tri_to_write != 0):
            with open(out_file, "a") as f:
                t_left = (t+1) % num_tri_to_write
                for tt in range(t_left):
                    f.write('%d %e %e %e %e %e %e %e %e %e %e\n' 
                                % (triangle_index_list[t-t_left+1+tt],
                                    self.b['k_mean'][-t_left+tt][0],
                                    self.b['k_mean'][-t_left+tt][1],
                                    self.b['k_mean'][-t_left+tt][2],
                                   triangle_bin_edge_list[t-t_left+1+tt][0],
                                   triangle_bin_edge_list[t-t_left+1+tt][1],triangle_bin_edge_list[t-t_left+1+tt][2],
                                   triangle_bin_edge_list[t-t_left+1+tt][3],triangle_bin_edge_list[t-t_left+1+tt][4],
                                   triangle_bin_edge_list[t-t_left+1+tt][5],self.b['N_tri'][-t_left+tt]))

        # Return dict of measurements made in this function call
        b_return = { 'index': self.b['index'][-len(triangle_index_list):],
                      'k_edge': self.b['k_edge'][-len(triangle_index_list):], 
                      'k_mean': self.b['k_mean'][-len(triangle_index_list):],
                      'B': self.b['B'][-len(triangle_index_list):],
                      'N_tri': self.b['N_tri'][-len(triangle_index_list):] }
        return b_return
    
    
    def save_bispectrum(self,
                        out_file):
        """Save bispectrum measurements to ASCII file.
        
        This routine takes the `b` dict, storing the bispectrum measurements and computed
        grid info in this FFTBispectrum objects, and write this information to disk.
        
        Parameters
        ----------
        out_file : str
            output file name
        """
        b = self.b
        if b is None:
            return
        
        header = 'k1_mean k2_mean k3_mean k1_low k1_high k2_low k2_high k3_low k3_high ' \
                    + '[all h Mpc^-1] B [h^-6 Mpc^6] N_tri'
        with open(out_file,"w") as f:
            np.savetxt(out_file,
                       zip(b['k_mean'][:,0],b['k_mean'][:,1],b['k_mean'][:,2],
                           b['k_edge'][:,0],b['k_edge'][:,1],b['k_edge'][:,2],
                           b['k_edge'][:,3],b['k_edge'][:,4],b['k_edge'][:,5],
                           b['B'],b['N_tri']),
                       header=header)
        
        return
    
        
