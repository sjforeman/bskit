"""Estimate sample-variance uncertainties on a measured power spectrum ratio using subboxes of the same simulation.

This script takes files produced by the measure_ps_from_bigfile_subboxes.py script
and computes an estimate of the sample-variance uncertainty on the power spectrum ratio
in each k bin, based on taking the spread of the subbox measurements and scaling it by the
subbox-to-full-box volume ratio. The specific spread we compute is the 68% range centered
on the median, because this is not as sensitive to outliers as the standard deviation.
"""

import numpy as np
import argparse

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('subbox_ps_numerator_file',
                    help='file with numerator power spectrum measurements from subboxes')
parser.add_argument('subbox_ps_denominator_file',
                    help='file with denominator power spectrum measurements from subboxes')
parser.add_argument('n_sub',type=int,
                    help='number of subboxes')
parser.add_argument('out_file',
                    help='output file for estimated power spectrum ratio errorbars')

# Parse arguments
args = parser.parse_args()
subbox_ps_num_file = args.subbox_ps_numerator_file
subbox_ps_denom_file = args.subbox_ps_denominator_file
n_sub = args.n_sub
out_file = args.out_file

# Load first 3 columns (subbox index, k, Pk) of subbox PS files
ps_num = np.loadtxt(subbox_ps_num_file)[:,:3]
ps_denom = np.loadtxt(subbox_ps_denom_file)[:,:3]

# Make sure that input files have same lengths
if ps_num.shape != ps_denom.shape:
    raise Exception('Input files have different numbers of entries!')

# Make sure that lengths of input files match specified number of subboxes
if (ps_num.shape[0] % n_sub != 0) or (ps_denom.shape[0] % n_sub != 0):
    raise Exception('Input number of subboxes disagrees with file lengths!')

# Divide Pk from numerator file by that from denominator file
ps_ratio = ps_num.copy()
ps_ratio[:,2] /= ps_denom[:,2]

# Reshape and transpose array to have shape (N_bins,3,N_sub) where 3 is for
# our (subbox index, k, Pk) columns
ps_ratio = ps_ratio.reshape((n_sub,ps_ratio.shape[0]/n_sub,3))
ps_ratio = np.transpose(ps_ratio,(1,2,0))

# For each k bin, compute central (around the median) the 68% range amongst the subboxes,
# divide by 2 (to convert to the analog of a standard deviation), and divide by sqrt(N_sub)
# as a way of accounting for the volume scaling
p_iqr = np.array([ [pk[1,0],(np.percentile(pk[2],84)-np.percentile(pk[2],16))/(2.*np.sqrt(n_sub))]
          for pk in ps_ratio])

# Save to file
np.savetxt(out_file,p_iqr,fmt='%.6e %.6e',
           header='k_mean_over_subbox_bin [h Mpc^-1] sig_ratio [dimensionless]')

print('Wrote output to %s' % out_file)
