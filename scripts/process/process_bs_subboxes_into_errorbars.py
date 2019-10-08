"""Estimate sample-variance uncertainties on a measured bispectrum ratio using subboxes of the same simulation.

This script takes files produced by the measure_subbox_bs_fast.py script
and computes an estimate of the sample-variance uncertainty on the power spectrum ratio
in each k bin, based on taking the spread of the subbox measurements and scaling it by the
subbox-to-full-box volume ratio. The specific spread we compute is the 68% range centered
on the median, because this is not as sensitive to outliers as the standard deviation.

The script also multiplies this spread-based estimate by 2, to correct for the empirical
difference between subbox and full-box bispectrum measurements.
"""

import numpy as np
import argparse

# Multiplier that relates final bispectrum error estimates to spread in subbox measurements.
# For our hydro bispectrum paper, we set this multiplier to 2, based on comparisons between
# subbox measurements and measurements from independent simulations (see Appendix C).
bs_error_mult = 2.

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('subbox_bs_numerator_file_prefix',
                    help='file prefix (before subbox #) for numerator bispectrum measurements from subboxes')
parser.add_argument('subbox_bs_denominator_file_prefix',
                    help='file prefix (before subbox #) for denominator bispectrum measurements from subboxes')
parser.add_argument('in_file_suffix',
                    help='file suffix (after subbox #) for input bispectrum files')
parser.add_argument('n_sub',type=int,
                    help='number of subboxes')
parser.add_argument('out_file',
                    help='output file for estimated bispectrum ratio errorbars')

# Parse arguments
args = parser.parse_args()
subbox_bs_num_file_prefix = args.subbox_bs_numerator_file_prefix
subbox_bs_denom_file_prefix = args.subbox_bs_denominator_file_prefix
subbox_bs_file_suffix = args.in_file_suffix
n_sub = args.n_sub
out_file = args.out_file

# Load bin indices and B values from subbox bispectrum measurement files
bs_num = []
bs_denom = []
for i in range(n_sub):
    bs_num.append(np.loadtxt(subbox_bs_num_file_prefix+str(i)+subbox_bs_file_suffix)[:,(0,10)])
    bs_denom.append(np.loadtxt(subbox_bs_denom_file_prefix+str(i)+subbox_bs_file_suffix)[:,(0,10)])

# Convert to numpy arrays
bs_num = np.array(bs_num)
bs_denom = np.array(bs_denom)
   
# Make sure that input files have same lengths
if bs_num.shape != bs_denom.shape:
    raise Exception('Input files have different numbers of entries!')
    
# Divide B from numerator file by that from denominator file, and
# transpose ratio array to have shape (N_bin,2,N_sub) where 2 is our (index,B) columns
bs_ratio = bs_num.copy()
bs_ratio[:,:,1] /= bs_denom[:,:,1]
bs_ratio = np.transpose(bs_ratio,(1,2,0))

# For each k bin, compute central (around the median) the 68% range amongst the subboxes,
# divide by 2 (to convert to the analog of a standard deviation), and divide by sqrt(N_sub)
# as a way of accounting for the volume scaling. Finally, multiply by bs_error_mult.
b_iqr = np.array([ [bk[0,0],
                    bs_error_mult*(np.percentile(bk[1],84)-np.percentile(bk[1],16))/(2.*np.sqrt(n_sub))]
           for bk in bs_ratio])

# Save to file
np.savetxt(out_file,b_iqr,fmt='%d %.6e',
           header='triangle index [dimensionless] sig_ratio [dimensionless]')

print('Wrote output to %s' % out_file)
