"""Combine files of bispectrum binning info and unnormalized bispectrum values.

The most efficient way to measure the bispectrum on several snapshots is to
first compute the binning info (mean k_i values per bin, and number of triangles
per bin) a single time for a given box size and binning scheme, and then, for
each snapshot, compute unnormalized bispectrum values corresponding to this binning.
The final step, accomplished by this script, is to combine the binning info file
and the bispectrum value file into a single file that includes all binning info
and normalized (by the N_tri) bispectrum values.
"""

import numpy as np
import argparse

# Tolerance for comparing bin edges in different files to combine
tol = 0.01

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('bin_info_file',
                    help='file with precomputed k_mean and N_bin values')
parser.add_argument('b_file',
                    help='file with B values')
parser.add_argument('k_max',type=float,
                    help='k_max to use for processed file')
parser.add_argument('out_file',
                    help='output file')

# Parse arguments
args = parser.parse_args()
bin_info_file = args.bin_info_file
b_file = args.b_file
k_max = args.k_max
out_file = args.out_file

# Load files with bin info and unnormalized bispectrum measurements
bin_info = np.loadtxt(bin_info_file)
b_vals = np.loadtxt(b_file)

# Determine maximum bin index to process
k1min_arr = bin_info.T[4]
# Consider the k1_min values for all bins. If the highest such value is
# greater than the user-specified k_max, set i_from_k_max (the highest
# bin index to process) to the lowest index such that the highest k1_min
# is below k_max. (There's probably a nicer way to do it than what I've written
# below...)
if k1min_arr[-1] > k_max:
    i_from_k_max = np.where(k1min_arr == k1min_arr[k1min_arr > k_max][0])[0][0] - 1
else:
    i_from_k_max = len(k1min_arr)
i_max = np.min([len(bin_info),len(b_vals),i_from_k_max])

# Make empty array of combined files
full_arr = np.zeros((i_max,12))

# Loop over bispectrum bins we want to process
for i in range(i_max):
    
    # Compare bin indices in bin_info and b_vals files, and raise error if
    # mismatch is found
    if bin_info[i][0] != b_vals[i][0]:
        raise ValueError('Bin index mismatch between bin_info and B files: line %d: %d vs %d'
                         % (i,bin_info[i][0],b_vals[i][0]))
      
    # Also compare bin edges - make sure they're the same in both input files, or else
    # the output will be wrong
    for j in range(6):
        if np.abs((bin_info[i][4+j]-b_vals[i][1+j])/bin_info[i][4+j]) > tol:
            raise ValueError('Bin bound mismatch between bin_info and B files: line %d: %e vs %e'
                            % (i,bin_info[i][4+j],b_vals[i][1+j]))
    
    # Set the first 10 columns of the output file to those from the bin_info file
    full_arr[i,0:10] = bin_info[i,0:10]
    # Divide the unnormalized bispectrum value from the b_vals file by N_tri
    # from the bin_info file
    full_arr[i,10] = b_vals[i,7]/bin_info[i,10]
    # Set the final output column to N_tri in each bin
    full_arr[i,11] = bin_info[i,10]

# Save the processed file!
np.savetxt(out_file,full_arr,fmt='%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e',
          header='triangle index [dimensionless] k1_mean [h Mpc^-1] k2_mean [h Mpc^-1] k3_mean [h Mpc^-1] k1_min [h Mpc^-1] k1_max [h Mpc^-1] k2_min [h Mpc^-1] k2_max [h Mpc^-1] k3_min [h Mpc^-1] k3_max [h Mpc^-1] B_bin [h^-6 Mpc^6] N_modes')

print('Finished processing:')
print('\tout file: %s' % out_file)
print('\ti_max: %d' % i_max)
