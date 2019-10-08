"""Transform bispectrum bin info file to a new box size.

If the bispectrum binning scheme is specified using multiples of the
fundamental wavenumber k_f = 2pi/Lbox, then to translate a previously-computed
bin info file from one box size to another, the k values in the file
merely need to be multiplied by the ratio of the original and new
box side lengths. This script applies this operation to a bin info file.
"""

import numpy as np
import argparse

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('in_file',
                    help='input bin info file to transform')
parser.add_argument('out_file',
                    help='output bin info file')
parser.add_argument('old_Lbox',type=float,
                    help='Lbox of input fule')
parser.add_argument('new_Lbox',type=float,
                    help='desired Lbox')

# Parse arguments
args = parser.parse_args()
in_file = args.in_file
out_file = args.out_file
old_Lbox = args.old_Lbox
new_Lbox = args.new_Lbox

bin_info = np.loadtxt(in_file)

for i in range(1,10):
    bin_info.T[i] *= old_Lbox/new_Lbox
    
np.savetxt(out_file,bin_info,fmt='%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e')