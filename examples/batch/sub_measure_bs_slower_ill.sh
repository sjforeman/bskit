#!/bin/bash -l
#PBS -l nodes=1:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=24:00:00
#PBS -N fftb_100_1_b_512_kf_15k_to_20k

###
# Example script (for CITA Sunnyvale cluster):
#
# Measuring bispectrum from specified density grid in bigfile format, using slower
# (but less memory-intensive) routine.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export SNAP_PREFIX=/mnt/raid-cita/sforeman/illustrisTNG/sjf_grids/TNG100-1_delta_b_512_z0.0.bigfile
export NGRID=512
export NGRIDCIC=2048
export LBOX=75.
export STARTI=15000
export ENDI=20000
export OUTFILE=/mnt/raid-cita/sforeman/illustrisTNG/bs/TNG100-1_delta_b_512_z0.0_bs_kf_${STARTI}to${ENDI}.dat

source activate nbodykit-env
cd $PBS_O_WORKDIR
source set_script_dir.sh
mpirun -np 16 python ${SCRIPT_DIR}measure_bs_slow.py bigfile_grid $SNAP_PREFIX $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE --k_bin_info 0.0418879 30. 0.0837758 -cic
