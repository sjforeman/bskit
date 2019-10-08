#!/bin/bash -l
#PBS -l nodes=1:ppn=16
#PBS -q sandyq
#PBS -r n
#PBS -l walltime=12:00:00
#PBS -N fftb_ill_100-1_cb_bDM_z0

###
# Example script (for CITA Sunnyvale cluster):
#
# Measuring cross bispectrum between two different density grids in bigfile format,
# on equilateral triangles only, using slower measurement routine.
###

export OMP_NUM_THREADS=1

export PYTHONUSERBASE=/mnt/raid-cita/$USER/.local
export PYTHONPATH=/mnt/raid-cita/$USER/.local/lib/python2.7/site-packages:/mnt/raid-cita/$USER:$PYTHONPATH
export PATH="/mnt/raid-cita/sforeman/anaconda2/bin:$PATH"

export Z=0.0

export SNAP_PREFIX_1=/mnt/raid-cita/sforeman/illustrisTNG/sjf_grids/TNG100-1_delta_b_1024_z${Z}.bigfile
export SNAP_PREFIX_2=/mnt/raid-cita/sforeman/illustrisTNG/sjf_grids/TNG_DM100-1_delta_m_1024_z${Z}.bigfile
export NGRID=1024
export NGRIDCIC=2048
export LBOX=75.
export STARTI=0
export ENDI=200
export OUTFILE=/mnt/raid-cita/sforeman/illustrisTNG/bs/TNG100-1_delta_b_cross_delta_DMO_1024_z${Z}_bseq_kf_6kf_40lowbins.dat

export DK=0.0837758
export DK_HIGH=0.50265
export KMIN=0.0418879
export KMAX=23.5
export NUM_LOWK_BINS=40

source activate nbodykit-env
cd $PBS_O_WORKDIR
source set_script_dir.sh
mpirun -np 16 python ${SCRIPT_DIR}measure_bs_slow.py bigfile_grid $SNAP_PREFIX_1 $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE --k_bin_info $KMIN $KMAX $DK -cic --high_k_bin_info \
$NUM_LOWK_BINS $DK_HIGH  --triangle_type equilateral --snap_prefix2 $SNAP_PREFIX_2