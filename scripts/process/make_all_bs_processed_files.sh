#!/bin/sh
#
# Script to process all bispectrum bin info and unnormalized-value files into
# files containing bin info and normalized values


BASE_DIR=/mnt/raid-cita/sforeman/

BIN_INFO_DIR=bs/ipy/bs_grid_info/

TNG_DIR=illustrisTNG/bs/
ILL_DIR=original_illustris/bs_orig/
EAGLE_DIR=eagle/bs_eagle/
BAH_DIR=bahamas/bs_bahamas/

PROC_DIR=paper_NEW/

for SIM_DIR in $TNG_DIR $ILL_DIR $EAGLE_DIR $BAH_DIR
do
    DIR=$BASE_DIR$SIM_DIR$PROC_DIR
    if [ -d $DIR ]
    then
        echo "Directory "$DIR" exists"
    else
        mkdir $DIR
        echo "Created directory "$DIR
    fi
done

for Z in 0 1 2 3
do
    for PART in m dm b
    do
        echo "--------------------------"
        echo "z = "$Z", "$PART
        echo "--------------------------"
        
        python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox205_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${TNG_DIR}TNG300-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 7.8 ${BASE_DIR}${TNG_DIR}${PROC_DIR}TNG300-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat
        
        python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox75_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${TNG_DIR}TNG100-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 21.5 ${BASE_DIR}${TNG_DIR}${PROC_DIR}TNG100-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat
        
        python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox75_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${ILL_DIR}Illustris-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 21.5 ${BASE_DIR}${ILL_DIR}${PROC_DIR}Illustris-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat

        python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox67.7_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${EAGLE_DIR}Eagle100_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 23.8 ${BASE_DIR}${EAGLE_DIR}${PROC_DIR}Eagle100_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat

        python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox400_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${BAH_DIR}Bahamas_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 4.0 ${BASE_DIR}${BAH_DIR}${PROC_DIR}Bahamas_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat

        python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox400_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${BAH_DIR}Bahamas_lowAGN_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 4.0 ${BASE_DIR}${BAH_DIR}${PROC_DIR}Bahamas_lowAGN_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat

        python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox400_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${BAH_DIR}Bahamas_hiAGN_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 4.0 ${BASE_DIR}${BAH_DIR}${PROC_DIR}Bahamas_hiAGN_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat
                
    done
done

PART=m
for Z in 0 1 2 3
do
    echo "--------------------------"
    echo "z = "$Z", DMO"
    echo "--------------------------"

    python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox205_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${TNG_DIR}TNG_DM300-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 7.8 ${BASE_DIR}${TNG_DIR}${PROC_DIR}TNG_DM300-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat
    
    python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox75_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${TNG_DIR}TNG_DM100-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 21.5 ${BASE_DIR}${TNG_DIR}${PROC_DIR}TNG_DM100-1_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat

    python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox75_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${ILL_DIR}Illustris-1-Dark_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 21.5 ${BASE_DIR}${ILL_DIR}${PROC_DIR}Illustris-1-Dark_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat

    python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox67.7_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${EAGLE_DIR}Eagle_DM100_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 23.8 ${BASE_DIR}${EAGLE_DIR}${PROC_DIR}Eagle_DM100_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat

    python process_fast_bs_measurement.py ${BASE_DIR}${BIN_INFO_DIR}Lbox400_1024_kf_6kf_40lowkbins.dat ${BASE_DIR}${BAH_DIR}Bahamas_DM2_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat 4.0 ${BASE_DIR}${BAH_DIR}${PROC_DIR}Bahamas_DM2_delta_${PART}_1024_z${Z}.0_bs_kf_6kf_40lowbins.dat
                
done