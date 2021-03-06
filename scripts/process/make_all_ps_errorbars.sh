#!/bin/sh
#
# Script to process all power spectrum subbox measurements into sample-variance errorbars


BASE_DIR=/mnt/raid-cita/sforeman/

TNG_DIR=illustrisTNG/ps/
ILL_DIR=original_illustris/ps_orig/
EAGLE_DIR=eagle/ps_eagle/
BAH_DIR=bahamas/ps_bahamas/

EBAR_DIR=processed_errorbars/

N_SUB=8

for SIM_DIR in $TNG_DIR $ILL_DIR $EAGLE_DIR $BAH_DIR
do
    DIR=$BASE_DIR$SIM_DIR$EBAR_DIR
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
        
        python process_ps_subboxes_into_errorbars.py ${BASE_DIR}${TNG_DIR}TNG300-1_delta_${PART}_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat ${BASE_DIR}${TNG_DIR}TNG_DM300-1_delta_m_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat $N_SUB ${BASE_DIR}${TNG_DIR}${EBAR_DIR}TNG300-1_delta_${PART}_1024_z${Z}.0_ps_kf_0.5kf_errorbars.dat
        
        python process_ps_subboxes_into_errorbars.py ${BASE_DIR}${TNG_DIR}TNG100-1_delta_${PART}_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat ${BASE_DIR}${TNG_DIR}TNG_DM100-1_delta_m_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat $N_SUB ${BASE_DIR}${TNG_DIR}${EBAR_DIR}TNG100-1_delta_${PART}_1024_z${Z}.0_ps_kf_0.5kf_errorbars.dat
        
        python process_ps_subboxes_into_errorbars.py ${BASE_DIR}${ILL_DIR}Illustris-1_delta_${PART}_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat ${BASE_DIR}${ILL_DIR}Illustris-1-Dark_delta_m_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat $N_SUB ${BASE_DIR}${ILL_DIR}${EBAR_DIR}Illustris-1_delta_${PART}_1024_z${Z}.0_ps_kf_0.5kf_errorbars.dat
        
        python process_ps_subboxes_into_errorbars.py ${BASE_DIR}${EAGLE_DIR}Eagle100_delta_${PART}_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat ${BASE_DIR}${EAGLE_DIR}Eagle_DM100_delta_m_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat $N_SUB ${BASE_DIR}${EAGLE_DIR}${EBAR_DIR}Eagle100_delta_${PART}_1024_z${Z}.0_ps_kf_0.5kf_errorbars.dat

        python process_ps_subboxes_into_errorbars.py ${BASE_DIR}${BAH_DIR}Bahamas_delta_${PART}_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat ${BASE_DIR}${BAH_DIR}Bahamas_DM2_delta_m_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat $N_SUB ${BASE_DIR}${BAH_DIR}${EBAR_DIR}Bahamas_delta_${PART}_1024_z${Z}.0_ps_kf_0.5kf_errorbars.dat
        
        python process_ps_subboxes_into_errorbars.py ${BASE_DIR}${BAH_DIR}Bahamas_hiAGN_delta_${PART}_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat ${BASE_DIR}${BAH_DIR}Bahamas_DM2_delta_m_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat $N_SUB ${BASE_DIR}${BAH_DIR}${EBAR_DIR}Bahamas_hiAGN_delta_${PART}_1024_z${Z}.0_ps_kf_0.5kf_errorbars.dat
                
        python process_ps_subboxes_into_errorbars.py ${BASE_DIR}${BAH_DIR}Bahamas_lowAGN_delta_${PART}_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat ${BASE_DIR}${BAH_DIR}Bahamas_DM2_delta_m_1024_z${Z}.0_8subboxes_ps_kf_0.5kf.dat $N_SUB ${BASE_DIR}${BAH_DIR}${EBAR_DIR}Bahamas_lowAGN_delta_${PART}_1024_z${Z}.0_ps_kf_0.5kf_errorbars.dat

    done
done