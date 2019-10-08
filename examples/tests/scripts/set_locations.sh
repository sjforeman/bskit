#!/bin/bash
#
# Script to set environment variables with locations of scripts
# and output, and amend PYTHONPATH with bskit location

export SCRIPT_DIR=../../scripts/
export GRID_DIR=../grids/
export OUT_DIR=../output_user/

export PYTHONPATH=../../../:$PYTHONPATH