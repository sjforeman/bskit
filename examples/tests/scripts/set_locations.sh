#!/bin/bash
#
# Script to set environment variables with locations of scripts
# and output, and amend PYTHONPATH with bskit location
#
# All paths are relative to bskit/examples/tests/scripts.
# If bskit is installed such that you don't have write access
# to the original grids or output_user directories, consider changing
# these paths to other directories you can write to.

export BSKIT_SCRIPT_DIR=../../../scripts/
export BSKIT_GRID_DIR=../grids/
export BSKIT_OUT_DIR=../output_user/

export PYTHONPATH=../../../:$PYTHONPATH
