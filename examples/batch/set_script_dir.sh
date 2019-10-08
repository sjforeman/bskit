#!/bin/bash
#
# One-liner to set environment variable with directory containing scripts,
# to be called by batch queue submission scripts (so that each submission script
# doesn't need to be changed if the directory structure is changed)

export SCRIPT_DIR=/mnt/raid-cita/sforeman/bs/bskit/scripts/