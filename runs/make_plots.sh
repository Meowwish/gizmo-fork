#!/bin/bash

MAXPROCS=8
OUTPUT_DIR=output_lowres

find $OUTPUT_DIR -name '*.hdf5' | parallel --bar --max-procs $MAXPROCS python3 plot_meshless.py 

