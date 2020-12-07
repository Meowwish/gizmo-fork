#!/bin/bash

MAXPROCS=4
OUTPUT_DIR=output_agora_medres_virial

find $OUTPUT_DIR -name '*.0.hdf5' | parallel --bar --max-procs $MAXPROCS python3 plot_meshless.py 

