#!/bin/bash

module load python3/3.7.4
module load parallel

MAXPROCS=8
OUTPUT_DIR=output

mv plots/*.png .
find $OUTPUT_DIR -name '*.hdf5' | parallel --bar --max-procs $MAXPROCS python3 plot_meshless.py 
mv *.png plots/
