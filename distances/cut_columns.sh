#!/bin/bash

# input directory
indir="4_column"
# output directory
outdir="2_column"

# make sure output directory exists
mkdir -p "$outdir"

# loop through all files in the input directory
for f in "$indir"/*; do
    # get just the filename (without path)
    fname=$(basename "$f")
    
    # cut the first two columns and save in output dir
    cut -f3,4 "$f" > "$outdir/$fname"
done
