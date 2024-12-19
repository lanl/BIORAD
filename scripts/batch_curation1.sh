#!/bin/sh
#
# Perform machine learning tests for the data curation manuscript
# All clustering performed using the "--cluster.align" flag (first time only)

# Create the output directory if it does not already exist
mkdir -p output

NUM_TRIAL=100

for TARGET in all pdbbind units no_kd different_heterodimer
do
    src/biorad -t ${NUM_TRIAL} --affinity data/curation_${TARGET}.csv \
        --pdb.ignore-charge --pdb heterodimer_pdb_hydrogen \
        --cluster cluster/cluster.align.${TARGET} \
        -o output/curation_${TARGET}.out \
        --hist.min 0.0 \
        --hist.max 12.0 \
        --hist.bins 1
done
