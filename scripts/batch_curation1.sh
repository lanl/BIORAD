#!/bin/sh
#
# Perform machine learning tests for the data curation paper
# All clustering performed using the "--cluster.align" flag (first time only)

NUM_TRIAL=100

for TARGET in all pdbbind units no_kd different_heterodimer
do
    # Note that random forest predictions are pulled to the mean of the regression output
    ../src/biorad -t ${NUM_TRIAL} --affinity data/curation_${TARGET}.csv \
        --pdb.ignore-charge --pdb ../data/PDBBind_hydrogen \
        --cluster cluster/cluster.align.${TARGET} \
        -o output_hetatm/curation_${TARGET}.out \
        --hist.min 0.0 \
        --hist.max 12.0 \
        --hist.bins 1
        
        # Adjacency features
        #--hist.min 0.1 \
        #--hist.max 20.0 \
        #--hist.bins 25
done
