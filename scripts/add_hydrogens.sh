#!/usr/bin/sh

# Use open babel to add hydrogens to PDB files
INPUT_DIR="heterodimer_pdb"
OUTPUT_DIR="heterodimer_pdb_hydrogen"

# Create the output directory if it does not already exist
mkdir -p $OUTPUT_DIR

for INPUT_FILE in ${INPUT_DIR}/*.pdb
do
    OUTPUT_FILE=$(basename $INPUT_FILE)

    echo "Adding hydrogens to $OUTPUT_FILE"
    
    obabel $INPUT_FILE -O ${OUTPUT_DIR}/${OUTPUT_FILE} -h
done
