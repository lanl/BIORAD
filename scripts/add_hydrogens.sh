#!/usr/bin/sh

# Use open babel to add hydrogens to PDB files
INPUT_DIR="heterodimer_pdb"
OUTPUT_DIR="heterodimer_pdb_hydrogen"

for INPUT_FILE in ${INPUT_DIR}/*.pdb
do
    OUTPUT_FILE=$(basename $INPUT_FILE)

    echo "Adding hydrogens to $OUTPUT_FILE"
    
    obabel $INPUT_FILE -O ${OUTPUT_DIR}/${OUTPUT_FILE} -h
done
