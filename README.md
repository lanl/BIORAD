# BIORAD
Source code and input data files for the manuscript "The Impact of Curation Errors in the PDBBind Database on Machine Learning Predictions of Protein-Protein Binding Affinity" by Gans et al.

This repository contains the necessary computer code and input files to perform the following tasks:
1) Download Open Access full-text publications from PubMed Central for protein heterodimer PDB records in the PDBBind.
2) Download the PDB files for the protein heterodimer PDBBind records in PDBBind.
3) Add hydrogens to the downloaded PDB files using the [Open Babel](https://openbabel.org/docs/Introduction/intro.html) program.
4) Train and test a simple random forest-based machine learning algorithm for predicting the equilibrium dissociation constant (actually log<sub>10</sub>(K<sub>D</sub>) ) from 3D protein structure coordinates.

## Downloading Open Access full-text publication from PubMed Central
- This step is only required for obtaining a local copy of the full-text publications associated with the protein heterodimer interactions in the PDBBind. This information was used to manually curate the extraction of K<sub>D</sub> values from the scientific literature.
- Run the python script `scripts/download_oa.py` using the provided CSV file of PDB and PubMed Central accessions (`data/pdbbind_oa.csv`) to download all of the corresponding files from PubMed Central. Using the output directory specified by the user (`-d <output directory>`), this script will create a separate subdirectory (labeled by PDB accession) to store the full-text information associated with each PDB accession. Note that the same manuscript may be associated with more than one PDB accession and, as a result, will be downloaded multiple times.

## Download the PDB files for the protein heterodimer PDBBind records


## Build the random forest machine learning algorithm for predicting protein-protein binding affinity
- This program is implemented in C++ and has the following dependencies that must be installed before the program can be compiled:
  -  [MPI](https://www.open-mpi.org/) (message passing interface). OpenMPI was used to develop and test the code, but any standard flavor of MPI should work.
  -  The [GSL - Gnu Scientific Library](https://www.gnu.org/software/gsl/). Please edit the provided `Makefile` to specify the locations of the GSL `include` directory and the library files (`libgsl.a` and `libgslcblas.a`)
  -  [optional] A C++ compiler that supports OpenMP
    -  If your C++ compiler does not support OpenMP (which is the case for the default OSX C++ compiler), please edit the provided `Makefile` to comment out the OpenMP flag; `OPENMP = -fopenmp` &rarr; `OPENMP = #-fopenmp`
- After the dependencies have been satisfied, run the `make` command from the `src` directory.
- 
