# BIORAD
Source code and input data files for the manuscript "The Impact of Curation Errors in the PDBBind Database on Machine Learning Predictions of Protein-Protein Binding Affinity" by Gans et al.
If you run into any problems, please describe the problem on the GitHub BIORAD issues tab and we'll get back to you as soon as we can.

This repository contains the necessary computer code and input files to perform the following tasks:
1) Download Open Access full-text publications from PubMed Central for protein heterodimer PDB records in the PDBBind.
2) Download the PDB files for the protein heterodimer PDBBind records in PDBBind.
3) Add hydrogens to the downloaded PDB files using the [Open Babel](https://openbabel.org/docs/Introduction/intro.html) program.
4) Train and test a simple random forest-based machine learning algorithm for predicting the equilibrium dissociation constant (actually log<sub>10</sub>(K<sub>D</sub>) ) from 3D protein structure coordinates.

## Downloading Open Access full-text publication from PubMed Central
- This step is only required for obtaining a local copy of the full-text publications associated with the protein heterodimer interactions in the PDBBind. This information was used to manually curate the extraction of K<sub>D</sub> values from the PubMed Central Open Access scientific literature.
- Run the python script `scripts/download_oa.py` using the provided CSV file of PDB and PubMed Central accessions (`data/pdbbind_oa.csv`) to download all of the corresponding files from PubMed Central. Using the output directory specified by the user (`-d <output directory>`), this script will create a separate subdirectory (labeled by PDB accession) to store the full-text information associated with each PDB accession. Note that the same manuscript may be associated with more than one PDB accession and, as a result, will be downloaded multiple times.

## Download the PDB files for the protein heterodimer PDBBind records
- The PDB files for each protein heterodimer record in the PDBBind can be downloaded using the `pdb_batch_download.sh` shell script [available](https://www.rcsb.org/scripts/batch_download.sh) from the RCSB PDB. This script accepts a comma-separated list of PDB accessions (using the `-f`) to download. The file `data/pdbbind_heterodimer_pdb.csv` contains a comma-separated list of PDB accessions for the protein heterodimer records in the PDBBind.
  - After obtaining the `pdb_batch_download.sh` script, running the command `pdb_batch_download.sh -f data/pdbbind_heterodimer_pdb.csv -o heterodimer_pdb` (from the root directory of the BIORAD repository) will download all of the protein heterodimer PDB files to the directory `heterodimer_pdb`. The name of the destination directory for storing the PDB files, `heterodimer_pdb`, is suggested, but not required.

## Add hydrogen atoms to the downloaded PDB files
- After downloading all of the PDB files for the protein heterodimer records in the PDBBind, hydrogen atoms can be added using the [Open Babel](https://openbabel.org/docs/Introduction/intro.html) program.
- After downloading and installing the Open Babel program, the `scripts/add_hydrogens.sh` script can be used to add hydrogen atoms. This script may be edited to specify the input directory of downloaded PDB files (e.g., `heterodimer_pdb`) and the desired output directory for storing the modified PDB files (e.g., `heterodimer_pdb_hydrogen`)

## Build the `biorad` software for predicting protein-protein binding affinity
- The `bioread` program is implemented in C++, using a Linux development environment (Ubuntu), and has the following dependencies that must be installed before the program can be compiled:
  -  MPI (message passing interface). [OpenMPI](https://www.open-mpi.org/) was used to develop and test the code, but any standard flavor of MPI should work.
  -  The [GSL - Gnu Scientific Library](https://www.gnu.org/software/gsl/). After installing and building the GSL, please edit the provided BIORAD `src/Makefile` to specify the locations of the GSL `include` directory and the library files (`libgsl.a` and `libgslcblas.a`)
  -  [optional] A C++ compiler that supports OpenMP
    -  If your C++ compiler does not support OpenMP (which is the case for the default OSX C++ compiler), please edit the provided `src/Makefile` to comment out the OpenMP flag; `OPENMP = -fopenmp` &rarr; `OPENMP = #-fopenmp`
- After the dependencies have been satisfied, run the `make` command from the `src` directory. This will create the `biorad` executable.

# Running the BIORAD software to train and test random forest models for predicting protein-protein binding affinity
- Running the script `scripts/batch_curration1.sh` in the root directory of the BIORAD respository will automatically run the `biorad` program on the different subsets of PDBBind and manually curated data.
  - This script assumes that the input PDB files are stored in the `heterodimer_pdb_hydrogen` directory. The script will need to be modified if this directory name is not being used.
  - This script will attempt to write output files (summarizing the machine learning test results) to a directory called `output`. The script creates the `output` directory if it does not already exist.
  - Please note that output files generated by `scripts/batch_curration1.sh` (with slightly different directory names) are already included in the `output` directory. If you don't want to overwrite these files, they should be renamed.
- Running the `biorad -h` command will display the set of allowed command line arguments for the `biorad` program.