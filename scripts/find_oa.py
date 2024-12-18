#!/usr/bin/env python3

import sys
import csv
import argparse

VERSION = '0.1'

def main():

    # Command line arguments (via argparse)
    parser = argparse.ArgumentParser( 
        description='Label PubMed IDs that correspond to OpenAccess records and add the PMC ID v.{}'.
        format(VERSION),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter )

    parser.add_argument('-i', dest='pdbbind_file', help='Input CSV file of PDB accessions with matching PubMed IDs')
    parser.add_argument('--inventory', dest='oa_file', help='Input inventory file of Open Access publications from https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv')
    
    args = parser.parse_args()

    if not args.pdbbind_file:
        sys.stderr.write('Please specify an input CSV file of PDB accession with matching PubMed IDs\n')
        sys.exit(0)

    if not args.oa_file:
        sys.stderr.write('Please specify an inventory file of Open Access publications\n')
        sys.exit(0)

    oa_pubmed_id = dict() # PMID -> PMCID
    header = []
    pubmed_id_col = -1

    with open(args.oa_file) as f:

        reader = csv.reader(f)

        for data in reader:

            if len(header) == 0:
                header = data

                pubmed_id_col = header.index('PMID')
                pubmed_central_id_col = header.index('Accession ID')
            else:

                if len(data) != len(header):
                    sys.stderr.write('Unexpected number of columns: {}\n'.format(line))
                    sys.exit(0)
                
                if (len(data[pubmed_id_col]) != 0) and (len(data[pubmed_central_id_col]) != 0):
                    oa_pubmed_id[ data[pubmed_id_col] ] = data[pubmed_central_id_col]


    print('Found {} OpenAccess PubMed IDs'.format(len(oa_pubmed_id)), file=sys.stderr)

    try:
        fin = open(args.pdbbind_file)
    except OSError:
        sys.stderr.write('Unable to open {} for reading\n'.format(pdbbind_file))
        sys.exit(0)

    header = []
    pubmed_id_col = -1
    num_oa_match = 0
    num_pdbbind = 0

    for line in fin:
         
        line = line.strip()
        data = line.split(',')

        if len(header) == 0:
            header = data
            pubmed_id_col = header.index('PubMedId')
        else:

            num_pdbbind += 1

            if (len(data[pubmed_id_col]) != 0) and (data[pubmed_id_col] in oa_pubmed_id):
                print('{},yes,{}'.format(line, oa_pubmed_id[ data[pubmed_id_col] ]))
                num_oa_match += 1
            else:
                print('{},no,'.format(line))

    fin.close()

    print('{}/{} PDBBind records have OpenAccess PubMed Ids'.format(num_oa_match, num_pdbbind), file=sys.stderr)

main()
