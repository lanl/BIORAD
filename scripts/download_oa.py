#!/usr/bin/env python3

import sys
import os
import argparse
import time
from urllib.request import urlopen
from urllib.parse import urlparse
from ftplib import FTP
import subprocess
import xml.etree.ElementTree as ET

VERSION = '0.1'
NCBI_OA_URL = 'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id='
SLEEP_DURATION = 1.0

def main():

    # Command line arguments (via argparse)
    parser = argparse.ArgumentParser( 
        description='Download full text OpenAccess records v.{}'.
        format(VERSION),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter )

    parser.add_argument('-i', dest='input_file', help='Input CSV file of PDB accessions (with matching PMC IDs) to download')
    parser.add_argument('-d', dest='output_directory', help='Output directory to store downloaded text')
    parser.add_argument('--clobber', dest='clobber', action='store_true', help='Overwrite existing output (default is to skip)')

    args = parser.parse_args()
    
    if not args.input_file:
        sys.stderr.write('Please specify a CSV input file of PDB accessions and PMC IDs to read (-i)\n')
        sys.exit(0)
    
    if not args.output_directory:
        sys.stderr.write('Please specify an output directory (-d)\n')
        sys.exit(0)

    try:
        fin = open(args.input_file)
    except OSError:
        sys.stderr.write('Unable to open {} for reading\n'.format(args.input_file))
        sys.exit(0)

    header = []
    pdb_to_pmcid = dict()
    pdb_col = -1
    pmcid_col = -1

    for line in fin:
        
        line = line.strip()

        # Skip comments
        if '#' in line:
            line = line[:line.find('#')]

        data = line.split(',')

        if len(data) == 0:
            continue

        if len(header) == 0:
            header = data

            pdb_col = header.index('PDB')
            pmcid_col = header.index('PMCID')

        else:

            if (len(data[pdb_col]) != 0) and (len(data[pmcid_col]) != 0):

                if data[pdb_col] in pdb_to_pmcid:
                    sys.stderr.write('Found duplicate PDB accession: {}\n'.format(data[pdb_col]))
                    sys.exit(1)
                
                pdb_to_pmcid[ data[pdb_col] ] = data[pmcid_col]

    fin.close()

    print('Found a total of {} PDB records will associated PMC ID values'.format(len(pdb_to_pmcid)), file=sys.stderr)

    if os.path.exists(args.output_directory):
        if not os.path.isdir(args.output_directory):
            sys.stderr.write('Please specify a valid output directory\n')
            sys.exit(0)
    else:
        os.makedirs(args.output_directory)

    for pdb,pmcid in pdb_to_pmcid.items():

        output_dir = os.path.join(args.output_directory, pdb)

        download_full_text(pdb, pmcid, output_dir, args.clobber)

    print('Downloading is complete!', file=sys.stderr)

def download_full_text(m_pdb, m_pmcid, m_output_dir, m_clobber):

    print('Downloading {} ({}) to {}'.format(m_pdb, m_pmcid, m_output_dir), file=sys.stderr)

    # If the output directory does not exist, create it now
    if not os.path.exists(m_output_dir):
        os.makedirs(m_output_dir)
    elif not m_clobber:
        sys.stderr.write('Output directory for {} already exists, skipping download\n'.format(m_pdb))
        return

    try:
        response = urlopen(NCBI_OA_URL + m_pmcid).read()
    except:
        sys.stderr.write('\nUnable to download text for {}\n'.format(m_pdb))
        return

    tree = ET.fromstring(response)

    # How many files did we download?
    num_download = 0

    for records in tree.findall('records'):

        for record in records.findall('record'):

            for link in record.findall('link'):

                format = link.get('format')
                href = link.get('href')

                url = urlparse(href)

                ftp_hostname = url.hostname
                ftp_full_path = url.path
                ftp_filename = os.path.split(ftp_full_path)[1]
                output_filename = os.path.join(m_output_dir, ftp_filename)

                print('\tStoring {} as {}'.format(href, output_filename), file=sys.stderr)
                
                ftp = FTP(ftp_hostname) # connect to host, default port
                ftp.login() # user anonymous, passwd anonymous@

                if format == 'pdf':

                    ftp.retrbinary('RETR ' + ftp_full_path, open(output_filename, 'wb').write)
                elif format == 'tgz':

                    ftp.retrbinary('RETR ' + ftp_full_path, open(output_filename, 'wb').write)

                    # Untar the file to the output directory
                    subprocess.run(['tar', '-xf', output_filename, '-C', m_output_dir])
                else:
                    sys.stderr.write('\tUnknown file format: {}\n'.format(format))
                
                ftp.quit()
                num_download += 1
                time.sleep(SLEEP_DURATION) # Don't overwhelm the NCBI servers


    if num_download == 0:
        sys.stderr.write('\tWarning! Failed to download any documents!\n')

main()