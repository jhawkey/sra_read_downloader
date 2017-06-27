#!/usr/bin/env python3

import asyncio
import argparse
import logging
import pathlib
import sys

from holtlib import slurm_job
from holtlib import slurm_modules

###
# Argument parser
###
def get_arguments():

    parser = ArgumentParser(description='Download reads from NCBI')

    parser.add_argument()

    return parser.parse_args()

def main():

    ###
    # Initilization
    ###

    # get arguments
    args = get_arguments()

    # initialize logging file
    logging.basicConfig(
        filename=logfile, # name of log file
        level=logging.DEBUG, # set logging level to debug
        filemode='w', # write to log file (so will overwrite on subsequent runs)
        format='%(asctime)s %(message)s', # format the logfile
        datefmt='%m/%d/%Y %H:%M:%S') # format the date and time
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv))) # print the command sent to the command line

    # list of accessions to pass to fastq-dump
    acc_list = []

    # key: accession, value: reason for failure
    failed_acc = {}

    ###
    # Parse file with list of accession IDs
    ###
    if args.accession_list:
        # check the file exists
        if not args.accesion_list.exists():
            logging.error('Accession list file, %s, does not exist. Quitting.' % args.accession_list)
            # quit
            sys.exit(1)
        # read in file
        logging.info('Reading in accession list file %s ...' % args.accesion_list)

        # initilize counters
        skipped_acc = 0
        checked_acc = 0

        # open accession list file and parse
        with args.accesion_list.open() as f:
            # parse each line
            for line in f:
                acc_no = line.strip()
                # accessions should
                # check that the accession number has a valid prefix
                if not acc_no.startswith('SRR') or not acc_no.startswith('ERR') or not acc_no.startswith('DRR'):
                    logging.info('Accession %s is not valid (must being with SRR, ERR or DRR). Not downloading.' % acc_no)
                    skipped_acc += 1
                    failed_acc[acc_no] = 'Invalid accession number'
                else:
                    # add each accession to master list
                    acc_list.append(acc_no)
                    checked_acc += 1

    ###
    # Locate accessions for all reads in a project ID
    ###
    if args.project_no:
         pass
         # To consider:
         # - which set of Illumina reads to download if there are multiple
         # - which set of long reads to download if there are multiple
         # - if long reads should be downloaded at all (perhaps a flag that is automatically set to true?)
         # - in all cases, write out decision made to log file and record what other options could have occurred
         # - create an output file for each sample with multiple read sets for the user to inspect

    ###
    # Locate accessions for all reads from a GenomeTrackr species
    ###
    if args.genome_trackr:
        pass
        # To consider:
        # - only downloading reads uploaded after a specified date
        # - only downloading reads with a particular value in a column of the metadata table
        # - if the sample ID should be checked to see if there are long reads associatd with it
        # - if there are multiple entries per sample (don't want duplicates)
        # - should take a look at some metadata files and see what columns are present
        # - look at Zoe's script


    ###
    # Use SRA Toolkit to download each accession ID from acc_list
    ###

    # To do:
    # - use asyncio to launch only some jobs at once


    # Example command for fastq-dump (Illumina reads)
    # fastq-dump --split-3 --gzip --readids <acc>
    # Example command for fastq-dump (long reads)
    # fastq-dump --gzip --readids <acc>

    # Example of how to set up slurm job object
    '''
    # initalize job
    new_job = slurm_job.SlurmJob(job_name=NAME, partition='sysgen', time='0-01:00:00')
    # set up modules
    modules = []
    for module in modules:
        new_job.modules.append(slurm_modules.get_module('helix', module))
    # add commands
    new_job.commands.append(download_read)
    new_job.commands.append(remove_sra_file)

    # run job and write out script
    new_job.submit_sbatch_job()
    new_job.write_sbatch_script(reference.stem + '_jobscript.sh')
    '''

    ###
    # Clean up and write output files
    ###

    # To consider:
    # - a file containing all failed accessions with their reasons for failure
    # - a file containing all SAMPLES with multiple READ SETS and the IDs of the read sets that
    #   were not downloaded, as well as the read set that was downloaded (and whether these read
    #   sets are Illumina or long, etc)
    # - a file of all successfully downloaded accessions and their locations


if __name__ == '__main__':
    main()
