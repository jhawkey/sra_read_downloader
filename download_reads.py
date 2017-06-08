#!/usr/bin/env python3

import asyncio
import argparse
import logging
import pathlib

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

    # parse a list of accessions
    if args.accession_list:
        pass

    # download all reads in a project
    if args.project_no:
         pass

    # download all reads from a genome trackr species
    if args.genome_trackr:
        pass

# Example command for fastq-dump
# fastq-dump --split-3 --gzip --readids <acc>.sra 

# Example on how to set up slurm job object
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

if __name__ == '__main__':
    main()
