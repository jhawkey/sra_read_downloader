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

def bioproject_uids_from_bioproject_accessions(bioproject_accs):
    esearch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi' + \
                  '?db=bioproject&term=' + ','.join(bioproject_accs)
    with urllib.request.urlopen(esearch_url) as esearch_response:
        esearch_xml = esearch_response.read()
        esearch_root = ET.fromstring(esearch_xml)
        return [x.text for x in esearch_root.findall('./IdList/Id')]


def biosample_uids_from_bioproject_uids(bioproject_uids):
    # TO DO: if there are too many BioProject UIDs, we should probably do the following stuff in chunks (e.g. 1000 at a time).

    elink_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi' + \
                '?dbfrom=bioproject&db=biosample&id=' + ','.join(bioproject_uids)
    with urllib.request.urlopen(elink_url) as elink_response:
        elink_xml = elink_response.read()
        elink_root = ET.fromstring(elink_xml)
        for link_set_db in elink_root.findall('./LinkSet/LinkSetDb'):
            if link_set_db.find('./LinkName').text == 'bioproject_biosample_all':
                return [x.text for x in link_set_db.findall('./Link/Id')]
    return []


def biosamples_from_biosample_uids(biosample_uids):
    # TO DO: if there are too many BioSample UIDs, we should probably do the following stuff in chunks (e.g. 1000 at a time).

    # First we build the BioSamples.
    biosamples = []
    efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' + \
                 '?dbfrom=biosample&db=biosample&id=' + ','.join(biosample_uids)

    # TO DO: Using '&retmode=json' would give me JSON results, which would be a lot nicer to parse!

    with urllib.request.urlopen(efetch_url) as efetch_response:
        efetch_xml = efetch_response.read()
        efetch_root = ET.fromstring(efetch_xml)
        for biosample_xml in efetch_root.findall('./BioSample'):
            biosamples.append(BioSample(biosample_xml))

    # Then we build the SRA experiments that are linked to those biosamples.
    elink_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi' + \
                '?dbfrom=biosample&db=sra&id=' + ','.join(b.uid for b in biosamples)
    with urllib.request.urlopen(elink_url) as elink_response:
        elink_xml = elink_response.read()
        elink_root = ET.fromstring(elink_xml)
        sra_experiment_uids = [x.text for x in elink_root.findall('./LinkSet/LinkSetDb/Link/Id')]
        sra_experiments = sra_experiments_from_sra_experiment_uids(sra_experiment_uids)

    # Now we have to associate SRA experiments with BioSamples.
    biosample_dict = {b.accession: b for b in biosamples}
    for experiment in sra_experiments:
        biosample = biosample_dict[experiment.biosample_accession]
        platform = experiment.platform.lower()
        if 'illumina' in platform:
            biosample.illumina_experiments.append(experiment)
        elif 'nanopore' in platform:
            biosample.long_read_experiments.append(experiment)
        elif 'pacbio' in platform:
            biosample.long_read_experiments.append(experiment)
        else:
            biosample.other_experiments.append(experiment)
    return biosamples


def sra_experiments_from_sra_experiment_uids(sra_experiment_uids):
    # TO DO: if there are too many SRA UIDs, we should probably do the following stuff in chunks (e.g. 1000 at a time).

    efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' + \
                 '?dbfrom=sra&db=sra&id=' + ','.join(sra_experiment_uids)
    sra_experiments = []
    with urllib.request.urlopen(efetch_url) as efetch_response:
        efetch_xml = efetch_response.read()
        efetch_root = ET.fromstring(efetch_xml)
        for sra_experiment_xml in efetch_root.findall('./EXPERIMENT_PACKAGE'):
            sra_experiments.append(SraExperiment(sra_experiment_xml))
    return sra_experiments


def get_sra_run_accession_for_biosamples(biosamples):
    sra_run_accessions = []
    sra_run_to_sample_dict = {}
    sra_run_warnings = []
    for biosample in biosamples:
        run_accessions, naming_dict, warnings = biosample.get_sra_run_accessions()
        sra_run_accessions += run_accessions
        sra_run_to_sample_dict.update(naming_dict)
        sra_run_warnings += warnings
    return sra_run_accessions, sra_run_to_sample_dict, sra_run_warnings


def get_multiple_run_warning_message(runs, run_type, biosample):
    return 'There were multiple ' + run_type + ' runs for sample ' + biosample.accession + \
           '. Only the most recent (' + runs[0].accession + ') was downloaded. These ' \
           'additional runs were ignored: ' + ', '.join(x.accession for x in runs[1:])


class BioSample(object):
    def __init__(self, biosample_xml):

        self.uid = biosample_xml.attrib.get('id')
        self.accession = biosample_xml.attrib.get('accession')
        self.submission_date = biosample_xml.attrib.get('submission_date')
        self.last_update = biosample_xml.attrib.get('last_update')
        self.sra_sample_accession = None
        for id_node in biosample_xml.iter('Id'):
            if id_node.attrib.get('db') == 'SRA':
                self.sra_sample_accession = id_node.text
        description = biosample_xml.find('Description')
        self.title = description.find('Title').text
        organism = description.find('Organism').attrib
        self.taxonomy_id = organism.get('taxonomy_id')
        self.taxonomy_name = organism.get('taxonomy_name')
        self.illumina_experiments = []
        self.long_read_experiments = []
        self.other_experiments = []

    def __repr__(self):
        biosample_repr = str(self.accession) + ' (' + self.taxonomy_name
        if self.illumina_experiments and self.long_read_experiments:
            biosample_repr += ', hybrid'
        elif self.illumina_experiments:
            biosample_repr += ', Illumina'
        elif self.long_read_experiments:
            biosample_repr += ', long read'
        else:
            biosample_repr += ', unknown'
        biosample_repr += ')'
        return biosample_repr

    def get_sra_run_accessions(self):
        """
        This function returns the SRA run accessions for this BioSample. It will include both
        Illumina and long read SRA runs. If there are multiple runs in a category (e.g. more than
        one Illumina run), it only returns the most recent.

        It returns:
          * a list of SRA run accessions
          * a dictionary of SRA run accessions -> SRA sample accessions (for naming files)
          * a list of warning messages
        """
        run_accessions = []
        run_to_sample_dict = {}
        warnings = []

        illumina_runs, long_read_runs, other_runs = [], [], []
        for experiment in self.illumina_experiments:
            illumina_runs += experiment.runs
        for experiment in self.long_read_experiments:
            long_read_runs += experiment.runs
        for experiment in self.other_experiments:
            other_runs += experiment.runs

        illumina_runs = sorted(illumina_runs, key=lambda x: x.published_date, reverse=True)
        long_read_runs = sorted(long_read_runs, key=lambda x: x.published_date, reverse=True)

        if illumina_runs:
            illumina_run = illumina_runs[0]
            run_accessions.append(illumina_run)
            run_to_sample_dict[illumina_run] = self.sra_sample_accession
            if len(illumina_runs) > 1:
                warnings.append(get_multiple_run_warning_message(illumina_runs, 'Illumina', self))

        if long_read_runs:
            long_read_run = long_read_runs[0]
            run_accessions.append(long_read_run)
            run_to_sample_dict[long_read_run] = self.sra_sample_accession
            if len(long_read_runs) > 1:
                warnings.append(get_multiple_run_warning_message(illumina_runs, 'long read', self))

        if other_runs:
            warnings.append('There were runs associated with sample ' + self.accession + ' which '
                            'were neither Illumina reads nor long reads. They were ignored: ' +
                            ', '.join(x.accession for x in other_runs))

        return run_accessions, run_to_sample_dict, warnings


class SraExperiment(object):
    def __init__(self, sra_experiment_xml):
        experiment = sra_experiment_xml.find('EXPERIMENT')
        self.accession = experiment.attrib.get('accession')
        self.alias = experiment.attrib.get('alias')
        design = experiment.find('DESIGN')
        self.biosample_accession = None
        for external_id in design.iter('EXTERNAL_ID'):
            if external_id.attrib.get('namespace') == 'BioSample':
                self.biosample_accession = external_id.text
        library_descriptor = design.find('LIBRARY_DESCRIPTOR')
        self.library_name = library_descriptor.find('LIBRARY_NAME').text
        self.library_strategy = library_descriptor.find('LIBRARY_STRATEGY').text
        self.library_source = library_descriptor.find('LIBRARY_SOURCE').text
        self.library_selection = library_descriptor.find('LIBRARY_SELECTION').text
        self.library_layout = library_descriptor.find('LIBRARY_LAYOUT')[0].tag
        platform_node = experiment.find('PLATFORM')[0]
        self.platform = platform_node.tag
        self.instrument_model = platform_node.find('INSTRUMENT_MODEL').text
        self.runs = []
        for run_xml in sra_experiment_xml.findall('RUN_SET/RUN'):
            self.runs.append(SraRun(run_xml))

    def __repr__(self):
        return str(self.accession)


class SraRun(object):
    def __init__(self, sra_run_xml):
        self.accession = sra_run_xml.attrib.get('accession')
        self.alias = sra_run_xml.attrib.get('alias')
        self.total_spots = int(sra_run_xml.attrib.get('total_spots'))
        self.total_bases = int(sra_run_xml.attrib.get('total_bases'))
        self.size = int(sra_run_xml.attrib.get('size'))
        self.published_date = sra_run_xml.attrib.get('published')
        statistics = sra_run_xml.find('Statistics')
        self.read_file_count = int(statistics.attrib.get('nreads'))
        self.read_counts = []
        self.read_average_lengths = []
        self.read_stdevs = []
        for read_file in statistics.findall('Read'):
            self.read_counts.append(int(read_file.attrib.get('count')))
            self.read_average_lengths.append(float(read_file.attrib.get('average')))
            self.read_stdevs.append(float(read_file.attrib.get('stdev')))
        assert len(self.read_counts) == self.read_file_count
        assert len(self.read_average_lengths) == self.read_file_count
        assert len(self.read_stdevs) == self.read_file_count

    def __repr__(self):
        return str(self.accession)
if __name__ == '__main__':
    main()
