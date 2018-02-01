#!/usr/bin/env python3
import argparse
import asyncio
import datetime
from ftplib import FTP
import logging
import os
import pathlib
import re
import sys
import urllib.request
import shutil
import subprocess
import xml.etree.ElementTree as ET
import pandas as pd
from holtlib import slurm_modules
from holtlib import slurm_job
import time

# TODOS:
    # catch "ftplib.error_perm: 550 pathogen/Results/Eschericia_coli_Shigella: No such file or directory" in parse_genome_trackr
    # catch "TimeoutError: [Errno 60] Operation timed out" in parse_genome_trackr


class BadAccession(Exception):
    pass


def sra_runs_from_bioproject_accessions(bioproject_accs):
    sra_runs = []
    bioproject_uids = uids_from_accession(bioproject_accs, 'bioproject')
    biosample_uids = biosample_uids_from_bioproject_uids(bioproject_uids)
    biosamples = biosamples_from_biosample_uids(biosample_uids)
    for biosample in biosamples:
        sra_runs += biosample.get_sra_runs()
    return sra_runs


def sra_runs_from_biosample_accessions(biosample_accs):
    sra_runs = []
    biosample_uids = uids_from_accession(biosample_accs, 'biosample')
    biosamples = biosamples_from_biosample_uids(biosample_uids)
    for biosample in biosamples:
        sra_runs += biosample.get_sra_runs()
    return sra_runs


def sra_runs_from_sra_accessions(sra_accs):
    sra_runs = []
    sra_run_uids = uids_from_accession(sra_accs, 'sra')
    biosample_uids = biosample_uids_from_sra_run_uids(sra_run_uids)
    biosamples = biosamples_from_biosample_uids(biosample_uids)
    for biosample in biosamples:
        sra_runs += biosample.get_sra_runs(sra_accs)
    return sra_runs


def sra_runs_from_sra_run_accessions(sra_run_accs):
    sra_runs = sra_runs_from_sra_accessions(sra_run_accs)
    return [x for x in sra_runs if x.accession in sra_run_accs]


def sra_runs_from_sra_experiment_accessions(sra_experiment_accs):
    sra_runs = sra_runs_from_sra_accessions(sra_experiment_accs)
    return [x for x in sra_runs if x.experiment.accession in sra_experiment_accs]


def sra_runs_from_sra_sample_accessions(sra_sample_accs):
    sra_runs = sra_runs_from_sra_accessions(sra_sample_accs)
    return [x for x in sra_runs if x.sample.sra_sample_accession in sra_sample_accs]


def uids_from_accession(accessions, database):
    """
    Takes a list of accessions for any NCBI database, returns uids in no particular order.
    """
    if not isinstance(accessions, list):
        accessions = [accessions]
    # Format URL
    esearch_template_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=%s&term=%s'
    esearch_url = esearch_template_url % (database, '+OR+'.join(accessions))

    # Make GET request
    with urllib.request.urlopen(esearch_url) as esearch_response:
        esearch_xml = esearch_response.read()
        esearch_root = ET.fromstring(esearch_xml)
        uids = [x.text for x in esearch_root.findall('./IdList/Id')]
        if len(accessions) != len(uids):
            raise BadAccession
        return uids


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


def biosample_uids_from_sra_run_uids(sra_run_uids):
    # TO DO: if there are too many SRA run UIDs, we should probably do the following stuff in chunks (e.g. 1000 at a time).

    elink_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi' + \
                '?dbfrom=sra&db=biosample&id=' + ','.join(sra_run_uids)
    with urllib.request.urlopen(elink_url) as elink_response:
        elink_xml = elink_response.read()
        elink_root = ET.fromstring(elink_xml)
        for link_set_db in elink_root.findall('./LinkSet/LinkSetDb'):
            if link_set_db.find('./LinkName').text == 'sra_biosample':
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

    # Make sure all of the runs nested under this sample have the SRA sample ID.
    for biosample in biosamples:
        biosample.add_sra_sample_to_runs()

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


def get_multiple_run_warning(runs, run_type, biosample):
    return 'There were multiple ' + run_type + ' runs for sample ' + biosample.accession + \
           '. Only the most recent (' + runs[0].accession + ') was downloaded. These ' \
           'additional runs were ignored: ' + ', '.join(x.accession for x in runs[1:])


def construct_accession_validators(type_suffices):
    '''Generate regex and associate with appropriate NCBI SRA API functions'''
    # Construct validators
    validators = list()

    # DRP/DRS/DRX/DRR, ERP/ERS/ERX/ERR, SRP/SRS/SRX/SRR
    source_prefix = {'DR', 'ER', 'SR'}
    for data_type, type_suffix in type_suffices.items():
        # Generate regex
        prefix_string = '|'.join(source_prefix)
        sra_validator = re.compile(r'^(?:%s)%s[0-9]+$' % (prefix_string, type_suffix))

        # Record
        validators.append((data_type, sra_validator))

    # PRJ and SAMN
    bioproject_validator = re.compile(r'^PRJ[A-Z]{2}[0-9]+$')
    biosample_validator = re.compile(r'^SAMN[0-9]+$')

    validators.append((sra_runs_from_bioproject_accessions, bioproject_validator))
    validators.append((sra_runs_from_biosample_accessions, biosample_validator))

    return validators


def validate_accessions(input_accessions, validators, type_suffices):
    '''Take a list of accessions and sort them using regex validators'''
    # Return variable
    validated_accessions = {k: list() for k in type_suffices.keys()}

    # Sort
    for input_accession in input_accessions:
        for accession_type, validator in validators:
            if validator.match(input_accession):
                validated_accessions[accession_type].append(input_accession)
                break
        else:
            logging.error('Could not determine accession type for %s' % input_accession)

    return validated_accessions


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

    def add_sra_sample_to_runs(self):
        experiments = self.illumina_experiments + self.long_read_experiments + \
                      self.other_experiments
        runs = []
        for experiment in experiments:
            runs += experiment.runs

        for run in runs:
            run.sample = self

    def get_sra_runs(self, sra_accs=None):
        """
        This function returns the SRA run accessions for this BioSample. It will include both
        Illumina and long read SRA runs. If there are multiple runs in a category (e.g. more than
        one Illumina run), it only returns the most recent.

        It returns a list of SraRun objects
        """
        runs = []

        illumina_runs, long_read_runs, other_runs = [], [], []
        for experiment in self.illumina_experiments:
            illumina_runs += experiment.runs
        for experiment in self.long_read_experiments:
            long_read_runs += experiment.runs
        for experiment in self.other_experiments:
            other_runs += experiment.runs

        illumina_runs = sorted(illumina_runs, key=lambda x: x.published_date, reverse=True)
        long_read_runs = sorted(long_read_runs, key=lambda x: x.published_date, reverse=True)

        # If the user asked for SRA accessions explicitly, then give the exact ones they asked for.
        if sra_accs is not None:
            all_runs = illumina_runs + long_read_runs + other_runs
            return [r for r in all_runs if r.accession in sra_accs]

        if illumina_runs:
            runs.append(illumina_runs[0])
            if len(illumina_runs) > 1:
                logging.warning(get_multiple_run_warning(illumina_runs, 'Illumina', self))
        if long_read_runs:
            runs.append(long_read_runs[0])
            if len(long_read_runs) > 1:
                logging.warning(get_multiple_run_warning(long_read_runs, 'long read', self))
        if other_runs:
            logging.warning('There were runs associated with sample ' + self.accession + ' which '
                            'were neither Illumina reads nor long reads. They were ignored: ' +
                            ', '.join(x.accession for x in other_runs))
        return runs


class SraExperiment(object):
    def __init__(self, sra_experiment_xml):
        experiment = sra_experiment_xml.find('EXPERIMENT')
        self.accession = experiment.attrib.get('accession')
        self.alias = experiment.attrib.get('alias')
        design = experiment.find('DESIGN')
        self.biosample_accession = None
        sample = sra_experiment_xml.find('SAMPLE')
        if sample is not None:
            for external_id in sample.iter('EXTERNAL_ID'):
                if external_id.attrib.get('namespace') == 'BioSample':
                    self.biosample_accession = external_id.text
        if self.biosample_accession is None:
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
            run = SraRun(run_xml)
            self.runs.append(run)
            run.experiment = self

    def __repr__(self):
        return str(self.accession)

    def get_platform_short(self):
        if self.platform == 'ILLUMINA':
            return 'illumina'
        elif self.platform == 'OXFORD_NANOPORE':
            return 'nanopore'
        elif self.platform == 'PACBIO_SMRT':
            return 'pacbio'
        else:
            return 'other'


class SraRun(object):

    running_downloads = 0

    def __init__(self, sra_run_xml):
        self.accession = sra_run_xml.attrib.get('accession')
        self.sample = None
        self.experiment = None
        self.warnings = []
        self.error = None
        self.alias = sra_run_xml.attrib.get('alias')
        self.total_spots = int(sra_run_xml.attrib.get('total_spots'))
        self.total_bases = int(sra_run_xml.attrib.get('total_bases'))
        self.size = int(sra_run_xml.attrib.get('size'))
        self.published_date = sra_run_xml.attrib.get('published')
        statistics = sra_run_xml.find('Statistics')
        #print(statistics.attrib.get('nreads'))
        #self.read_file_count = int(statistics.attrib.get('nreads'))
        self.read_counts = []
        self.read_average_lengths = []
        self.read_stdevs = []
        for read_file in statistics.findall('Read'):
            self.read_counts.append(int(read_file.attrib.get('count')))
            self.read_average_lengths.append(float(read_file.attrib.get('average')))
            self.read_stdevs.append(float(read_file.attrib.get('stdev')))
        #assert len(self.read_counts) == self.read_file_count
        #assert len(self.read_average_lengths) == self.read_file_count
        #assert len(self.read_stdevs) == self.read_file_count

        self.attempts = 0
        self.max_attempts = 3
        self.output_fps = list()

    def __repr__(self):
        return self.sample.sra_sample_accession + '_' + \
               self.accession + '_' + self.experiment.platform

    def get_filename_base(self):
        return self.sample.accession + '_' + self.accession

    def rename_output(self):
        file_renames = []
        if self.experiment.platform == 'ILLUMINA':
            old_name_1 = self.accession + '_1.fastq.gz'
            old_name_2 = self.accession + '_2.fastq.gz'
            new_name_1 = '%s_1.fastq.gz' % self.get_filename_base()
            new_name_2 = '%s_2.fastq.gz' % self.get_filename_base()
            file_renames.append((old_name_1, new_name_1))
            file_renames.append((old_name_2, new_name_2))
        else:
            old_name = self.accession + '.fastq.gz'
            new_name = '%s_%s_%s.fastq.gz' % (self.sample.accession, self.accession, self.experiment.get_platform_short())
            file_renames.append((old_name, new_name))

        for old_name, new_name in file_renames:
            # Rename files and update final filepaths
            os.rename(old_name, new_name)
            self.output_fps.append(new_name)

    async def download_task(self, simultaneous_downloads):
        # Wait until there are free job slots
        logging.info('Starting %s' % self.accession)
        while SraRun.running_downloads >= simultaneous_downloads:
            await asyncio.sleep(10)

        # Consume a job slot when this job is executed
        SraRun.running_downloads += 1

        while self.attempts < self.max_attempts:
            # Increment download attempt
            self.attempts += 1

            # Construct command
            fastq_dump_template = 'fastq-dump --gzip %s'

            fastq_dump_args = list()
            if self.experiment.platform == 'ILLUMINA':
                fastq_dump_args.append('--split-3')
            fastq_dump_args.append('--readids')
            fastq_dump_args.append(self.accession)

            fastq_dump_cmd = fastq_dump_template % ' '.join(fastq_dump_args)

            # Get a subprocess coroutine and execute
            logging.info('Downloading %s' % self.accession)


            process = await asyncio.create_subprocess_shell(fastq_dump_cmd,
                            stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)

            # Wait for subprocess to complete; blocking of current task
            stdout, stderr = await process.communicate()
            returncode = process.returncode

            # Check returncode; if we were successful break loop else attempt
            # to rescue
            if returncode == 0:
                break
            else:
                # TODO: check if error is resulting from a connect error
                # Example, don't retry if there was an issue other than a connection error:
                connection_errors = [
                    'connection busy while validating within network system module',
                    'timeout exhausted while reading file within network system module',
                    'Reading information from the socket failed',
                    'transfer interrupted while reading file within network system module',
                    'SSL - Connection requires a read call',
                    'item not found while constructing within virtual database module'
                    ]
                if not any(e in stderr.decode() for e in connection_errors):
                    self.error = 'fastq-dump error: %s' % stderr.decode()
                    break


        # File renaming/ moving synchronous but this won't be an issue
        # ...unless we're somehow renaming/ moving files between filesystems
        try:
            self.rename_output()
        except OSError as e:
            self.error = 'Unable to move file (%s)' % e

        # Release job count and return
        SraRun.running_downloads -= 1


def validate(date_text):
    '''
    Check that date specified to script is valid
    '''
    try:
        datetime.datetime.strptime(date_text, '%Y-%m-%d')
    except ValueError:
        logging.error('Date supplied (%s) for GenomeTrackr is in inccorrect format, should be YYYY-MM-DD.' % date_text)
        raise ValueError("Incorrect data format, should be YYYY-MM-DD")


# check that the table isn't empty
def check_dataframe_status(data_frame):
    if data_frame.empty:
        logging.error('GenomeTrackr table is empty, please check that any subsetting commands given are correct.')

def parse_genome_trackr(species, date, genome_trackr_col, genome_trackr_col_value):
    logging.info('Connecting to GenomeTrackr FTP server...')

    # Putting a try statement here as sometimes there are issues connecting to the FTP server
    try:
        ftp = FTP('ftp-trace.ncbi.nih.gov')
        ftp.login()
        # move to the genometrackr folder
        ftp.cwd('pathogen/Results/' + species)
        # get list of directories in species folder and identify most recent
        # most recent folder will be PDG directory with biggest number
        logging.info('Inspecting ' + species + ' directories...')
        dir_list = ftp.nlst()
        # intialise number to check against to see if we have biggest number
        check_num = 0
        correct_dir = ''
        for directory in dir_list:
            if directory.startswith('PDG'):
                folder_num = int(directory.split('.')[1])
                if folder_num > check_num:
                    # set check_num to be largest value we've seen
                    check_num = folder_num
                    # save the name of directory for later
                    correct_dir = directory
        # navigate into this directory and the Metadata folder
        logging.info('Most recent directory = ' + correct_dir)
        ftp.cwd(correct_dir + '/Metadata')
        # open a temp file to download to
        temp_tsv = 'genome_trackr_temp_' + str(os.getpid()) + '.tsv'
        genome_trackr_file = open(temp_tsv, 'wb')
        # download tsv file
        metadata_path = correct_dir + '.metadata.tsv'
        logging.info('Downloading metadata file ' + metadata_path + ' to ' + temp_tsv + '...')
        ftp.retrbinary('RETR ' + metadata_path, genome_trackr_file.write)
        # close the ftp connection and the genome trackr file
        ftp.close()
        genome_trackr_file.close()

        logging.info('Reading in GenomeTrackr data for ' + species + ' from ' + temp_tsv + ' ...')
        genome_trackr_table = pd.read_csv(temp_tsv, sep='\t')
        # if there's a date, extract any samples on or after that date
        if date:
            logging.info('Selecting only genomes that were added to the GenomeTrackr database on or after this date: %s' % date)
            genome_trackr_table = genome_trackr_table[genome_trackr_table['target_creation_date'] >= date]
            # check table not empty
            check_dataframe_status(genome_trackr_table)

        if genome_trackr_col:
            # extract only rows with column equaling value of interest
            ## TO DO: Make this more robust
            logging.info('Subsetting GenomeTrackr table on column %s with value of %s ...' % (genome_trackr_col, genome_trackr_col_value))
            genome_trackr_table = genome_trackr_table[genome_trackr_table[genome_trackr_col] == genome_trackr_col_value]
            # check table not empty
            check_dataframe_status(genome_trackr_table)
        # get a list of all the biosample accessions for the entries of interest
        genome_trackr_biosample_accessions = list(genome_trackr_table['biosample_acc'])
        logging.info('Have a list of %s Biosamples for download from GenomeTrackr.' % str(len(genome_trackr_biosample_accessions)))
        return genome_trackr_biosample_accessions

    except OSError:
        logging.info('Unable to connect to the GenomeTrackr FTP. You may want to try using a different computer/cluster or internet connection. Alternately, you can download the GenomeTrackr metadata file yourself from ftp-trace.ncbi.nih.gov/pathogen/Results/ and select your accessions of interest. These can then be passed to the script using --accession_list.')


def initialize_logging_file(logfile):
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        filemode='w',                      # write to log file (will overwrite on subsequent runs)
        format='%(asctime)s %(message)s')  # format the logfile

    # Output log to stdout as well as log file.
    root = logging.getLogger()
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
    root.addHandler(handler)

    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))  # print the terminal command


def check_fastq_dump_version():
    logging.info('Checking fastq-dump presence and version')
    if shutil.which('fastq-dump') is None:
        logging.error('Could not find fastq-dump')
        sys.exit()
    try:
        version_string = subprocess.check_output(['fastq-dump', '--version']).decode().strip()
        version_string = version_string.split(' : ')[1]
    except (subprocess.CalledProcessError, IndexError):
        logging.error('Could not parse "fastq-dump --version" output')
        sys.exit()
    current_url = 'https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current.version'
    with urllib.request.urlopen(current_url) as current_version_response:
        current_version = current_version_response.read().decode().strip()
    if version_string != current_version:
        logging.error('fastq-dump version (' + version_string +
                      ') is not the current version (' + current_version + ')')
        sys.exit()
    logging.info('fastq-dump version (' + version_string + ') is good')


###
# Argument parser
###
def get_arguments():

    parser = argparse.ArgumentParser(description='Download reads from NCBI')

    parser.add_argument('--accession_list', required=False, type=pathlib.Path,
                        help='File of accessions (one per line)')
    parser.add_argument('--accessions', required=False, type=str, nargs='+',
                        help='NCBI accessions (BioProject, BioSample or SRA, separated by spaces)')
    parser.add_argument('--logfile', default='download_reads.log',
                        help='Log file')
    parser.add_argument('--download_jobs', type=int, default=8,
                        help='Number of simultaneous downloads to run at once')
    parser.add_argument('--species', required=False, type=str, help='Species from GenomeTrackr. Must be one of the following: Acinetobacter, Campylobacter, Citrobacter_freundii, Elizabethkingia, Enterobacter, Escherichia_coli_Shigella, Klebsiella, Klebsiella_oxytoca, Kluyvera_intermedia, Legionella_pneumophila, Listeria, Morganella, Mycobacterium_tuberculosis, Neisseria, Providencia, Pseudomonas_aeruginosa, Salmonella, Serratia, Staphylococcus_pseudintermedius, Vibrio_parahaemolyticus')
    parser.add_argument('--date', required=False, type=str, help='Only required when downloading reads from GenomeTrackr. Will only download reads uploaded on or after the date specified in this argument. The corresponding column in the GenomeTrackr table is "target_creation_date". Date MUST be in the following format: YYYY-MM-DD.')
    parser.add_argument('--genome_trackr_col', required=False, type=str, help='Name of column in GenomeTrackr table which will be used to select only rows which equal a particular value - this value can be set with --genome_trackr_col_value.')
    parser.add_argument('--genome_trackr_col_value', required=False, type=str, help='Value in column of GenomeTrackr to use to select specific rows. Column name can be set with --genome_trackr_col.')
    # Ensure that input files exist if specified

    return parser.parse_args()


def main():

    ###
    # Initilization
    ###
    args = get_arguments()
    initialize_logging_file(args.logfile)
    check_fastq_dump_version()
    sra_runs = []

    # key: accession, value: reason for failure
    failed_acc = {}

    ###
    # Check arguments
    ###

    # do a check to make sure date is in the right format
    if args.date:
        validate(args.date)

    input_accessions = set()
    if args.accessions:
        input_accessions |= set(args.accessions)
        plural = '' if len(args.accessions) == 1 else 's'
        logging.info('Successfully read in ' + str(len(args.accessions)) + ' accession' + plural)
    if args.accession_list:
        logging.info('Reading in accession list file %s' % args.accession_list)
        with args.accession_list.open('r') as fh:
            input_accessions_from_file = {line.rstrip() for line in fh}

        plural = '' if len(input_accessions_from_file) == 1 else 's'
        logging.info('Successfully read in ' + str(len(input_accessions_from_file)) + ' accession' + plural + ' from file ' + str(args.accession_list))
        input_accessions |= input_accessions_from_file

    if args.accessions or args.accession_list:
        # Construct validators
        type_suffices = {sra_runs_from_bioproject_accessions: 'P',
                         sra_runs_from_biosample_accessions: 'S',
                         sra_runs_from_sra_experiment_accessions: 'X',
                         sra_runs_from_sra_run_accessions: 'R'}
        validators = construct_accession_validators(type_suffices)

        # Validate and sort input accessions
        validated_accessions = validate_accessions(input_accessions, validators, type_suffices)

        # Init SraRun objects from sorted accessions using appropriate function
        for func, accessions in validated_accessions.items():
            if accessions:
                sra_runs.extend(func(accessions))

    ###
    # Locate accessions for all reads from a GenomeTrackr species
    ###
    if args.species:
        # add to the list of sra objects for each biosample
        sra_runs += sra_runs_from_biosample_accessions(parse_genome_trackr(args.species, args.date, args.genome_trackr_col, args.genome_trackr_col_value))


    # Use asyncio to download reads in parallel.
    loop = asyncio.get_event_loop()
    async_futures = asyncio.gather(*[x.download_task(args.download_jobs) for x in sra_runs])
    logging.info('Downloading reads...')
    loop.run_until_complete(async_futures)
    loop.close()

    # TODO: Should check if experiment is actually RNA sequencing, not DNA, at least put this in the output file
    with open('accession_master_list.csv', 'w') as master_list:
        header = ('run_accession', 'experiment_accession', 'biosample_accession', 'library_source', 'seq_platform')
        print(*header, sep='\t', file=master_list)
        # Write out some data associated with each accession and any error(s)
        for sra_run in sra_runs:
            if sra_run.error:
                logging.info('Error downloading %s: %s', sra_run.accession, sra_run.error)
            else:
                data = (sra_run.accession, sra_run.experiment.accession, sra_run.sample.accession, sra_run.experiment.library_source, sra_run.experiment.platform)
                print(*data, sep='\t', file=master_list)

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
