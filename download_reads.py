#!/usr/bin/env python3
import argparse
import asyncio
import datetime
from ftplib import FTP
import logging
import pathlib
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET


import pandas as pd


from holtlib import slurm_job
from holtlib import slurm_modules


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
        sra_runs += biosample.get_sra_runs()
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
    # Ensure accession argument is as a list
    if not isinstance(accessions, list):
        accessions = [accessions]
    # Format URL
    esearch_template_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=%s&term=%s'
    esearch_url = esearch_template_url % (database, '+OR+'.join(accessions))

    # Make GET request
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
        self.warnings = []

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

    def get_sra_runs(self):
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

        if illumina_runs:
            illumina_run = illumina_runs[0]
            if len(illumina_runs) > 1:
                self.warnings.append(get_multiple_run_warning_message(illumina_runs, 'Illumina',
                                                                      self))
            runs.append(illumina_run)

        if long_read_runs:
            long_read_run = long_read_runs[0]
            runs.append(long_read_run)
            if len(long_read_runs) > 1:
                self.warnings.append(get_multiple_run_warning_message(illumina_runs, 'long read',
                                                                      self))

        if other_runs:
            self.warnings.append('There were runs associated with sample ' + self.accession +
                                 ' which were neither Illumina reads nor long reads. They were '
                                 'ignored: ' + ', '.join(x.accession for x in other_runs))
        return runs


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
        return self.sample.sra_sample_accession + '_' + \
               self.accession + '_' + self.experiment.platform

    def get_filename_base(self):
        return self.sample.sra_sample_accession + '_' + self.accession

    async def download_task(self, simultaneous_downloads):

        # Wait until there are free job slots
        while SraRun.running_downloads >= simultaneous_downloads:
            await asyncio.sleep(10)

        # Consume a job slot when this job is executed
        SraRun.running_downloads += 1

        fastq_dump_cmd = ['fastq-dump', '--gzip']
        if self.experiment.platform == 'ILLUMINA':
            fastq_dump_cmd.append('--split-3')
        fastq_dump_cmd += ['--readids', self.accession]

        print(fastq_dump_cmd)  # TEMP
        # TO DO: actually run fastq-dump

        file_renames = []
        if self.experiment.platform == 'ILLUMINA':
            old_name_1 = self.accession + '_1.fastq.gz'
            old_name_2 = self.accession + '_2.fastq.gz'
            new_name_1 = self.get_filename_base() + '_1.fastq.gz'
            new_name_2 = self.get_filename_base() + '_2.fastq.gz'
            file_renames.append((old_name_1, new_name_1))
            file_renames.append((old_name_2, new_name_2))
        else:
            old_name = self.accession + '.fastq.gz'
            new_name = (self.get_filename_base() + '_' + self.experiment.get_platform_short() +
                        '.fastq.gz')
            file_renames.append((old_name, new_name))

        for old_name, new_name in file_renames:
            print(old_name, '->', new_name)  # TEMP
            # TO DO: actually do the file renaming

        print()  # TEMP

        # Release job count and return
        SraRun.running_downloads -= 1


def validate(date_text):
    '''
    Check that date specified to script is valid
    '''
    try:
        datetime.datetime.strptime(date_text, '%Y-%m-%d')
    except ValueError:
        logging.error('Date supplied (%s) for GenomeTrackr is in inccorrect format, should be YYYY-MM-DD.' % args.date)
        raise ValueError("Incorrect data format, should be YYYY-MM-DD")


# check that the table isn't empty
def check_dataframe_status(data_frame):
    if data_frame.empty:
        logging.error('GenomeTrackr table is empty, please check that any subsetting commands given are correct.')


###
# Argument parser
###
def get_arguments():

    parser = argparse.ArgumentParser(description='Download reads from NCBI')

    parser.add_argument('--accession_list', required=False, type=pathlib.Path,
                        help='File of accessions (one per line)')
    parser.add_argument('--bioprojects', required=False, nargs='+',
                        help='NCBI BioProject accessions')
    parser.add_argument('--biosamples', required=False, nargs='+',
                        help='NCBI BioSample accessions')
    parser.add_argument('--sra', required=False, nargs='+',
                        help='NCBI SRA accessions')
    parser.add_argument('--genome_trackr', required=False, type=str,
                        help='GenomeTrackr species')
    parser.add_argument('--logfile', default='download_reads.log',
                        help='Log file')
    parser.add_argument('--download_jobs', type=int, default=8,
                        help='Number of simultaneous downloads to run at once')
    parser.add_argument('--species', required=False, type=str, help='Species of interest for GenomeTrackr.')
    parser.add_argument('--date', required=False, type=str, help='Only required when downloading reads from GenomeTrackr. Will only download reads uploaded on or after the date specified in this argument. The corresponding column in the GenomeTrackr table is "target_creation_date". Date MUST be in the following format: YYYY-MM-DD.')
    parser.add_argument('--genome_trackr_col', required=False, type=str, help='Name of column in GenomeTrackr table which will be used to select only rows which equal a particular value - this value can be set with --genome_trackr_col_value.')
    parser.add_argument('--genome_trackr_col_value', required=False, type=str, help='Value in column of GenomeTrackr to use to select specific rows. Column name can be set with --genome_trackr_col.')
    # Ensure that input files exist if specified

    return parser.parse_args()


def main():

    ###
    # Initilization
    ###

    # get arguments
    args = get_arguments()

    # initialize logging file
    logging.basicConfig(
        filename=args.logfile, # name of log file
        level=logging.DEBUG, # set logging level to debug
        filemode='w', # write to log file (so will overwrite on subsequent runs)
        format='%(asctime)s %(message)s', # format the logfile
        datefmt='%m/%d/%Y %H:%M:%S') # format the date and time
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv))) # print the command sent to the command line

    # list of SRA runs
    sra_runs = []

    # key: accession, value: reason for failure
    failed_acc = {}

    ###
    # Check arguments
    ###

    # do a check to make sure date is in the right format
    if args.date:
        validate(args.date)

    ###
    # Parse file with list of accession IDs
    ###
    if args.accession_list:
        # Get accessions
        logging.info('Reading in accession list file %s' % args.accession_list)
        with args.accession_list.open('r') as fh:
            input_accessions = {line.rstrip() for line in fh}

        logging.info('Successfully read in %s accesions from file %s' % (str(len(input_accessions)), args.accession_list))
        # Construct validators
        validators = list()

        # DRP/DRS/DRX/DRR, ERP/ERS/ERX/ERR, SRP/SRS/SRX/SRR
        source_prefix = {'DR', 'ER', 'SR'}
        type_suffices = {sra_runs_from_bioproject_accessions: 'P',
                         sra_runs_from_biosample_accessions: 'S',
                         sra_runs_from_sra_experiment_accessions: 'X',
                         sra_runs_from_sra_run_accessions: 'R'}

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

        # Validate and sort accessions
        validated_accessions = {k: list() for k in type_suffices.keys()}
        for input_accession in input_accessions:
            for accession_type, validator in validators:
                if validator.match(input_accession):
                    validated_accessions[accession_type].append(input_accession)
                    break
            else:
                logging.error('Could not determine accession type for %s' % input_accession)

        # Init SraRun objects from sorted accessions using appropriate function
        for func, accession in validated_accessions.items():
            sra_runs.extend(func(accessions))

    ###
    # Get SraRun objects from NCBI accessions
    ###
    if args.bioprojects:
        sra_runs += sra_runs_from_bioproject_accessions(args.bioprojects)
    if args.biosamples:
        sra_runs += sra_runs_from_biosample_accessions(args.biosamples)
    if args.sra:
        sra_run_accessions = [x for x in args.sra if x[2] == 'R']
        sra_experiment_accessions = [x for x in args.sra if x[2] == 'X']
        sra_sample_accessions = [x for x in args.sra if x[2] == 'S']
        if sra_run_accessions:
            sra_runs += sra_runs_from_sra_run_accessions(sra_run_accessions)
        if sra_experiment_accessions:
            sra_runs += sra_runs_from_sra_experiment_accessions(sra_experiment_accessions)
        if sra_sample_accessions:
            sra_runs += sra_runs_from_sra_sample_accessions(sra_sample_accessions)

    ###
    # Locate accessions for all reads from a GenomeTrackr species
    ###
    if args.genome_trackr:
        # login to FTP server
        ftp = FTP('ftp-trace.ncbi.nih.gov')
        ftp.login()
        # move to the genometrackr folder
        ftp.cwd('pathogen/Results')
        # TODO: get the correct species folder
        ftp.cwd()
        # get list of directories in species folder and identify most recent
        # most recent folder will be PDG directory with biggest number
        dir_list = ftp.nlst()
        # intialise number to check against to see if we have biggest number
        check_num = 0
        for directory in dir_list:
            if directory.startswith('PDG'):
                folder_num = int(directory.split('.')[1])
                if folder_num > check_num:
                    # set check_num to be largest value we've seen
                    check_num = folder_num
                    # save the name of directory for later
                    correct_dir = directory
        # navigate into this directory and the Metadata folder
        ftp.cwd(correct_dir + '/Metadata')
        # open a temp file to download to
        genome_trackr_file = open('genome_trackr_temp.tsv', 'w')
        # download tsv file
        ftp.retrbinary(correct_dir + '.metadata.tsv', genome_trackr_file.write)
        # close the ftp connection
        ftp.close()

        #genome_trackr_table = pd.read_csv('genome_trackr_metadata_15082016.tsv', sep='\t')
        logging.info('Reading in GenomeTrackr data for %s ...' % args.species)
        # TEMP
        genome_trackr_table = pd.read_csv(args.genome_trackr_file, sep='\t')
        # if there's a date, extract anything after the date
        if args.date:
            # take all entires on or after specified date
            logging.info('Selecting only genomes that were added to the GenomeTrackr database on or after this date: %s' % args.date)
            genome_trackr_table = genome_trackr_table[genome_trackr_table['target_creation_date'] >= args.date]
            # check table not empty
            check_dataframe_status(genome_trackr_table)

        if args.genome_trackr_col:
            # extract only rows with column equaling value of interest
            ## TO DO: Make this more robust
            logging.info('Subsetting GenomeTrackr table on column %s with value of %s ...' % (args.genome_trackr_col, args.genome_trackr_col_value))
            genome_trackr_table = geome_trackr_table[genome_trackr_table[args.genome_trackr_col] == args.genome_trackr_col_value]
            # check table not empty
            check_dataframe_status(genome_trackr_table)


        # get a list of all the biosample accessions for the entires of interest
        genome_trackr_biosample_accessions = list(genome_trackr_table['biosample_acc'])
        logging.info('Have a list of %s Biosamples for download from GenomeTrackr.' % str(len(genome_trackr_biosample_accessions)))

        # add to the list of sra objects for each biosample
        sra_runs += sra_runs_from_biosample_accessions(genome_trackr_biosample_accessions)

    # Use asyncio to download reads in parallel.
    loop = asyncio.get_event_loop()
    async_futures = asyncio.gather(*[x.download_task(args.download_jobs) for x in sra_runs])
    logging.info('Downloading reads...')
    loop.run_until_complete(async_futures)
    loop.close()

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
