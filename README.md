# SRA Read Downloader

This will download reads from the NCBI SRA database. Its input are accessions of either the SRA, ENA or DDBJ format. Accessions can be either BioProject, BioSample, Experiments or Run.

Briefly, BioProjects are accession numbers linked to projects (prefix PRJ or SRP/ERP/DRP), which contain one or many BioSamples (prefix SAMN or SRS/ERS/DRS). A BioSample is a sample of an organism, and within a BioSample there can be multiple Experiments (SRX/ERX/DRX) and Runs (SRR/ERR/DRR). An Experiment is an experimental setup for sequencing. Each Run is a sequencing run performed on some sequencing platform, using a specific experimental setup specified in the Experiment. For more information and perhaps a better explanation of the different accession types, see here .

Ultimately, we are interested in the Runs, as these are our read files. All read files downloaded by this script will be named with the Run accession, regardless of whether you provide BioProject, BioSample or Experiment accessions. An output file (detailed below) will provide you with the necessary information to link your resulting Run accessions back to whichever accession you provided.

As input, this script can take any combination of the following options:

a text file of accessions (one per line) with --accesion_list
a list of accessions to the --accession argument, separated by spaces
a GenomeTrackr species
The script will then take the accession and identify the read files associated with this accession. In cases where there are multiple read files associated with a BioSample, the script will take only the most recent set of read files. If there are both short and long reads associated with the BioSample, it will download the most recent set of long read files and the most recent set of short read files.

To download reads from GenomeTrackr, the --species flag must be used with one of the following species. This option must be specified exactly, as these are the names of the folders in the GenomeTrackr database.

Acinetobacter
Campylobacter
Citrobacter_freundii
Elizabethkingia
Enterobacter
Escherichia_coli_Shigella
Klebsiella
Klebsiella_oxytoca
Kluyvera_intermedia
Legionella_pneumophila
Listeria
Morganella
Mycobacterium_tuberculosis
Neisseria
Providencia
Pseudomonas_aeruginosa
Salmonella
Serratia
Staphylococcus_pseudintermedius
Vibrio_parahaemolyticus
When the script is executed with the --species flag, the latest GenomeTrackr metadata table is downloaded from the GenomeTrackr project for that species. This metadata table can then be filtered by the script so only a subset of these reads are downloaded (using the options below).

 If only --species is provided, all reads from that species will be downloaded (which will take a VERY long time!)
 The --date flag can be used to download only reads which have been uploaded since the specified date (format must be YYYY-MM-DD).
 The --genome_trackr_col flag can be used to specify a particular column in the GenomeTrackr metadata table, and with the --genome_trackr_col_value flag, can filter the metadata table by only including reads which contain the 
--genome_trackr_col_value in the --genome_trackr_col.
Both the subsetted metadata table and original metadata table are left behind for the user to inspect.

Reads are downloaded in parallel using fastq-dump from SRAToolKit. Eight read sets are downloaded in parallel at once by default - this can be changed using the --download_jobs argument.  

Finally, the name of the log file can be altered using the --log flag (default is download_reads.log)

After all reads have been downloaded, a master_list.csv file will be left behind. This file details each read set which has been download successfully, as well as additional metadata (including Run, Experiment, BioSample and BioProject accessions and data type [Illumina/Oxford Nanopore etc]).

Any reads which could not be downloaded are listed in the failed_accessions.txt file, with their reasons for failure.

This script requires PYTHON 3. It should be run in screen, as it does NOT submit jobs to the slurm queue, so can be run on any computer/cluster. Reads will be downloaded into the folder where the script is run.




### Brief usage

Download reads with accessions in a text file:<br>
```
python /vlsci/SG0006/shared/code/holtlab/download_reads.py --accession_list accessions.txt
```

Download reads using a file and also some other accessions:<br>
```
python   /vlsci/SG0006/shared/code/holtlab/download_reads.py --accession_list accessions.txt --accessions SAMN06112191 SRR1023455
```

Download reads from GenomeTrackr:<br>
```
python   /vlsci/SG0006/shared/code/holtlab/download_reads.py --species Escherichia_coli_Shigella --date 2017-12-26
```

Use all options!<br>
```
python /vlsci/SG0006/shared/code/holtlab/download_reads.py --accession_list accessions.txt --accessions SAMN06112191 SRR1023455 --species Escherichia_coli_Shigella --date 2017-12-26
```




## Full usage

```
usage: sra_read_downloader.py [-h] [--accession_list ACCESSION_LIST]

                         [--accessions ACCESSIONS [ACCESSIONS ...]]

                         [--logfile LOGFILE] [--download_jobs DOWNLOAD_JOBS]

                         [--species SPECIES] [--date DATE]

                         [--genome_trackr_col GENOME_TRACKR_COL]

                         [--genome_trackr_col_value GENOME_TRACKR_COL_VALUE]



Download reads from NCBI



optional arguments:

  -h, --help            show this help message and exit

  --accession_list ACCESSION_LIST

                        File of accessions (one per line)

  --accessions ACCESSIONS [ACCESSIONS ...]

                        NCBI accessions (BioProject, BioSample or SRA,

                        separated by spaces)

  --logfile LOGFILE     Log file

  --download_jobs DOWNLOAD_JOBS

                        Number of simultaneous downloads to run at once

  --species SPECIES     Species from GenomeTrackr. Must be one of the

                        following: Acinetobacter, Campylobacter,

                        Citrobacter_freundii, Elizabethkingia, Enterobacter,

                        Escherichia_coli_Shigella, Klebsiella,

                        Klebsiella_oxytoca, Kluyvera_intermedia,

                        Legionella_pneumophila, Listeria, Morganella,

                        Mycobacterium_tuberculosis, Neisseria, Providencia,

                        Pseudomonas_aeruginosa, Salmonella, Serratia,

                        Staphylococcus_pseudintermedius,

                        Vibrio_parahaemolyticus

  --date DATE           Only required when downloading reads from

                        GenomeTrackr. Will only download reads uploaded on or

                        after the date specified in this argument. The

                        corresponding column in the GenomeTrackr table is

                        "target_creation_date". Date MUST be in the following

                        format: YYYY-MM-DD.

  --genome_trackr_col GENOME_TRACKR_COL

                        Name of column in GenomeTrackr table which will be

                        used to select only rows which equal a particular

                        value - this value can be set with

                        --genome_trackr_col_value.

  --genome_trackr_col_value GENOME_TRACKR_COL_VALUE

                        Value in column of GenomeTrackr to use to select

                        specific rows. Column name can be set with

                        --genome_trackr_col.
```




## License

GNU General Public License, version 3
