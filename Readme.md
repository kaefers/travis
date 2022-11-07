# TRAVIS - TRAnscriptome VIrus Scanner

The main purpose for TRAVIS is to scan samples directed towards a certain virus group.
TRAVIS will read user provided libraries and configuration files in comma separated value
format (CSV). It is possible to add all valuable information that can be important for
downstream analyses tailored to the project into the libraries. This information could
be taxonomy, host-ranges, symptoms etc. The specified reference sequences will be
downloaded from NCBI and then systematic searches at amino acid level will be performed.
The output will contain multiple text files, tables and visualizations.
The user has to specify a ’main’ gene that all of the viruses in the library have in common
and is more or less conserved. This gene could be e.g. a polymerase or a capsid gene. All
other genes will be regarded as ’company’ genes. By default, TRAVIS will search first for
the ’main’ genes first and only starts searching for ’company’ genes in samples that are
positive for the ’main’ gene. It acts as a marker gene and running ’company’ searches only
on ’main’ positive samples will reduce the whole search time. After identification of the
potential viral sequences (later called 'suspicious sequences'), they are checked against
a larger, comprehensive database to spot false positives more easily. I suggest to use the
non-redundant protein database from NCBI.

TRAVIS is separated into three parts that are executed subsequently: TRAVIS Henchman,
TRAVIS Core, and TRAVIS Scavenger. 

This guide is not comprehensive for all functionality as more features will be implemented
in the future. It assumes that you have basic knowlege about the use of the Unix command
line.

TRAVIS is still in beta and probably contains bugs. Feel free to use it at own risk
and please report any bugs to me (simonkaefer@outlook.com). 

## Concept and Workflow


Each of the TRAVIS subprograms (Henchman, Core and Scavenger) is called with a single
manifest file (TRAVIS Control Center, TCC). If you want to have a completely automated
run of TRAVIS without manual interaction, you can call them subsequently in e.g. a bash
script.
```
$perl TRAVIS_Henchman_pt1.pl TCC.csv
$perl TRAVIS_Henchman_pt2.pl TCC.csv
$perl TRAVIS_Core. pl TCC.csv
$perl TRAVIS_Scavenger.pl TCC.csv
```


The configuration of TRAVIS is done by creating manifest files that contain all necessary
information. These manifest files are actually plain-text comma separated value files
(CSV) that you can edit either in a text editor of your choice or spreadsheet software
like LibreOffice, OpenOffice or MS Excel. But when you export the CSVs make sure that
the export has been done properly. That means opening it in a text editor and check
whether the entries are actually separated by comma and not by semi-colon (the german
version of Excel does that!). Another problem are quotes around the entries. TRAVIS does
not like quotes. Also please only use alphanumeric characters, dashes and underscores for
whatever you enter in the manifest files. Note that TRAVIS internally uses double and
triple underscores as separators.
However, TRAVIS Henchman creates a manifest file called ’Troubling TRAVIS Table’ (TTT), where all intended 
searches for TRAVIS Core are listed. You can adjust
this table according to your specific needs. This can drastically reduce calculation time. It
is also possible to manually create a TTT or use an old one and skip TRAVIS Henchman completely.
Henchman was split into two parts. The first part downloads and sorts your references, and
the second part prepares the data for the searches. Depending on your reference library, 
this can take quite some time and computing resources. 
As some high performance clusters have ridiculously strict firewalls, you might want to 
download the references on a standard computer with a proper internet connection and 
let the cluster do the heavy lifting afterwards. 
Currently, TRAVIS uses [curl](https://curl.haxx.se/) to download files from NCBI. If the 
access is not restricted on your cluster, you should be able to run everyting on it.
I recently added [wget](https://www.gnu.org/software/wget/) support as well as I ran into trouble with the institute's firewall.
It's also the new default. 
If you have a large dataset, you can run TRAVIS Scavenger on the same TCC while
TRAVIS Core is still running in order to get the results that have already been generated.
Because TRAVIS runs all intended searches completely for each sample and logs the results,
it is also possible to resume calculations from the last processed sample.
Each sample gets an own set of output files based on the given ID in the sample library.
These output files will be fastas, tables (CSV) and visualizations of (SVG) the matches.
I recommend to open the SVGs in a web-browser because detailed information about the
matches will be displayed when you hover your cursor over certain elements. This has
been tested with Mozilla Firefox and Google Chrome under Windows 10, 11 and Ubuntu 16.04, 20.04.
These details can also be found in the corresponding CSV.







## Installation
TRAVIS is written in perl and should work out of the box on most UNIX systems.
If you have compiled versions of HMMER3, BLAST+, DIAMOND, MMSeqs2 and MAFFT, you
are good to go. You can specify the paths in the configuration file. However, if
you have the programs installed and working with shortcuts/aliases, you can also use
these. A combination of HMMER3 (v. 3.3,1), BLAST+ (v. 2.12.0), DIAMOND (v. v2.0.15.153) MMseqs2 (v.
45111b641859ed0ddd875b94d6fd1aef1a675b7e) and MAFFT (v. 7.490) worked well on
Ubuntu 20.04 but i guess, other versions won’t make problems as long the respective
developers do not change their parameter calls or output format...

On Ubuntu 20.04, installation of TRAVIS works as follows:
```
sudo apt install git curl wget
git clone https://github.com/kaefers/travis
```


Please use the respective documentation of the dependencies for installation.
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [HMMER3](http://hmmer.org/)
* [MMSeqs2](https://github.com/soedinglab/MMseqs2)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/linux.html)
* [DIAMOND](https://github.com/bbuchfink/diamond) (please check https://github.com/bbuchfink/diamond/wiki/2.-Installation for installation with BLAST database support)

or use the singularity image that can be found [here](https://doi.org/10.6084/m9.figshare.21510204)

You need also need a local copy of a general database for the false-positive check. I suggest the non-redundant protein database of NCBI. However, it
has gotten a little too big to be useful with blastp over the last few years, unless you have a massive amount of cores available. Good news is, that it is
easy to use a BLAST-formatted database with DIAMOND and DIAMOND is reasonably fast even on slower machines while getting most of the same results as BLAST.
There is a nice little tool for the download of NCBI databases provided by NCBI: [update_blastdb.pl](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl).


```
perl update_blastdb.pl --decompress nr

```
This will take a while and use about 350GB of hard drive space. This database can be indexed for the use with DIAMOND:

```
diamond prepdb --db /path/to/nr
```

If you think you can afford skipping this step, you can skip it completely by setting the 'blastp_db_full' and 'diamond_db_full' parameter in TCC to 'skip'.



## Tutorial

Here, i will walk you through running TRAVIS on the test data that is contained within the package.
We will start from the main directory, where you downloaded TRAVIS, which is /home/simon/travis/ for me.

```
cd /home/simon/travis/test_data/
unzip testfastas.zip
```

Then edit the 'reo_test_TCC.csv' to match your system, settings and paths. For all the options, see the 'Options' section.
I have been using pre-compiled versions of the dependencies and specified the
absoulte path to them. This seemed to work well on all occasions so far.


```
cd ..
perl TRAVIS_Henchman_pt1 test_data/reo_test_TCC.csv
```

This should have downloaded and sorted all reference data you specified in the reference library.
To cluster, split and align them, you need to run pt 2:

```
perl TRAVIS_Henchman_pt2 test_data/reo_test_TCC.csv
```

After that you can check the intended searches in the Troubling TRAVIS Table (TTT)
You can e.g.:
* switch off searches
* check and modify alignments that are the basis for hmmsearch
* change combination of search tools on certain groups/clusters
* add other proteins/groups/clusters to the ’main’ pool





## Options
These are the options that can be set in the TCC. Have a look at reo_test_TCC.csv in test_data/ for examples.

#### database_name
Contains the name of the database to be generated by TRAVIS Henchman. If you already
have a prepared database that you want to use as it is, you can skip TRAVIS Henchman,
modify TCC, and start TRAVIS Core.


#### resume_calculation
In case of crashes, you can resume the calculation based on the last save point. That save
point is the last completely searched sample. This is encoded by
1 or 0 encoding on/off

#### sample_dir
Specifies the path to the nucleotide data.


#### ORF_dir
Specifies the path to where the ORF data should be stored.


#### ORF lengths limits
Sets limits to the ORFs to be extracted in number of amino acids (AA).
I suggest to not use anything below 20 AA.

```
min_ORF_length
max_ORF_length
```

#### sample_library
Specifies the path to the sample library with ’filename’,’ID’,’factor1’,’factor2’,’factor3’...
You can add as many factors as you want depending on the information you need to be
associated with the results later on.

Required:
* 1 st column has to be the filename of the sample
* 2 nd column has to be a unique name or ID
* any number of columns containing any information

#### reference_library
Specifies the path to the reference library. As TRAVIS is relying on NCBI up to now, it is
necessary to specify either accession numbers (separated by ’&’) in the reference library or
the path to the respective assembly report on the NCBI-FTP server, if available. However,
you as well need to specify the NT accession number and the PID of the ’main’ gene. you
can add as many factors as you want. these can be used for naming the references and
sorting them into subgroups

Required:
* 1 st column has to be an acronym or ID
* 2 nd column has to be a unique name
* any number of columns containing any information
* the last three columns have to be ’all_NT_ACC,<main>_NT_ACC,<main>_PID’
where <main> can be replaced by a meaningful name

#### Local Reference Databases
Specifies the path to local reference databases. This database saves everything related to
the reference library so you do not have to download everything from NCBI over and over
again.

```
local_reference_database
reference_fastas
reference_gbx
```


#### header_names
Names of the columns that you want to drag through the whole analysis included in the
header of the reference sequences. They have to be identical to the column names in your
reference library.


#### split_references
Names of the columns that you want to split the references by. So you can e.g. create
subgroups by genus or family.

#### sample_subset
If this is set to ’main_positive’, only company sequences will be searched if main sequences
were found in the respective sample. This can still be changed in TTT before running
TRAVIS Core. Can be set to 'main_positive' or 'all'.

#### result_dir
All relevant results will be stored here if not declared otherwise.

#### TTT
Specifies the path to the Troubling TRAVIS Table.

#### nCPU
Specifies how many processors can be used.


#### max_references
Limits how many references will be plotted per matching ORF in the Scavenger output.


#### HMMER3
Specifies paths and settings of HMMER3.

```
hmmbuild
hmmsearch
hmmsearch_settings
jackhmmer
jackhmmer_settings
```

#### MAFFT
Specifies paths and settings of MAFFT. In my experience, if you want to use a portable
version of MAFFT, the proper $PATHs have to be configured. By specifying the location
of the MAFFT_BINARIES, i could easily solve issues regarding that.

```
mafft
mafft_settings
mafft_binaries
```

#### MMSeqs2
Specifies paths and settings of MMSeqs2. ’minimal_cluster_size’ is for the clustering of
the company sequences by TRAVIS Henchman.

```
mmseqs
mmseqs_cluster_settings
mmseqs_search_settings
minimal_cluster_size
```

#### BLASTP
Specifies paths and settings of BLASTP.
The 'blastp_db_full' parameter expects the path to your local 'full' protein database (i.e. the nr). Skipping this step can be achieved with using 'skip'. If you want to cheat with another database, you can use this here, too. 
It just needs to be a blast database that has been created with makeblastdb.

```
blastp
blastp_settings
makeblastdb
blastp_db_full
```

#### DIAMOND
Specifies paths and settings of DIAMOND.
The 'diamond_db_full' parameter works like the 'blastp_db_full' parameter. Please remember to use 'diamond prepdb' before running the
search.

```
diamond
diamond_settings
diamond_db_full
```
