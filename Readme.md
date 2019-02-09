# TRAVIS - TRAnscriptome VIrus Scanner

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

This guide is not comprehensive for all functionality as more features will be implemented
in the future. It assumes that you have basic knowlege about the use of the Unix command
line.

This pipeline is still in beta and probably contains bugs. Feel free to use it at own risk
and report any bugs to me. 

## Concept and Workflow


Each of the TRAVIS subprograms (Henchman, Core and Scavenger) is called with a single
manifest file as a parameter. If you want to have a completely automated
run of TRAVIS without manual interaction, you can call them subsequently in e.g. a bash
script.
```
$perl TRAVIS_Henchman_pt1.pl TCC.csv
$perl TRAVIS_Henchman_pt2.pl TCC.csv
$perl TRAVIS_Core. pl TCC.csv
$perl TRAVIS_Scavenger.pl TCC.csv
```

However, TRAVIS Henchman creates a manifest file called ’Troubling TRAVIS Table’, where all intended 
searches for TRAVIS Core are listed. You can adjust
this table according to your specific needs. This can drastically reduce calculation time. It
is also possible to manually create a TTT or use an old one and skip TRAVIS Henchman.
If you have a large dataset, you can run TRAVIS Scavenger on the same TCC while
TRAVIS Core is still running in order to get the results that have already been generated.
Because TRAVIS runs all intended searches completely for each sample and logs the results,
it is also possible to resume calculations from the last processed sample.
Each sample gets an own set of output files based on the given ID in the sample library.
These output files will be fastas, tables (CSV) and visualizations of (SVG) the matches.
I recommend to open the SVGs in a web-browser because detailed information about the
matches will be displayed when you hover your cursor over certain elements. This has
been tested with Mozilla Firefox and Google Chrome under Windows 10 and Ubuntu 16.04.
These details can also be found in the corresponding CSV










copy via 

```
git clone dingens
```

install requirered programs

[//]: #* [HMMer3](http://www.dropwizard.io/1.0.2/docs/) - hammer time!
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [HMMER3](http://hmmer.org/)
* [MMSeqs2](https://github.com/soedinglab/MMseqs2)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/linux.html)

or download package from url: blalbalba.dings.com/file5213789


### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo