# Compare Dante WGS and 23andme genotyping data

## What it does?

* Perform QC analyses on WGS data
  * [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [qualimap](http://qualimap.bioinfo.cipf.es/)
* Convert 23andme genotype raw data to vcf format
* Compare Dante WGS vcf to 23andme vcf using [rtg vcfeval](https://cdn.rawgit.com/RealTimeGenomics/rtg-tools/master/installer/resources/tools/RTGOperationsManual/rtg_command_reference.html#vcfeval)

## Requirements

* Supported OS:
  * macOS (Developed and tested)
  * Linux (Not tested but should work)
  * Windows (Not tested; coin toss)

* Tools
  * [Anaconda](https://www.anaconda.com/) (Miniconda might work as well)
  * [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda)
    * `conda install -c bioconda -c conda-forge snakemake`

## How to run?

* Clone or download this repository
  * Clone: `git clone xxx`
  * Download: just download this re
* cd into repo directory.
  * `cd /path/to/personal_genome_analysis`
* Set up config file as described below.
* Run the snakemake pipeline
  * `snakemake -k -p --use-conda`


## Set up `configs/configs.yaml`