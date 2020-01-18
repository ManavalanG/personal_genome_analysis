# Compare Dante WGS and 23andme genotyping data

Two major goals for this project:

* Perform some QC analysis on WGS data from Dante labs.
* Compare 23andme genotyping data to WGS vcf to ensure they actually match.

## What does it do?

* Perform QC analyses on WGS data
  * [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [qualimap](http://qualimap.bioinfo.cipf.es/)
* Convert 23andme genotype raw data to vcf format
* Compare Dante WGS vcf to 23andme vcf using [rtg vcfeval](https://cdn.rawgit.com/RealTimeGenomics/rtg-tools/master/installer/resources/tools/RTGOperationsManual/rtg_command_reference.html#vcfeval)

## Requirements

* OS
  * macOS (Developed and tested)
  * Linux (Not tested but should work)
  * Windows (Not tested; coin toss)

* Tools
  * [Anaconda](https://www.anaconda.com/) (Miniconda might work as well)
  * [Snakemake](https://snakemake.readthedocs.io)
    * `conda install -c bioconda -c conda-forge snakemake`

## How to run?

* Clone or download this repository
  * Clone: `git clone https://github.com/ManavalanG/personal_genome_analysis`
  * Download: See [here](https://stackoverflow.com/a/6466993/3998252) for how to.
* cd into repo directory.
  * `cd /path/to/personal_genome_analysis`
* Set up config file `configs/configs.yaml`.
  * Rename file [`configs/example_configs.yaml`](configs/example_configs.yaml) to `configs/configs.yaml`
  * Edit this file to specify file paths and sample names
* Run the snakemake pipeline
  * `snakemake -k -p --use-conda`


*Note: Feedbacks are welcome!!*