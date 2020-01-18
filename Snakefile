from pathlib import Path


################  Get Configs  ################
configfile: "configs/configs.yaml"

sample_generic = config['sample']['generic']
sample_dante = config['sample']['dante']
sample_dict = {sample_generic: sample_dante}

fastq_path = config['fastq_path']
bam_path = config['bam_path']

################################################


EXTERNAL_DIR = Path("data/external")
INTERIM_DIR = Path("data/interim")
PROCESSED_DIR = Path("data/processed")


READ_TYPES = ["R1", "R2"]

# rtg vcfeval params to be used based on filter status
VCFEVAL_PARAMS_OUTDIR_LIST = [  'pass_only',
                                'pass_only_noGT',
                                'pass_and_fail_noGT']
VCFEVAL_PARAMS_LIST = [ '',
                        '--squash-ploidy',
                        '--all-records --squash-ploidy']
VCFEVAL_PARAMS_DICT = dict(zip(VCFEVAL_PARAMS_OUTDIR_LIST, VCFEVAL_PARAMS_LIST))


##### Wildcard constraints #####
wildcard_constraints:
    sample= "|".join([sample_dante, sample_generic]),
    read="|".join(READ_TYPES)


################################################
include: "rules/common.smk"
include: "rules/qc.smk"
include: "rules/compare_vcfs.smk"
################################################


rule all:
    input:
        expand(str(PROCESSED_DIR / "qc/fastqc/{sample}-{read}.html"),
                    sample=sample_dante, read=READ_TYPES),
        expand(str(PROCESSED_DIR / "qc/qualimap/{sample}/qualimapReport.html"),
                    sample=sample_dante),
        expand(str(PROCESSED_DIR / "23andme/vcf/{sample}.vcf.gz"),
                    sample=sample_generic),
        expand(str(PROCESSED_DIR / "rtg_vcfeval_results/{sample}/{comparison_type}/done"),
                    sample=sample_generic,
                    comparison_type=VCFEVAL_PARAMS_DICT.keys())



