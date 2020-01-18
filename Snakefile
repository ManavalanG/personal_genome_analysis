from pathlib import Path

configfile: "configs/configs.yaml"


sample_generic = config['sample']['generic']
sample_dante = config['sample']['dante']
sample_dict = {sample_generic: sample_dante}

fastq_path = config['fastq_path']
bam_path = config['bam_path']



EXTERNAL_DIR = Path("data/external")
INTERIM_DIR = Path("data/interim")
PROCESSED_DIR = Path("data/processed")
read_types = ["R1", "R2"]


# params to be used based on filter status
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
    read="|".join(read_types)


rule all:
    input:
        expand(str(PROCESSED_DIR / "qc/fastqc/{sample}-{read}.html"),
                    sample=sample_dante, read=read_types),
        expand(str(PROCESSED_DIR / "qc/qualimap/{sample}/qualimapReport.html"),
                    sample=sample_dante),
        expand(str(PROCESSED_DIR / "23andme/vcf/{sample}.vcf.gz"),
                    sample=sample_generic),
        expand(str(PROCESSED_DIR / "rtg_vcfeval_results/{sample}/{comparison_type}/done"),
                    sample=sample_generic,
                    comparison_type=VCFEVAL_PARAMS_DICT.keys())



rule fastqc:
    input:
        lambda wildcards: str(list(Path(fastq_path).glob(f'*{wildcards.sample}*{wildcards.read}*.fastq.gz'))[0])
    output:
        html= PROCESSED_DIR / "qc/fastqc/{sample}-{read}.html",
        zip= PROCESSED_DIR / "qc/fastqc/{sample}-{read}.zip"
    message:
        "Run fastqc QC analysis on fastq file. Sample: {wildcards.sample}, Read: {wildcards.read}"
    wrapper:
        "0.27.1/bio/fastqc"
        

rule qualimap:
    input:
        Path(bam_path) / "{sample}.bam"
    output:
        PROCESSED_DIR / "qc/qualimap/{sample}/qualimapReport.html",
    message:
        "Run qualimap QC analysis on bam file. Sample: {wildcards.sample}"
    conda:
        "configs/envs/qualimap.yaml"
    params:
        outdir = lambda wildcards, output: Path(output[0]).parent,
        java_mem_size = "8G"
    shell:
        r"""
        qualimap bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            -outdir {params.outdir} \
            --java-mem-size={params.java_mem_size} 
        """


rule download_ref_genome:
    output:
        ref = config['ref_genome'],
    message:
        "Downloads reference genome. Note: Bcftools needs it in bgzip format."
    params:
        url = config['ref_genome_url']
    shell:
        r"""
        curl -o {output.ref}.gz "{params.url}" 
        gzip -d --quiet {output.ref}.gz || true  # ignore trailing garbage
        """


rule convert2vcf_23andme:
    input:
        tsv = config['23andme_file'],
        ref = config['ref_genome']
    output:
        PROCESSED_DIR / "23andme/vcf/{sample}.vcf.gz"
    message:
        "Convert 23andme tsv file to vcf format. Sample: {wildcards.sample}"
    conda:
        "configs/envs/bcftools.yaml"
    params:
        sample_name=sample_generic
    shell:
        r"""
        bcftools convert \
            --tsv2vcf {input.tsv} \
            --fasta-ref {input.ref} \
            --samples {params.sample_name} \
            --output-type z \
            --output {output}
        """


rule tabix_index:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    message:
        "Tabix indexing vcf"
    conda:
        "configs/envs/htslib.yaml"
    shell:
        """
        tabix -p vcf {input}
        """


rule ref_sdk:
    input:
        ref = config['ref_genome']
    output:
        directory(config['ref_genome'] + '.sdf')
    message:
        "Creates RTG Sequence Data File (SDF) for ref genome"
    conda:
        "configs/envs/rtg.yaml"
    shell:
        """
        rtg format -o {output} {input}
        """


rule vcfeval_compare:
    input:
        dante = lambda wildcards: str(Path(config['vcf_path']) / f"{sample_dict[wildcards.sample]}.filtered.snp.vcf.gz"),
        dante_index = lambda wildcards: str(Path(config['vcf_path']) / f"{sample_dict[wildcards.sample]}.filtered.snp.vcf.gz.tbi"),
        vcf_23andme = PROCESSED_DIR / "23andme/vcf/{sample}.vcf.gz",
        vcf_23andme_index = PROCESSED_DIR / "23andme/vcf/{sample}.vcf.gz.tbi",
        ref = config['ref_genome'] + '.sdf'
    output:
        weighted = str(PROCESSED_DIR / "rtg_vcfeval_results/{sample}/{comparison_type}/weighted_roc.tsv.gz"),
        fn = str(PROCESSED_DIR / "rtg_vcfeval_results/{sample}/{comparison_type}/fn.vcf.gz"),
        tp_baseline = str(PROCESSED_DIR / "rtg_vcfeval_results/{sample}/{comparison_type}/tp-baseline.vcf.gz"),
        tp = str(PROCESSED_DIR / "rtg_vcfeval_results/{sample}/{comparison_type}/tp.vcf.gz"),
        fp = str(PROCESSED_DIR / "rtg_vcfeval_results/{sample}/{comparison_type}/fp.vcf.gz"),
        done = str(PROCESSED_DIR / "rtg_vcfeval_results/{sample}/{comparison_type}/done"),
    conda:
        "configs/envs/rtg.yaml"
    params:
        extra_params = lambda wildcards, output: VCFEVAL_PARAMS_DICT[wildcards.comparison_type],
        outpath = lambda wildcards, output: str(Path(output.weighted).parent)
    message:
        "RTG vcfeval comparing '{wildcards.comparison_type}' vcf sample: {wildcards.sample}, Additional parameters used, if any: {params.extra_params}"
    shell:
        r"""
        rm -rf {params.outpath}  # RTG-vcfeval requires dir not to exist or it fails.
        rtg vcfeval \
            --baseline {input.dante} \
            --calls {input.vcf_23andme} \
            --output {params.outpath} \
            --template {input.ref} \
            --vcf-score-field QUAL  {params.extra_params}
        """
