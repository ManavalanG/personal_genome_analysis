from pathlib import Path

configfile: "configs/configs.yaml"

sample = config['sample']
fastq_path = config['fastq_path']
bam_path = config['bam_path']



EXTERNAL_DIR = Path("data/external")
INTERIM_DIR = Path("data/interim")
PROCESSED_DIR = Path("data/processed")
read_types = ["R1", "R2"]

##### Wildcard constraints #####
wildcard_constraints:
    sample=sample,
    read="|".join(read_types)


rule all:
    input:
        # expand(str(PROCESSED_DIR / "qc/fastqc/{sample}-{read}.html"),
        #             sample=sample, read=read_types),
        # expand(str(PROCESSED_DIR / "qc/qualimap/{sample}/qualimapReport.html"),
        #             sample=sample),
        EXTERNAL_DIR / "ref_genome/hs37d5.fa.gz"

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
        "Run qualimap QC analysis on bam file. Sample: {sample}"
    params:
        outdir = lambda wildcards, output: Path(output[0]).parent,
        java_mem_size = "10G"
    shell:
        r"""
        qualimap bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            -outdir {params.outdir}
            --java-mem-size={params.java_mem_size} \
        """


rule download_ref_genome:
    output:
        EXTERNAL_DIR / "ref_genome/hs37d5.fa.gz"
    message:
        "Downloads reference genome"
    params:
        url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    shell:
        r"""
        curl -o {output} \
            {params.url}
        """