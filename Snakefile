from pathlib import Path

configfile: "configs/configs.yaml"

sample = config['sample']
fastq_path = config['fastq_path']
bam_path = config['bam_path']



INTERIM_DIR = Path("data/interim")
PROCESSED_DIR = Path("data/processed")


rule all:
    input:
        expand(str(PROCESSED_DIR / "qc/fastqc/{sample}.html"),
                    sample=sample),
        expand(str(PROCESSED_DIR / "qc/qualimap/{sample}/qualimapReport.html"),
                    sample=sample),

rule fastqc:
    input:
        r1 = Path(fastq_path).glob('*R1*.fastq.gz'),
        r2 = Path(fastq_path).glob('*R2*.fastq.gz')
    output:
        html= PROCESSED_DIR / "qc/fastqc/{sample}.html",
        zip= PROCESSED_DIR / "qc/fastqc/{sample}zip"
    message:
        "Run fastqc QC analysis on fastq file"
    wrapper:
        "0.27.1/bio/fastqc"
        

rule qualimap:
    input:
        Path(bam_path).glob('*.bam')
    output:
        PROCESSED_DIR / "qc/qualimap/{sample}/qualimapReport.html",
    message:
        "Run qualimap QC analysis on bam file"
    params:
        outdir = lambda wildcards, output: Path(output).parent,
        java_mem_size = "10G"
    shell:
        """
        qualimap bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            -outdir {params.outdir}
            --java-mem-size={params.java_mem_size} \
        """