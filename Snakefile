from pathlib import Path

configfile: "configs/configs.yaml"

sample = config['sample']
fastq_path = config['fastq_path']


INTERIM_DIR = Path("data/interim")
PROCESSED_DIR = Path("data/processed")


rule all:
    input:
        expand(str(PROCESSED_DIR / "qc/fastqc/{sample}.html"),
                    sample=sample)


rule fastqc:
    input:
        r1 = Path(fastq_path).glob('*R1*.fastq.gz'),
        r2 = Path(fastq_path).glob('*R2*.fastq.gz')
    output:
        html= PROCESSED_DIR / "qc/fastqc/{sample}.html",
        zip= PROCESSED_DIR / "qc/fastqc/{sample}zip"
    wrapper:
        "0.27.1/bio/fastqc"
        

# rule qualimap:
#     input:
#         Path(fastq_path).glob('*R1*.fastq.gz')
#     output:
#     # message:
#     shell: