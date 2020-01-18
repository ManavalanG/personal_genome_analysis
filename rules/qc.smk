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
