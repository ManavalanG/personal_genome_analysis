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

