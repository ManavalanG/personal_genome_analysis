rule convert2vcf_23andme:
    input:
        tsv = config['23andme_file'],
        ref = config['ref_genome']
    output:
        PROCESSED_DIR / "23andme/vcf/{sample}.vcf.gz"
    message:
        "Convert 23andme tsv file to vcf format. Sample: {wildcards.sample}"
    conda:
        "../configs/envs/bcftools.yaml"
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


rule vcfeval_compare:
    input:
        dante = lambda wildcards: str(wgs_vcf_path / f"{sample_dict[wildcards.sample]}.filtered.snp.vcf.gz"),
        dante_index = lambda wildcards: str(wgs_vcf_path / f"{sample_dict[wildcards.sample]}.filtered.snp.vcf.gz.tbi"),
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
        "../configs/envs/rtg.yaml"
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
