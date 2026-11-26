from pathlib import Path
from rispy_cli.io import build_fq_by_sample
from rispy_cli.cli_utils import deepvariant_shell_activation

DV_ACTIVATE = deepvariant_shell_activation()

root = Path(config["paths"]["output_dir"])

fastq = Path(config["paths"]["fastq"])

ref = config["paths"]["reference"]

fq_by_sample = build_fq_by_sample(fastq=fastq)

rule all:
    input:
        NanoStats_UNIQUE = expand (str(root / '{sample}' / 'nanoplot_after_filtering' / 'NanoStats.txt'), sample = fq_by_sample.keys()),
        NanoStats = expand (str(root / '{sample}' / 'nanoplot_raw' / 'NanoStats.txt'), sample = fq_by_sample.keys()),
        Filtered = expand (str(root / '{sample}' / 'deep_variant' / 'Filtering' / 'filtered_variants.vcf.gz'), sample = fq_by_sample.keys())

rule nanoplot_raw:
    """
    Quality control of raw reads 
    """
    input:
        FQ = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
    output:
        NanoStats = root / '{sample}' / 'nanoplot_raw' / 'NanoStats.txt',
    params:
        NanoStats_dir = lambda wildcards: root / wildcards.sample / 'nanoplot_raw',
    threads: config["resources"]["nanoplot_raw"]["threads"]
    resources:
        mem_mb = config["resources"]["nanoplot_raw"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'nanoplot_raw.txt'
    shell:
        """
        NanoPlot --fastq {input.FQ} -t {threads} -o {params.NanoStats_dir}
        """

rule extract_headers:
    """
    Extract the read headers of simplex and duplex reads
    """
    input:
        FQ = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
    output:
        TSV_ALL = root / '{sample}' / 'read-ids' / 'all.tsv', #tsv of headers of all reads
        TSV_ids = root / '{sample}' / 'read-ids' / 'idstofilter.tsv', #tsv of read identifiers of simplex reads used to generate duplex reads
        TSV_filter = root / '{sample}' / 'read-ids' / 'filter.tsv' #tsv of full headers of simplex reads used to generate duplex reads
    threads: config["resources"]["extract_headers"]["threads"]
    resources:
        mem_mb = config["resources"]["extract_headers"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'extract_headers.txt'
    shell:
        """
        seqkit fx2tab -n {input.FQ} > {output.TSV_ALL}
        seqkit fx2tab -n {input.FQ} | cut -f 1 | grep ';' | tr ';' '\n' > {output.TSV_ids} || true
        cat {output.TSV_ALL} | grep -f {output.TSV_ids} | grep -v ';' > {output.TSV_filter} || true
        """

rule filter_reads:
    """
    Filter simplex reads that were used to generate duplex reads from fastq 
    """
    input:
        FQ = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
        TSV = root / '{sample}' / 'read-ids' / 'filter.tsv'
    output:
        UNIQUE = root / '{sample}' / 'unique_fastq' / 'unique.fastq.gz',
    threads: config["resources"]["filter_reads"]["threads"]
    resources:
        mem_mb = config["resources"]["filter_reads"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'filter_reads.txt'
    shell:
        """
        seqkit grep -n -v -f {input.TSV} {input.FQ} -o {output.UNIQUE} || true
        """

rule trim_reads:
    """
    Trimming of ONT reads
    """
    input:
        FQ = root / '{sample}' / 'unique_fastq' / 'unique.fastq.gz',
    output:
        ONT_FILT = root / '{sample}' / 'chopper' / 'unique_filter.fastq.gz',
    params:
        q = config["parameters"]["chopper"]["quality"],
        l = config["parameters"]["chopper"]["min_len"]
    threads: config["resources"]["trim_reads"]["threads"]
    resources:
        mem_mb = config["resources"]["trim_reads"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'trim_reads.txt'
    shell:
        """
        gunzip -c {input.FQ} | chopper -q {params.q} -l {params.l} --threads {threads} | gzip > {output.ONT_FILT}
        """

rule nanoplot_after_filtering:
    """
    Quality control of trimmed reads
    """
    input:
        ONT_FILT = root / '{sample}' / 'chopper' / 'unique_filter.fastq.gz',
    output:
        NanoStats_UNIQUE = root / '{sample}' / 'nanoplot_after_filtering' / 'NanoStats.txt',
    params:
        NanoStats_dir = lambda wildcards: root / wildcards.sample / 'nanoplot_after_filtering' ,
    threads: config["resources"]["nanoplot_after_filtering"]["threads"]
    resources:
        mem_mb = config["resources"]["nanoplot_after_filtering"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'nanoplot_after_filtering.txt'
    shell:
        """
        NanoPlot --fastq {input.ONT_FILT} -t {threads} -o {params.NanoStats_dir}
        """

rule index_reference:
    """
    Indexing the reference fasta
    """
    input:
        FASTA = ref
    output:
        out_ref = root / 'reference' / 'reference.fasta',
        MMI = root / 'reference' / 'refgenomes.mmi'
    params:
        minimap2_preset = 'map-ont'
    threads: config["resources"]["index_reference"]["threads"]
    resources:
        mem_mb = config["resources"]["index_reference"]["memory"]
    benchmark: root / 'reference' / 'benchmarks' / 'index_reference.txt'
    shell:
        """
        cp {input.FASTA} {output.out_ref}
        minimap2 -I4G -t {threads} -x {params.minimap2_preset} -d {output.MMI} {output.out_ref}
        samtools faidx {output.out_ref}
        """

rule alignment:
    """
    Mapping the long reads to reference and convert to bam, then sorting and indexing the bam file
    """
    input:
        ONT_FQ = root / '{sample}' / 'chopper' / 'unique_filter.fastq.gz',
        MMI = root / 'reference' / 'refgenomes.mmi'
    output:
        BAM = root / '{sample}' / 'mapping' / 'mapped_reads_sorted.bam',
    params:
        minimap2_preset = 'map-ont',
        tmp_dir= root / 'tmp'
    threads: config["resources"]["alignment"]["threads"]
    resources:
        mem_mb = config["resources"]["alignment"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'alignment.txt'
    shell:
        """
        minimap2 -a -x {params.minimap2_preset} {input.MMI} {input.ONT_FQ} -t {threads} \
        | samtools view -b - -@ {threads} | samtools sort - -o {output.BAM} -T {params.tmp_dir} -@ {threads}
        samtools index {output.BAM}
        """

rule make_examples:
    """
    Running make_examples to generate examples for deepvariant
    """
    input:
        BAM = rules.alignment.output.BAM,
        REF = rules.index_reference.output.out_ref
    output:
        EXAMPLES = root  / '{sample}' / 'deep_variant' / 'examples.tfrecord.gz'
    params:
        CHECKPOINT = config["parameters"]["call_variants"]["CHECKPOINT"],
        ACTIVATE = DV_ACTIVATE,
        REGIONS = config["parameters"]["call_variants"]["regions"]
    threads: config["resources"]["make_examples"]["threads"]
    resources:
        mem_mb = config["resources"]["make_examples"]["memory"]
    benchmark: root / '{sample}' / 'deep_variant' / 'benchmarks' / 'make_examples.txt'
    shell:
        """
        bash -c "{params.ACTIVATE} && \
        make_examples --mode calling --call_small_model_examples=False \
        --multi_allelic_mode='include_het_alt_images' --select_variant_types='all' --regions={params.REGIONS} \
        --vsc_min_count_snps='1' --vsc_min_fraction_snps='0.001' --vsc_min_count_indels='1' --vsc_min_fraction_indels='0.001' \
        --checkpoint={params.CHECKPOINT} --ref={input.REF} --reads={input.BAM} --examples={output.EXAMPLES}"
        """

rule call_variants:
    """
    Running call_variants to call variants with the examples generated from make_examples
    """
    input:
        EXAMPLES = rules.make_examples.output.EXAMPLES
    output:
        CALL = root  / '{sample}' / 'deep_variant' / 'call_variants_output-00000-of-00001.tfrecord.gz'
    params:
        CHECKPOINT = config["parameters"]["call_variants"]["CHECKPOINT"],
        Call = lambda wildcards: root / wildcards.sample / 'deep_variant' / 'call_variants_output.tfrecord.gz',
        ACTIVATE = DV_ACTIVATE
    threads: config["resources"]["call_variants"]["threads"]
    resources:
        mem_mb = config["resources"]["call_variants"]["memory"]
    benchmark: root / '{sample}' / 'deep_variant' / 'benchmarks' / 'call_variants.txt'
    shell:
        """
        bash -c "{params.ACTIVATE} && \
        call_variants --examples={input.EXAMPLES} --checkpoint={params.CHECKPOINT} --outfile={params.Call} --writer_threads={threads}"
        """

rule postprocess_variants:
    """
    Running postprocess_variants from deepvariant
    """
    input:
        CALL = rules.call_variants.output.CALL,
        REF = rules.index_reference.output.out_ref
    output:
        VCF = root / '{sample}' / 'deep_variant' / 'variants.vcf.gz'
    params:
        Call = lambda wildcards: root / wildcards.sample / 'deep_variant' / 'call_variants_output.tfrecord.gz',
        ACTIVATE = DV_ACTIVATE,
        REGIONS = config["parameters"]["call_variants"]["regions"]
    threads: config["resources"]["postprocess_variants"]["threads"]
    resources:
        mem_mb = config["resources"]["postprocess_variants"]["memory"]
    benchmark: root / '{sample}' / 'deep_variant' / 'benchmarks' / 'postprocess_variants.txt'
    shell:
        """
        bash -c "{params.ACTIVATE} && \
        postprocess_variants --infile={params.Call} --multi_allelic_qual_filter='0.0' --vcf_stats_report='true' \
        --use_multiallelic_model='true' --regions={params.REGIONS} --ref={input.REF} --outfile={output.VCF} --cpus={threads}"
        """

rule variant_filtering:
    """
    Filtering the variants based on quality and depth, and getting the stats
    """
    input:
        VCF_gz = rules.postprocess_variants.output.VCF
    output:
        STATS = root / '{sample}' / 'deep_variant' / 'stats' / 'filtered_variants.stats.vchk',
        Filtered = root / '{sample}' / 'deep_variant' / 'Filtering' / 'filtered_variants.vcf.gz',
    params:
        VCF = lambda wildcards: root / wildcards.sample / 'deep_variant' / 'variants.vcf',
        STATS_RAW = lambda wildcards: root / wildcards.sample / 'deep_variant' / 'stats' / 'raw_variants.stats.vchk',
        QUAL = config["parameters"]["variant_filtering"]["min_QUAL"],
        DP = config["parameters"]["variant_filtering"]["min_DP"],
    threads: config["resources"]["variant_filtering"]["threads"]
    resources:
        mem_mb = config['resources']['variant_filtering']['memory']
    benchmark: root / '{sample}' / 'deep_variant' / 'benchmarks' / 'variant_filtering.txt'
    shell:
        """
        gunzip -kf {input.VCF_gz}
        bcftools stats {params.VCF} > {params.STATS_RAW}
        bcftools filter -i 'QUAL>={params.QUAL} & FORMAT/DP>={params.DP}' {params.VCF} -o {output.Filtered} -O z --write-index="tbi"
        gunzip -kf {output.Filtered}
        bcftools stats {output.Filtered} > {output.STATS}
        """


