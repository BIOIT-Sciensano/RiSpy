from pathlib import Path
from rispy_cli.io import build_fq_by_sample
from rispy_cli.cli_utils import deepvariant_shell_activation

DV_ACTIVATE = deepvariant_shell_activation()

root = Path(config["paths"]["output_dir"])

fastq_R1 = Path(config["paths"]["fastq_R1"])

fastq_R2 = Path(config["paths"]["fastq_R2"])

ref = config["paths"]["reference"]

fq_by_sample = build_fq_by_sample(fastq_R1=fastq_R1,fastq_R2=fastq_R2)

rule all:
    input:
        FASTQC_Trimmed_1P = expand(str(root / '{sample}' / 'fastqc_after_trimming' / 'trimmed_reads_1P_fastqc.html'),sample=fq_by_sample.keys()),
        FASTQC_Raw = expand(str(root / '{sample}' / 'fastqc_before_trimming'),sample=fq_by_sample.keys()),
        Filtered = expand(str(root / '{sample}' / 'deep_variant' / 'Filtering' / 'filtered_variants.vcf.gz'),sample=fq_by_sample.keys()),

rule fastqc_raw:
    """
    Quality control of raw reads
    """
    input:
        FQ_fwd = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
        FQ_rev = lambda wildcards: str(fq_by_sample[wildcards.sample][1])
    output:
        outdir = directory(root / "{sample}" / "fastqc_before_trimming")
    threads: config["resources"]["fastqc_raw"]["threads"]
    resources:
        mem_mb = config["resources"]["fastqc_raw"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'fastqc_raw.txt'
    shell:
        """
        mkdir -p {output.outdir}
        fastqc -t {threads} {input.FQ_fwd} {input.FQ_rev} -o {output.outdir}
        """

rule trim_reads:
    """
    Trimming of paired end reads
    """
    input:
        FQ_fwd = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
        FQ_rev = lambda wildcards: str(fq_by_sample[wildcards.sample][1])
    output:
        FQ_1P = root / '{sample}' / 'trimming' / 'trimmed_reads_1P.fastq.gz',
        FQ_2P = root / '{sample}' / 'trimming' / 'trimmed_reads_2P.fastq.gz',
    params:
        window_size = config["parameters"]["fastp"]["window_size"],
        quality = config["parameters"]["fastp"]["quality"],
        min_len = config["parameters"]["fastp"]["min_len"],
        json = lambda wildcards: root / wildcards.sample / 'trimming' / 'fastp.json',
        html = lambda wildcards: root / wildcards.sample / 'trimming' / 'fastp.html'
    threads: config["resources"]["trim_reads"]["threads"]
    resources:
        mem_mb = config["resources"]["trim_reads"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'trim_reads.txt'
    shell:
        """
        fastp --cut_right --cut_right_window_size {params.window_size} --cut_right_mean_quality {params.quality} \
        --length_required {params.min_len} --detect_adapter_for_pe --thread {threads} -j {params.json} -h {params.html} \
        -i {input.FQ_fwd} -I {input.FQ_rev} -o {output.FQ_1P} -O {output.FQ_2P}
        """

rule fastqc_trimmed:
    """
    Quality control of trimmed reads
    """
    input:
        FQ_1P = rules.trim_reads.output.FQ_1P,
        FQ_2P = rules.trim_reads.output.FQ_2P,
    output:
        FASTQC_Trimmed_1P = root / '{sample}' / 'fastqc_after_trimming' / 'trimmed_reads_1P_fastqc.html',
        FASTQC_Trimmed_2P = root / '{sample}' / 'fastqc_after_trimming' / 'trimmed_reads_2P_fastqc.html',
    params:
        outdir = lambda wildcards: root / wildcards.sample / 'fastqc_after_trimming',
    threads: config["resources"]["fastqc_trimmed"]["threads"]
    resources:
        mem_mb = config["resources"]["fastqc_trimmed"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'fastqc_trimmed.txt'
    shell:
        """
        fastqc -t {threads} {input.FQ_1P} {input.FQ_2P} -o {params.outdir}
        """

rule consensus_fasta_index:
    """
    Indexing the reference fasta
    """
    input:
        FASTA = ref,
    output:
        REFERENCE = root / 'reference' / 'reference.fasta'
    threads: config["resources"]["consensus_fasta_index"]["threads"]
    resources:
        mem_mb = config["resources"]["consensus_fasta_index"]["memory"]
    benchmark: root / 'reference' / 'benchmarks' / 'consensus_fasta_index.txt'
    shell:
        """
        cp {input.FASTA} {output.REFERENCE}
        bwa index {output.REFERENCE}
        samtools faidx -@ {threads} {output.REFERENCE}
        """

rule alignment:
    """
    Mapping the short reads to reference and convert to bam, then sorting and indexing the bam file
    """
    input:
        FASTA = rules.consensus_fasta_index.output.REFERENCE,
        FQ_1P = rules.trim_reads.output.FQ_1P,
        FQ_2P = rules.trim_reads.output.FQ_2P,
    output:
        BAM = root / '{sample}' / 'mapping' / 'mapped_reads_sorted.bam'
    threads: config["resources"]["alignment"]["threads"]
    params:
        outdir = root
    resources:
        mem_mb = config["resources"]["alignment"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'alignment.txt'
    shell:
        """
        bwa mem -K 100000000 -v 3 -t {threads} -k 19 -r 1.5 -M {input.FASTA} {input.FQ_1P} {input.FQ_2P} \
        | samtools view -b - -@ {threads} | samtools sort - -o {output.BAM} -T {params.outdir} -@ {threads}
        samtools index -@ {threads} {output.BAM}
        """

rule make_examples:
    """
    Running make_examples to generate examples for deepvariant
    """
    input:
        BAM = rules.alignment.output.BAM,
        REF = rules.consensus_fasta_index.output.REFERENCE
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
        --multi_allelic_mode='include_het_alt_images' --select_variant_types='all' \
        --regions={params.REGIONS} --vsc_min_count_snps='1' --vsc_min_fraction_snps='0.001' \
        --vsc_min_count_indels='1' --vsc_min_fraction_indels='0.001' \
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
        REF = rules.consensus_fasta_index.output.REFERENCE
    output:
        VCF = root  / '{sample}' / 'deep_variant' / 'variants.vcf.gz'
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
        mem_mb = config["resources"]["variant_filtering"]["memory"]
    benchmark: root / '{sample}' / 'deep_variant' / 'benchmarks' / 'variant_filtering.txt'
    shell:
        """
        gunzip -kf {input.VCF_gz}
        bcftools stats {params.VCF} > {params.STATS_RAW}
        bcftools filter -i 'QUAL>={params.QUAL} & FORMAT/DP>={params.DP}' {params.VCF} -o {output.Filtered} -O z --write-index="tbi"
        gunzip -kf {output.Filtered}
        bcftools stats {output.Filtered} > {output.STATS}
        """

