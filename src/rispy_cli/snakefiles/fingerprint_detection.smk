from pathlib import Path
import os
from rispy_cli.fingerprint_detection_pipeline import (
    detect_unique_snps,
    bed_maker,
    chromosome_name_conversion,
    validation_fingerprints_alt,
    validation_fingerprints_ref,
    barcodes_positions_new_ref,
    chromosome_arms,
    summarize_fingerprint_matches
)
from rispy_cli.io import build_fq_by_sample
from rispy_cli.generate_html_fingerprint_detection import generate_report_html

fastq_R1 = config["paths"].get("fastq_R1")
fastq_R2 = config["paths"].get("fastq_R2")
fastq    = config["paths"].get("fastq")

fq_by_sample = build_fq_by_sample(fastq_R1=fastq_R1,fastq_R2=fastq_R2,fastq=fastq)

root = Path(config["paths"]["output_dir"])

chain_file_old_nip_to_new_nip_ref = config["paths"]["chain_file"]

chromosome_mapping_file = config["paths"]["chromosome_mapping"]

centromere = config["paths"]["centromere"]

query_vcf = config["paths"].get("query_vcf", None)

bam = config["paths"].get("bam", None)

fingerprints = config["paths"]["fingerprints"]  # {name: {"barcodes":..., "ref_vcf":...}}

fp_names = list(fingerprints.keys())

top_n = config["parameters"]["top_n"]

report_html_file = root / 'report.html'

# Inside illugen/ontgen the detection results go into the combined generation report, so no separate detection report
make_report = not str((config.get("run_info") or {}).get("command", "")).endswith("gen")

# File the generation report reads: the summary (several fingerprints) or chr_arms.txt (one fingerprint, completes its config_results.yaml)
detection_result_files = ([str(root / s / 'Fingerprint_Detection' / 'fingerprint_summary.tsv') for s in fq_by_sample.keys()] if len(fp_names) > 1
                          else [str(root / s / 'Fingerprint_Detection' / fp_names[0] / 'Final_Barcodes' / 'chr_arms.txt') for s in fq_by_sample.keys()])

def barcodes_for(fp): return fingerprints[fp]["barcodes"]

def ref_vcf_for(fp):  return fingerprints[fp]["ref_vcf"]

rule all:
    input:
        Detected_Barcodes = expand(str(root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Final_Barcodes' / 'detected_barcodes.tsv'),sample=fq_by_sample.keys(),fp=fp_names),
        Detection_Result = detection_result_files,
        HTML = [str(report_html_file)] if make_report else []

rule unique_SNPs:
    """
    Detecting the unique SNPs in the combinations of 2-SNV barcodes
    """
    input:
        Fingerprints = lambda wildcards: barcodes_for(wildcards.fp)
    output:
        Unique_SNPs = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'unique_snps.tsv',
    threads: config["resources"]["unique_SNPs"]["threads"]
    resources:
        mem_mb = config["resources"]["unique_SNPs"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'benchmarks' / 'unique_SNPs.txt'
    run:
        detect_unique_snps(input.Fingerprints,output.Unique_SNPs)


rule bed_maker:
    """
    Making a bed file out of the unique SNPs 
    """
    input:
        FINGERPRINTS = rules.unique_SNPs.output.Unique_SNPs,
    output:
        FINGER_BED = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'fingerprints.bed'
    threads: config["resources"]["bed_maker"]["threads"]
    resources:
        mem_mb = config["resources"]["bed_maker"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'benchmarks' / 'bed_maker.txt'
    run:
        bed_maker(input.FINGERPRINTS,output.FINGER_BED)

rule liftover:
    """
    Executing liftOver to get the SNP positions in the new reference genome
    """
    input:
        FINGER_BED = rules.bed_maker.output.FINGER_BED,
        CHAINFILE = chain_file_old_nip_to_new_nip_ref
    output:
        FINGER_BED_LIFTOVER = root / '{sample}' / 'Fingerprint_Detection'  / '{fp}' / 'fingerprints_liftover.bed',
        UNMAPPED = root / '{sample}' / 'Fingerprint_Detection'  / '{fp}' / 'unmapped.bed'
    threads: config["resources"]["liftover"]["threads"]
    resources:
        mem_mb = config["resources"]["liftover"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'benchmarks' / 'liftover.txt'
    shell:
        """
        liftOver {input.FINGER_BED} {input.CHAINFILE} {output.FINGER_BED_LIFTOVER} {output.UNMAPPED}
        """

rule chromosome_name_conversion:
    """
    Converting the chromosome names to match the new reference genome
    """
    input:
        CHROM_MAP = chromosome_mapping_file,
        FINGER_BED_LIFTOVER = rules.liftover.output.FINGER_BED_LIFTOVER
    output:
        FINGER_BED_CONVERT_CHR = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'fingerprints_liftover_converted_chr_names.bed'
    threads: config["resources"]["chromosome_name_conversion"]["threads"]
    resources:
        mem_mb = config["resources"]["chromosome_name_conversion"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'benchmarks' / 'chromosome_name_conversion.txt'
    run:
        chromosome_name_conversion(input.CHROM_MAP, input.FINGER_BED_LIFTOVER, output.FINGER_BED_CONVERT_CHR)


rule filter_fingerprints:
    """
    Making two separate bed files for reference and alternative alleles
    Note: The decision of a SNP being reference or alternative allele is made based on the new reference genome 
    """
    input:
        FINGER_BED_CONVERT_CHR = rules.chromosome_name_conversion.output.FINGER_BED_CONVERT_CHR,
        VCF_STRICT = lambda wildcards: ref_vcf_for(wildcards.fp)
    output:
        REF_BED = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'reference_allele_fingerprints.bed',
        ALT_BED = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'alternative_allele_fingerprints.bed',
    threads: config["resources"]["filter_fingerprints"]["threads"]
    resources:
        mem_mb = config["resources"]["filter_fingerprints"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'benchmarks' / 'filter_fingerprints.txt'
    shell:
        """
        bedtools subtract -a {input.FINGER_BED_CONVERT_CHR} -b {input.VCF_STRICT} > {output.REF_BED}
        bedtools intersect -a {input.FINGER_BED_CONVERT_CHR} -b {input.VCF_STRICT} > {output.ALT_BED}
        """

rule validation_fingerprints_alt:
    """
    Checking if we can find the alternative SNPs in the samples
    """
    input:
        FINGER_BED_CONVERT_CHR = rules.filter_fingerprints.output.ALT_BED,
        VCF = lambda wildcards: query_vcf or root / wildcards.sample / "deep_variant" / "Filtering" / "filtered_variants.vcf.gz",
        Barcodes = lambda wildcards: barcodes_for(wildcards.fp)
    output:
        Detected_snps = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Alternative_Alleles' / 'detected_snps.txt',
        Undetected_snps = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Alternative_Alleles' / 'undetected_snps.txt',
        Detected_Barcodes = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Alternative_Alleles' / 'tmp_detected_barcodes.tsv',
    params:
        config = lambda wildcards: root / wildcards.sample / 'Fingerprint_Detection' / wildcards.fp / 'config_results.yaml',
        config_dir = root
    threads: config["resources"]["validation_fingerprints_alt"]["threads"]
    resources:
        mem_mb = config["resources"]["validation_fingerprints_alt"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Alternative_Alleles' / 'benchmarks' / 'validation_fingerprints_alt.txt'
    run:
        os.makedirs(str(params.config_dir), exist_ok=True)
        validation_fingerprints_alt(input.FINGER_BED_CONVERT_CHR, input.VCF, input.Barcodes, output.Detected_snps, output.Undetected_snps, output.Detected_Barcodes, params.config)


rule validation_fingerprints_ref:
    """
    Checking if we can find the reference SNPs in the samples
    """
    input:
        FINGER_BED_CONVERT_CHR = rules.filter_fingerprints.output.REF_BED,
        VCF = lambda wildcards: query_vcf or root / wildcards.sample / "deep_variant" / "Filtering" / "filtered_variants.vcf.gz",
        BAM = lambda wildcards: bam or root / wildcards.sample / 'mapping' / 'mapped_reads_sorted.bam',
        Barcodes = rules.validation_fingerprints_alt.output.Detected_Barcodes
    output:
        Detected_snps = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Reference_Alleles' / 'detected_snps.txt',
        Undetected_snps = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Reference_Alleles' / 'undetected_snps.txt',
        Detected_Barcodes = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Final_Barcodes' / 'detected_barcodes.tsv'
    params:
        config = lambda wildcards: root / wildcards.sample / 'Fingerprint_Detection' / wildcards.fp / 'config_results.yaml'
    threads: config["resources"]["validation_fingerprints_ref"]["threads"]
    resources:
        mem_mb = config["resources"]["validation_fingerprints_ref"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Reference_Alleles' / 'benchmarks' / 'validation_fingerprints_ref.txt'
    run:
        validation_fingerprints_ref(input.FINGER_BED_CONVERT_CHR, input.VCF, input.BAM, input.Barcodes, output.Detected_snps, output.Undetected_snps, output.Detected_Barcodes, params.config)


rule barcodes_positions_new_ref:
    """
    Generating the final set of fingerprint that could be detected in the samples with their positions in the new reference genome
    """
    input:
        BARCODES = rules.validation_fingerprints_ref.output.Detected_Barcodes,
        SNPs = rules.liftover.output.FINGER_BED_LIFTOVER
    output:
        BARCODES_LIFTOVER = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Final_Barcodes' / 'detected_barcodes.liftover.tsv'
    threads: config["resources"]["barcodes_positions_new_ref"]["threads"]
    resources:
        mem_mb = config["resources"]["barcodes_positions_new_ref"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Final_Barcodes' / 'benchmarks' / 'barcodes_positions_new_ref.txt'
    run:
        barcodes_positions_new_ref(input.BARCODES, input.SNPs, output.BARCODES_LIFTOVER)


rule chromosome_arms:
    """
    Checking which chromosome arms the SNPs are located in
    """
    input:
        CENTROMERE = centromere,
        BARCODES = rules.barcodes_positions_new_ref.output.BARCODES_LIFTOVER,
    output:
        OUT_TXT = root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Final_Barcodes' / 'chr_arms.txt'
    params:
        config = lambda wildcards: root / wildcards.sample / 'Fingerprint_Detection' / wildcards.fp / 'config_results.yaml'
    threads: config["resources"]["chromosome_arms"]["threads"]
    resources:
        mem_mb = config["resources"]["chromosome_arms"]["memory"]
    benchmark: root / '{sample}' / 'Fingerprint_Detection' / '{fp}' / 'Final_Barcodes' / 'benchmarks' / 'chromosome_arms.txt'
    run:
        chromosome_arms(input.CENTROMERE, input.BARCODES, output.OUT_TXT, params.config)

if len(fp_names) > 1:
    rule aggregate_fingerprint_matches:
        input:
            Detected = lambda wildcards: expand(str(root / wildcards.sample / 'Fingerprint_Detection' / '{fp}' / 'Final_Barcodes' / 'detected_barcodes.tsv'),fp=fp_names),
            CHR_ARMS = lambda wildcards: expand(str(root / wildcards.sample / 'Fingerprint_Detection' / '{fp}' / 'Final_Barcodes' / 'chr_arms.txt'),fp=fp_names),
        output:
            SUMMARY = root / '{sample}' / 'Fingerprint_Detection' / 'fingerprint_summary.tsv'
        params:
            fp_names = fp_names
        threads: config["resources"]["aggregate_fingerprint_matches"]["threads"]
        resources:
            mem_mb = config["resources"]["aggregate_fingerprint_matches"]["memory"]
        run:
            summarize_fingerprint_matches(list(zip(params.fp_names,input.Detected,[barcodes_for(fp) for fp in params.fp_names])),output.SUMMARY)

if len(fp_names) > 1:
    sample = list(fq_by_sample.keys())[0]
    rule report_html:
        input:
            SUMMARY = root / sample / 'Fingerprint_Detection' / 'fingerprint_summary.tsv'
        output:
            HTML = report_html_file
        params:
            config_qc = (p if p.exists() else None) if (p := root / sample / "config_results.yaml") else None,
            top_n = top_n
        threads: config["resources"]["report_html"]["threads"]
        resources:
            mem_mb = config["resources"]["report_html"]["memory"]
        run:
            generate_report_html(output.HTML,params.config_qc,Path(input.SUMMARY),params.top_n,config)

else:
    only_fp = fp_names[0]
    sample = list(fq_by_sample.keys())[0]
    rule report_html:
        input:
            Check = lambda wildcards: root / sample / 'Fingerprint_Detection' / only_fp / 'Final_Barcodes' / 'chr_arms.txt'
        output:
            HTML = report_html_file
        params:
            config_qc = (p if p.exists() else None) if (p := root / sample / "config_results.yaml") else None,
            top_n = top_n,
            SUMMARY = root / sample / 'Fingerprint_Detection' / only_fp / 'config_results.yaml'
        threads: config["resources"]["report_html"]["threads"]
        resources:
            mem_mb = config["resources"]["report_html"]["memory"]
        run:
            generate_report_html(output.HTML,params.config_qc,params.SUMMARY,params.top_n,config)


