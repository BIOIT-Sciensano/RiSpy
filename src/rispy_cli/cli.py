from pathlib import Path
import typer
from rispy_cli.cli_utils import (
    REFERENCE,
    CHAIN_FILE,
    CHROMOSOME_MAPPING,
    get_default_Illumina_model,
    get_default_ONT_model,
    Res,
    resource_options,
    call_command,
    is_test_data,
    fastp_params,
    chopper_params,
    calling_params,
    set_paths,
    run_pipeline,
    opt_output_dir,
    opt_cores,
    opt_mem_mb,
    opt_dry_run,
    opt_checkpoint,
    opt_fastq_r1,
    opt_fastq_r2,
    opt_fastq_ont,
    opt_max_barcodes,
    opt_hq_snp_positions,
    opt_pickle_file,
    opt_snp_info,
    opt_cluster_json,
    opt_representative_accessions,
    opt_barcodes,
    opt_ref_vcf,
    opt_top_n,
    opt_identification_threshold,
    opt_centromere,
    opt_fastp_window_size,
    opt_fastp_quality,
    opt_fastp_min_len,
    opt_chopper_quality,
    opt_chopper_min_len,
    opt_breadth_depth_cutoff,
    opt_min_qual,
    opt_min_dp,
    opt_regions,
    get_default_regions,
    get_default_generation_regions,
    run_info
)
from rispy_cli.validators.detection_illumina import validate_detection_illumina_inputs
from rispy_cli.validators.detection_ont import validate_detection_ont_inputs
from rispy_cli.validators.generation_illumina import validate_generation_illumina_inputs
from rispy_cli.validators.generation_ont import validate_generation_ont_inputs

app = typer.Typer(help="Fingerprint generation/detection from WGS data")

"""
-----------------------------
Per-rule resource options (threads / memory)
-----------------------------
"""

DEEPVARIANT_RES = [
    Res("make_examples", 16000, label="make-examples"),
    Res("call_variants", 16000, label="call-variants"),
    Res("postprocess_variants", 16000, 8, label="postprocess-variants"),
    Res("variant_filtering", 2000),
]
MAPPING_STATS_RES = [
    Res("calculate_mean_and_median_depth", 16000),
    Res("breadth_of_coverage", 2000),
    Res("count_bam", 2000, 8),
    Res("fold_80", 16000),
]
ILLUMINA_QC_RES = [
    Res("fastqc_raw", 8000, 4),
    Res("trim_reads", 16000, 8),
    Res("fastqc_trimmed", 8000, 4),
    Res("consensus_fasta_index", 4000, 2),
    Res("alignment", 16000, 8),
    Res("qualimap", 16000, 8),
    Res("samtools_depth", 2000, 4),
] + MAPPING_STATS_RES + DEEPVARIANT_RES
ONT_QC_RES = [
    Res("nanoplot_raw", 2000, 4),
    Res("extract_headers", 2000, 2),
    Res("filter_reads", 2000, 4),
    Res("trim_reads", 2000, 4),
    Res("nanoplot_after_filtering", 2000, 4),
    Res("index_reference", 5000, 4),
    Res("alignment", 16000, 8),
    Res("qualimap", 16000, 8),
    Res("samtools_depth", 2000, 2),
] + MAPPING_STATS_RES + DEEPVARIANT_RES
GENERATION_RES = [
    Res("hq_positions", 16000),
    Res("liftover", 2000),
    Res("chromosome_name_conversion", 2000),
    Res("hq_variants", 16000),
    Res("genotype_detection", 16000),
    Res("pickle_to_fasta", 16000),
    Res("clustering_cd_hit", 16000, 8),
    Res("update_cluster_json", 16000),
    Res("feature_selection", 16000),
    Res("add_info_to_indices", 2000),
    Res("list_remaining_accessions", 16000),
    Res("screen_2_snp_combinations", 16000),
    Res("add_info_to_snp_combinations", 2000),
    Res("unique_SNPs", 2000, flag="--unique-snps"),
    Res("bed_maker", 2000),
    Res("vcf_maker", 2000),
    Res("report_html", 2000),
]

def detection_res(liftover: int = 4000, filter_fingerprints: int = 4000) -> list:
    """Detection resources"""
    return [
        Res("unique_SNPs", 4000),
        Res("bed_maker", 2000),
        Res("liftover", liftover),
        Res("chromosome_name_conversion", 4000),
        Res("filter_fingerprints", filter_fingerprints),
        Res("validation_fingerprints_alt", 4000),
        Res("validation_fingerprints_ref", 4000),
        Res("barcodes_positions_new_ref", 4000),
        Res("chromosome_arms", 4000),
        Res("aggregate_fingerprint_matches", 4000),
        Res("report_html", 2000),
    ]

TEST_DATA_MESSAGE = "Test data detected: skipping fingerprint detection, running generation only."

"""
-----------------------------
Fingerprint Generation Illumina Full Pipeline Subcommand
-----------------------------
"""

@app.command("illugen")
@resource_options(ILLUMINA_QC_RES, GENERATION_RES)
def full_pipeline_illumina_generate(
    fastq_R1: Path = opt_fastq_r1(),
    fastq_R2: Path = opt_fastq_r2(),
    output_dir: Path = opt_output_dir(),
    call_variants_checkpoint_model: Path = opt_checkpoint(get_default_Illumina_model()),
    max_barcodes: int = opt_max_barcodes(),
    identification_threshold: int = opt_identification_threshold(),
    hq_snp_positions: Path = opt_hq_snp_positions(),
    pickle_file: Path = opt_pickle_file(),
    snp_info: Path = opt_snp_info(),
    cluster_json: Path = opt_cluster_json(),
    representative_accessions: Path = opt_representative_accessions(),
    cores: int = opt_cores(),
    mem_mb: int = opt_mem_mb(),
    dry_run: bool = opt_dry_run(),

    fastp_window_size: int = opt_fastp_window_size(),
    fastp_quality: int = opt_fastp_quality(),
    fastp_min_len: int = opt_fastp_min_len(),
    variant_filtering_min_QUAL: int = opt_min_qual(),
    variant_filtering_min_DP: int = opt_min_dp(),
    regions: Path = opt_regions(get_default_generation_regions()),
    breadth_of_coverage_depth_cutoff: int = opt_breadth_depth_cutoff(),
    resources: dict = None,
):
    """
    Running the complete Illumina pipeline from FASTQ files to fingerprint detection first and then generation.
    """
    if is_test_data(fastq_R1, fastq_R2):
        typer.echo(TEST_DATA_MESSAGE)
    else:
        call_command(full_pipeline_Illumina_detect, fastq_R1=fastq_R1, fastq_R2=fastq_R2, output_dir=output_dir,
                     call_variants_checkpoint_model=call_variants_checkpoint_model,
                     cores=cores, mem_mb=mem_mb, dry_run=dry_run,
                     fastp_window_size=fastp_window_size, fastp_quality=fastp_quality, fastp_min_len=fastp_min_len,
                     breadth_of_coverage_depth_cutoff=breadth_of_coverage_depth_cutoff,
                     variant_filtering_min_QUAL=variant_filtering_min_QUAL,
                     variant_filtering_min_DP=variant_filtering_min_DP,
                     identification_threshold=identification_threshold,
                     regions=regions,
                     **{f"{r.rule}_{k}": v for r in ILLUMINA_QC_RES for k, v in resources[r.rule].items()})
    validate_generation_illumina_inputs(
        fastq_R1=fastq_R1,
        fastq_R2=fastq_R2,
        reference=REFERENCE,
        hq_snp_positions=hq_snp_positions,
        chain_file=CHAIN_FILE,
        chromosome_mapping=CHROMOSOME_MAPPING,
        pickle_file=pickle_file,
        snp_info=snp_info,
        cluster_json=cluster_json,
        representative_accessions=representative_accessions
    )
    overrides = {
        "paths": set_paths(output_dir, fastq_R1=fastq_R1, fastq_R2=fastq_R2,
                           generation=(hq_snp_positions, pickle_file, snp_info, cluster_json, representative_accessions)),
        "parameters": {
            **fastp_params(fastp_window_size, fastp_quality, fastp_min_len),
            **calling_params(breadth_of_coverage_depth_cutoff, call_variants_checkpoint_model,
                             variant_filtering_min_QUAL, variant_filtering_min_DP, regions),
            "screen_2_snp_combinations": {"max_fingerprints": int(max_barcodes)},
            "identification_threshold": int(identification_threshold),
        },
        "resources": resources,
        "run_info": run_info(),
    }
    run_pipeline("Generation_Illumina_full_pipeline_config_template.yaml", overrides, output_dir,
                 ["illumina_qc_mapping_variant_calling.smk", "mapq.smk", "fingerprint_generation.smk"],
                 cores, mem_mb, dry_run, "config_generation.yaml")

"""
-----------------------------
Fingerprint Generation ONT Full Pipeline Subcommand
-----------------------------
"""

@app.command("ontgen")
@resource_options(ONT_QC_RES, GENERATION_RES)
def full_pipeline_ont_generate(
    fastq: Path = opt_fastq_ont(),
    output_dir: Path = opt_output_dir(),
    call_variants_checkpoint_model: Path = opt_checkpoint(get_default_ONT_model()),
    max_barcodes: int = opt_max_barcodes(),
    identification_threshold: int = opt_identification_threshold(),
    hq_snp_positions: Path = opt_hq_snp_positions(),
    pickle_file: Path = opt_pickle_file(),
    snp_info: Path = opt_snp_info(),
    cluster_json: Path = opt_cluster_json(),
    representative_accessions: Path = opt_representative_accessions(),
    cores: int = opt_cores(),
    mem_mb: int = opt_mem_mb(),
    dry_run: bool = opt_dry_run(),

    chopper_quality: int = opt_chopper_quality(),
    chopper_min_len: int = opt_chopper_min_len(),
    breadth_of_coverage_depth_cutoff: int = opt_breadth_depth_cutoff(),
    variant_filtering_min_QUAL: int = opt_min_qual(),
    variant_filtering_min_DP: int = opt_min_dp(),
    regions: Path = opt_regions(get_default_generation_regions()),
    resources: dict = None,
):
    """
    Running the complete ONT pipeline from FASTQ files to fingerprint detection first and then generation.
    """
    if is_test_data(fastq):
        typer.echo(TEST_DATA_MESSAGE)
    else:
        call_command(full_pipeline_ONT_detect, fastq=fastq, output_dir=output_dir,
                     call_variants_checkpoint_model=call_variants_checkpoint_model,
                     cores=cores, mem_mb=mem_mb, dry_run=dry_run,
                     chopper_quality=chopper_quality, chopper_min_len=chopper_min_len,
                     breadth_of_coverage_depth_cutoff=breadth_of_coverage_depth_cutoff,
                     variant_filtering_min_QUAL=variant_filtering_min_QUAL,
                     variant_filtering_min_DP=variant_filtering_min_DP,
                     identification_threshold=identification_threshold,
                     regions=regions,
                     **{f"{r.rule}_{k}": v for r in ONT_QC_RES for k, v in resources[r.rule].items()})
    validate_generation_ont_inputs(
        fastq=fastq,
        reference=REFERENCE,
        hq_snp_positions=hq_snp_positions,
        chain_file=CHAIN_FILE,
        chromosome_mapping=CHROMOSOME_MAPPING,
        pickle_file=pickle_file,
        snp_info=snp_info,
        cluster_json=cluster_json,
        representative_accessions=representative_accessions
    )
    overrides = {
        "paths": set_paths(output_dir, fastq=fastq,
                           generation=(hq_snp_positions, pickle_file, snp_info, cluster_json, representative_accessions)),
        "parameters": {
            **chopper_params(chopper_quality, chopper_min_len),
            **calling_params(breadth_of_coverage_depth_cutoff, call_variants_checkpoint_model,
                             variant_filtering_min_QUAL, variant_filtering_min_DP, regions),
            "screen_2_snp_combinations": {"max_fingerprints": int(max_barcodes)},
            "identification_threshold": int(identification_threshold),
        },
        "resources": resources,
        "run_info": run_info(),
    }
    run_pipeline("Generation_ONT_full_pipeline_config_template.yaml", overrides, output_dir,
                 ["ont_qc_mapping_variant_calling.smk", "mapq.smk", "fingerprint_generation.smk"],
                 cores, mem_mb, dry_run, "config_generation.yaml")

"""
-----------------------------
Fingerprint Detection Illumina Full Pipeline Subcommand
-----------------------------
"""

@app.command("illudetect")
@resource_options(ILLUMINA_QC_RES, detection_res())
def full_pipeline_Illumina_detect(
    fastq_R1: Path = opt_fastq_r1(),
    fastq_R2: Path = opt_fastq_r2(),
    output_dir: Path = opt_output_dir(),
    barcodes: Path = opt_barcodes(),
    ref_vcf: Path = opt_ref_vcf(),
    top_n: int = opt_top_n(),
    identification_threshold: int = opt_identification_threshold(),
    call_variants_checkpoint_model: Path = opt_checkpoint(get_default_Illumina_model()),
    dry_run: bool = opt_dry_run(),
    centromere: Path = opt_centromere(),
    cores: int = opt_cores(),
    mem_mb: int = opt_mem_mb(),

    fastp_window_size: int = opt_fastp_window_size(),
    fastp_quality: int = opt_fastp_quality(),
    fastp_min_len: int = opt_fastp_min_len(),
    breadth_of_coverage_depth_cutoff: int = opt_breadth_depth_cutoff(),
    variant_filtering_min_QUAL: int = opt_min_qual(),
    variant_filtering_min_DP: int = opt_min_dp(),
    regions: Path = opt_regions(get_default_regions()),
    resources: dict = None,
):
    """
    Running the complete Illumina pipeline from FASTQ files to fingerprint detection.
    """
    fingerprints = validate_detection_illumina_inputs(
        fastq_R1=fastq_R1,
        fastq_R2=fastq_R2,
        reference=REFERENCE,
        chain_file=CHAIN_FILE,
        chromosome_mapping=CHROMOSOME_MAPPING,
        centromere=centromere,
        barcodes=barcodes,
        ref_vcf=ref_vcf
    )
    overrides = {
        "paths": set_paths(output_dir, fastq_R1=fastq_R1, fastq_R2=fastq_R2,
                           centromere=centromere, fingerprints=fingerprints),
        "parameters": {
            **fastp_params(fastp_window_size, fastp_quality, fastp_min_len),
            **calling_params(breadth_of_coverage_depth_cutoff, call_variants_checkpoint_model,
                             variant_filtering_min_QUAL, variant_filtering_min_DP, regions),
            "top_n": int(top_n),
            "identification_threshold": int(identification_threshold),
        },
        "resources": resources,
        "run_info": run_info(),
    }
    run_pipeline("Detection_Illumina_full_pipeline_config_template.yaml", overrides, output_dir,
                 ["illumina_qc_mapping_variant_calling.smk", "mapq.smk", "fingerprint_detection.smk"],
                 cores, mem_mb, dry_run, "config_detect.yaml")

"""
-----------------------------
Fingerprint Detection ONT Full Pipeline Subcommand
-----------------------------
"""

@app.command("ontdetect")
@resource_options(ONT_QC_RES, detection_res(liftover=8000, filter_fingerprints=8000))
def full_pipeline_ONT_detect(
    fastq: Path = opt_fastq_ont(),
    output_dir: Path = opt_output_dir(),
    barcodes: Path = opt_barcodes(),
    ref_vcf: Path = opt_ref_vcf(),
    top_n: int = opt_top_n(),
    identification_threshold: int = opt_identification_threshold(),
    call_variants_checkpoint_model: Path = opt_checkpoint(get_default_ONT_model()),
    dry_run: bool = opt_dry_run(),
    centromere: Path = opt_centromere(),
    cores: int = opt_cores(),
    mem_mb: int = opt_mem_mb(),

    chopper_quality: int = opt_chopper_quality(),
    chopper_min_len: int = opt_chopper_min_len(),
    breadth_of_coverage_depth_cutoff: int = opt_breadth_depth_cutoff(),
    variant_filtering_min_QUAL: int = opt_min_qual(),
    variant_filtering_min_DP: int = opt_min_dp(),
    regions: Path = opt_regions(get_default_regions()),
    resources: dict = None,
):
    """
    Running the complete ONT pipeline from FASTQ files to fingerprint detection.
    """
    fingerprints = validate_detection_ont_inputs(
        fastq=fastq,
        reference=REFERENCE,
        chain_file=CHAIN_FILE,
        chromosome_mapping=CHROMOSOME_MAPPING,
        centromere=centromere,
        barcodes=barcodes,
        ref_vcf=ref_vcf
    )
    overrides = {
        "paths": set_paths(output_dir, fastq=fastq,
                           centromere=centromere, fingerprints=fingerprints),
        "parameters": {
            **chopper_params(chopper_quality, chopper_min_len),
            **calling_params(breadth_of_coverage_depth_cutoff, call_variants_checkpoint_model,
                             variant_filtering_min_QUAL, variant_filtering_min_DP, regions),
            "top_n": int(top_n),
            "identification_threshold": int(identification_threshold),
        },
        "resources": resources,
        "run_info": run_info(),
    }
    run_pipeline("Detection_ONT_full_pipeline_config_template.yaml", overrides, output_dir,
                 ["ont_qc_mapping_variant_calling.smk", "mapq.smk", "fingerprint_detection.smk"],
                 cores, mem_mb, dry_run, "config_detect.yaml")
