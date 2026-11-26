from .common import check_exists, check_tool

def validate_generation_illumina_inputs(
    fastq_R1,
    fastq_R2,
    reference,
    hq_snp_positions,
    chain_file,
    chromosome_mapping,
    pickle_file,
    snp_info,
    cluster_json,
    representative_accessions
):
    # Required files
    check_exists(fastq_R1, "FASTQ R1")
    check_exists(fastq_R2, "FASTQ R2")
    check_exists(reference, "Reference genome")
    check_exists(hq_snp_positions, "High-quality SNP positions file")
    check_exists(chain_file, "Chain file for liftover")
    check_exists(chromosome_mapping, "Chromosome mapping file")
    check_exists(pickle_file, "Pickle file")
    check_exists(snp_info, "SNP information file")
    check_exists(cluster_json, "Cluster JSON file")
    check_exists(representative_accessions, "Representative accessions file")

    # Required tools
    for tool in ["fastqc", "fastp", "bwa", "samtools", "bcftools", "qualimap", "liftOver", "bedtools", "cd-hit-est-2d"]:
        check_tool(tool)
