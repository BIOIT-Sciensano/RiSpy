from .common import check_exists, check_tool

def validate_generation_ont_inputs(
    fastq,
    reference,
    hq_snp_positions,
    chain_file,
    chromosome_mapping,
    pickle_file,
    snp_info,
    cluster_json,
    representative_accessions
):
    check_exists(fastq, "FASTQ")
    check_exists(reference, "Reference genome")
    check_exists(hq_snp_positions, "High-quality SNP positions file")
    check_exists(chain_file, "Chain file for liftover")
    check_exists(chromosome_mapping, "Chromosome mapping file")
    check_exists(pickle_file, "Pickle file")
    check_exists(snp_info, "SNP information file")
    check_exists(cluster_json, "Cluster JSON file")
    check_exists(representative_accessions, "Representative accessions file")

    # Required tools
    for tool in ["NanoPlot", "seqkit", "chopper", "minimap2", "samtools", "bcftools", "qualimap", "liftOver", "bedtools", "cd-hit-est-2d"]:
        check_tool(tool)
