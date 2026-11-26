import typer
from .common import check_exists, check_tool
from rispy_cli.cli_utils import list_fingerprint_pairs

def validate_detection_ont_inputs(
    fastq,
    reference,
    chain_file,
    chromosome_mapping,
    centromere,
    barcodes,
    ref_vcf
):
    check_exists(fastq, "FASTQ")
    check_exists(reference, "Reference genome")
    check_exists(chain_file, "Chain file")
    check_exists(chromosome_mapping, "Chromosome mapping file")
    check_exists(centromere, "Centromere TSV file")

    # Required tools
    for tool in ["NanoPlot", "seqkit", "chopper", "minimap2", "samtools", "bcftools", "qualimap", "liftOver", "bedtools"]:
        check_tool(tool)

    # Fingerprint argument pairing
    if (barcodes is None) != (ref_vcf is None):
        raise typer.BadParameter("--barcodes and --ref-vcf must be given together, or both omitted.")

    if barcodes is not None:
        fingerprints = {barcodes.stem: {"barcodes": barcodes, "ref_vcf": ref_vcf}}

    if barcodes is None:
        fingerprints = list_fingerprint_pairs()

    return fingerprints

