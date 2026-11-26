from Bio import SeqIO
from pathlib import Path
import re

def build_fq_by_sample(fastq_R1=None, fastq_R2=None, fastq=None):
    # Normalize "None" strings
    if fastq_R1 in ("None", "null", "NULL", ""):
        fastq_R1 = None
    if fastq_R2 in ("None", "null", "NULL", ""):
        fastq_R2 = None
    if fastq in ("None", "null", "NULL", ""):
        fastq = None
    # Case 1: nothing provided
    if fastq_R1 is None and fastq_R2 is None and fastq is None:
        raise ValueError(
            "No FASTQ input provided. Expected either:\n"
            "  • Illumina: fastq_R1 + fastq_R2\n"
            "  • ONT: fastq\n"
        )
    fq_by_sample = {}
    # ILLUMINA MODE
    if fastq_R1 is not None and fastq_R2 is not None:
        fastq_R1 = Path(fastq_R1)
        fastq_R2 = Path(fastq_R2)
        # Must end in .fastq.gz or .fq.gz
        gz_pattern = r"\.(?:fastq|fq)\.gz$"
        if not re.search(gz_pattern, fastq_R1.name, flags=re.IGNORECASE):
            raise ValueError(f"Invalid R1 filename (must end in .fastq.gz or .fq.gz): {fastq_R1.name}")
        if not re.search(gz_pattern, fastq_R2.name, flags=re.IGNORECASE):
            raise ValueError(f"Invalid R2 filename (must end in .fastq.gz or .fq.gz): {fastq_R2.name}")
        # read marker detection
        illumina_regex = re.compile(r"^(?P<prefix>.*?)(?:[_\-.])(?P<read>[Rr][12])(?P<suffix>.*?)\.(?:fastq|fq)\.gz$")
        m1 = illumina_regex.match(fastq_R1.name)
        m2 = illumina_regex.match(fastq_R2.name)
        if not m1:
            raise ValueError(f"Invalid Illumina R1 filename: {fastq_R1.name}")
        if not m2:
            raise ValueError(f"Invalid Illumina R2 filename: {fastq_R2.name}")
        # Validate correct read direction
        read1 = m1.group("read").lower()
        read2 = m2.group("read").lower()
        if read1 != "r1":
            raise ValueError(f"R1 file does not contain R1/r1 marker: {fastq_R1.name}")
        if read2 != "r2":
            raise ValueError(f"R2 file does not contain R2/r2 marker: {fastq_R2.name}")
        # Extract prefix/suffix
        prefix1 = m1.group("prefix")
        prefix2 = m2.group("prefix")
        suffix1 = m1.group("suffix")
        suffix2 = m2.group("suffix")
        # R1 sample
        if prefix1:
            sample_R1 = prefix1.strip("_-.")
        else:
            sample_R1 = suffix1.strip("_-.")
        if not sample_R1:
            sample_R1 = "sample"
        # R2 sample
        if prefix2:
            sample_R2 = prefix2.strip("_-.")
        else:
            sample_R2 = suffix2.strip("_-.")
        if not sample_R2:
            sample_R2 = "sample"
        # Validate sample consistency
        if sample_R1 != sample_R2:
            raise ValueError(
                f"Sample mismatch between R1 and R2:\n"
                f"  R1 sample: {sample_R1}\n"
                f"  R2 sample: {sample_R2}\n"
                "Both files must belong to the same sample."
            )
        fq_by_sample[sample_R1] = [fastq_R1, fastq_R2]
        return fq_by_sample
    # ONT MODE
    if fastq is not None:
        fastq = Path(fastq)
        # Must end in .fastq.gz or .fq.gz
        if not re.search(r"\.(?:fastq|fq)\.gz$", fastq.name, flags=re.IGNORECASE):
            raise ValueError(f"Invalid ONT filename (must end in .fastq.gz or .fq.gz): {fastq.name}")

        # ONT sample name = everything before extension
        ont_regex = re.compile(r"^(?P<sample>.*?)\.(?:fastq|fq)\.gz$", flags=re.IGNORECASE)
        m = ont_regex.match(fastq.name)
        if not m:
            raise ValueError(f"Invalid ONT filename: {fastq.name}")
        sample = m.group("sample").strip("_-.")
        if not sample:
            sample = "sample"
        fq_by_sample[sample] = [fastq]
        return fq_by_sample
    # Should never reach here
    raise ValueError("Invalid FASTQ input combination")

def rename_and_filter_fasta(input_fasta: str, output_fasta: str) -> None:
    """
    Define new chromosome names for the first 12 sequences
    """
    chr_names = [f"chr{i}" for i in range(1, 13)]
    with open(input_fasta, "r") as in_handle, open(output_fasta, "w") as out_handle:
        records = list(SeqIO.parse(in_handle, "fasta"))
        #Take only first 12 records (chromosomes)
        selected_records = records[:12]
        for rec, new_name in zip(selected_records, chr_names):
            rec.id = new_name
            #remove description to keep header clean
            rec.description = ""
        SeqIO.write(selected_records, out_handle, "fasta")

