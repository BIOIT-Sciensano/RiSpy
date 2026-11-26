from pathlib import Path
import sys
import functools
import inspect
import importlib.resources as resources
import os
from typing import NamedTuple, Optional
import yaml
import typer
from typer.models import ParameterInfo
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import rispy_cli
from rispy_cli.generate_html_fingerprint_detection import MATCH_THRESHOLD

SNAKEFILES = Path(__file__).parent / "snakefiles"

"""
-----------------------------
Packaged data files
-----------------------------
"""

def get_default_centromere_path() -> Path:
    return resources.files("rispy_cli.data") / "Centromere.tsv"

def get_default_hq_positions() -> Path:
    return resources.files("rispy_cli.data") / "3K-HDRA-snp-comm-miss5pc_index_position.tsv"

def get_default_reference() -> Path:
    return resources.files("rispy_cli.data") / "GCF_034140825.1_ASM3414082v1_genomic.fna"

def get_default_chain_file() -> Path:
    return resources.files("rispy_cli.data") / "CHAINFILE.chain"

def get_default_chromosome_mapping() -> Path:
    return resources.files("rispy_cli.data") / "chromosome_mapping.txt"

def get_default_pickle_file() -> Path:
    return resources.files("rispy_cli.data") / "hq_genotype.pkl"

def get_default_SNP_Info() -> Path:
    return resources.files("rispy_cli.data") / "3kall_snpposition_map.tsv"

def get_default_Cluster_json() -> Path:
    return resources.files("rispy_cli.data") / "clusters.json"

def get_default_representative_accessions() -> Path:
    return resources.files("rispy_cli.data") / "representatives.fasta"

def get_default_Illumina_model() -> Path:
    return resources.files("rispy_cli.deepvariant") / "1.10.0" / "checkpoints" / "wgs" / "deepvariant.wgs.ckpt"

def get_default_ONT_model() -> Path:
    return resources.files("rispy_cli.deepvariant") / "1.10.0" / "checkpoints" / "ont" / "deepvariant.ont.ckpt"

def get_default_fingerprints_dir() -> Path:
    return Path(str(resources.files("rispy_cli") / "fingerprints"))

def get_default_regions() -> Path:
    return resources.files("rispy_cli.data") / "Fingerprint_SNVs.bed"

def get_default_generation_regions() -> Path:
    return resources.files("rispy_cli.data") / "Generation_SNVs.bed"

# Fixed inputs: not exposed on the command line, the user cannot change them.
REFERENCE = get_default_reference()
CHAIN_FILE = get_default_chain_file()
CHROMOSOME_MAPPING = get_default_chromosome_mapping()
TEST_REFERENCE = resources.files("rispy_cli") / "test_data" / "test_reference.fasta"

def list_fingerprint_pairs(fingerprints_dir: Path | None = None) -> dict[str, dict[str, Path]]:
    """
    Scan the fingerprints directory and pair up name.tsv (barcodes) with
    name.vcf (ref_vcf) by shared basename.
    Returns {name: {"barcodes": Path, "ref_vcf": Path}}
    """
    fingerprints_dir = fingerprints_dir or get_default_fingerprints_dir()
    tsvs = {p.stem: p for p in fingerprints_dir.glob("*.tsv")}
    vcfs = {p.stem: p for p in fingerprints_dir.glob("*.vcf")}
    missing_vcf = sorted(set(tsvs) - set(vcfs))
    missing_tsv = sorted(set(vcfs) - set(tsvs))
    if missing_vcf:
        raise FileNotFoundError(f"No matching .vcf for barcode file(s): {missing_vcf}")
    if missing_tsv:
        raise FileNotFoundError(f"No matching .tsv for vcf file(s): {missing_tsv}")
    return {name: {"barcodes": tsvs[name], "ref_vcf": vcfs[name]} for name in sorted(tsvs)}

def deepvariant_shell_activation():
    """
    Puts DeepVariant on PATH so that:
      - make_examples, call_variants, postprocess_variants are on PATH
      - DeepVariant's venv (Python 3.10) comes first on PATH
      - PYTHONPATH is cleared so no outer environment leaks packages in
    The launchers are patched at install time to use the venv's python by
    absolute path, so no `source activate` is needed.
    """
    base = resources.files("rispy_cli.deepvariant") / "1.10.0"
    bin_dir = base / "bazel-out" / "k8-opt" / "bin" / "deepvariant"
    venv_bin = base / "venv" / "bin"
    return (
        f"export PATH={venv_bin}:{bin_dir}:$PATH && "
        f"unset PYTHONPATH"
    )

"""
-----------------------------
Typer helpers
-----------------------------
"""

def option_default(value):
    """Typer default object -> its plain value; anything else is returned unchanged."""
    return value.default if isinstance(value, ParameterInfo) else value

def call_command(func, **kwargs):
    """
    Call a Typer command function, filling any argument with the default declared in its typer.Option(...).
    """
    bound = {name: kwargs[name] if name in kwargs else option_default(param.default)
             for name, param in inspect.signature(func).parameters.items()}
    return func(**bound)

INPUT_FLAGS = {"--fastq-R1", "-f1", "--fastq-R2", "-f2", "--fastq", "-f", "--output-dir", "-o"}

def run_info() -> dict:
    """
    Name of the command and the options the user typed (input files and output directory excluded),
    read from the command line so it also works when illugen/ontgen run detection internally.
    Stored in the config for the HTML reports.
    """
    args = sys.argv[1:]
    if not args:
        return {}
    options = {}
    rest = args[1:]
    i = 0
    while i < len(rest):
        token = rest[i]
        i += 1
        if not token.startswith("-"):
            continue
        flag, _, value = token.partition("=")
        if not value and i < len(rest) and not rest[i].startswith("-"):
            value = rest[i]
            i += 1
        if flag not in INPUT_FLAGS:
            options[flag] = value
    return {"command": args[0], "non_default_parameters": options}

class Res(NamedTuple):
    rule: str                      # Snakemake rule name (also the Python parameter prefix)
    memory: int                    # default memory (MB)
    threads: Optional[int] = None  # default threads; None = rule has no threads option
    flag: Optional[str] = None     # CLI flag prefix if it differs from "--<rule-with-dashes>"
    label: Optional[str] = None    # name shown in --help if it differs from the rule name

def resource_options(*groups):
    """
    Appends one --<rule>-threads / --<rule>-memory option per Res to the command,
    then hands the collected values to the command as `resources={rule: {...}}`
    (ready to be written to the config as-is).
    """
    spec = [res for group in groups for res in group]
    def decorate(func):
        sig = inspect.signature(func)
        base = [p for p in sig.parameters.values() if p.name != "resources"]
        extra = []
        for r in spec:
            flag = r.flag or "--" + r.rule.replace("_", "-")
            label = r.label or r.rule
            if r.threads is not None:
                extra.append(inspect.Parameter(
                    f"{r.rule}_threads", inspect.Parameter.KEYWORD_ONLY, annotation=int,
                    default=typer.Option(r.threads, f"{flag}-threads", help=f"threads for {label}")))
            extra.append(inspect.Parameter(
                f"{r.rule}_memory", inspect.Parameter.KEYWORD_ONLY, annotation=int,
                default=typer.Option(r.memory, f"{flag}-memory", help=f"memory (MB) for {label}")))
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            resources_dict = {}
            for r in spec:
                entry = {}
                if r.threads is not None:
                    entry["threads"] = int(option_default(kwargs.pop(f"{r.rule}_threads", r.threads)))
                entry["memory"] = int(option_default(kwargs.pop(f"{r.rule}_memory", r.memory)))
                resources_dict[r.rule] = entry
            return func(*args, resources=resources_dict, **kwargs)
        wrapper.__signature__ = sig.replace(parameters=base + extra)
        wrapper.__annotations__ = {**{k: v for k, v in func.__annotations__.items() if k != "resources"},
                                   **{p.name: int for p in extra}}
        return wrapper
    return decorate

"""
-----------------------------
Shared options (identical flags/defaults/help in every command that uses them)
-----------------------------
"""

def opt_output_dir():
    return typer.Option(..., "--output-dir", "-o", help="Output directory")

def opt_cores():
    return typer.Option(8, "--threads", "-t", help="Number of CPU cores")

def opt_mem_mb():
    return typer.Option(16000, "--mem-mb", "-m", help="Memory in MB")

def opt_dry_run():
    return typer.Option(False, "--dry-run", help="Generate the config file and perform a dry-run without executing rules")

def opt_checkpoint(default):
    return typer.Option(default, "--call-variants-checkpoint-model", "-vcm", help="Call variants checkpoint model path")

def opt_fastq_r1():
    return typer.Option(..., "--fastq-R1", "-f1", help="Paired-end Illumina FASTQ file: forward reads (R1)")

def opt_fastq_r2():
    return typer.Option(..., "--fastq-R2", "-f2", help="Paired-end Illumina FASTQ file: reverse reads (R2)")

def opt_fastq_ont():
    return typer.Option(..., "--fastq", "-f", help="ONT FASTQ file")

# generation inputs
def opt_max_barcodes():
    return typer.Option(100, "--max-barcodes", "-mb", min=1, help="Maximum number of 2-SNV barcodes to generate (default: 100)")

def opt_hq_snp_positions():
    return typer.Option(get_default_hq_positions(), "--hq-snp-positions", "-hsp", help="Tab‑delimited file listing SNP index, chromosome, and genomic position for 160K high‑quality SNPs")

def opt_pickle_file():
    return typer.Option(get_default_pickle_file(), "--pickle-file", "-pf", help="Pickle file containing the high‑quality genotype matrix with a 'VARIETY' column and SNP genotype columns")

def opt_snp_info():
    return typer.Option(get_default_SNP_Info(), "--snp-info", "-si", help="Tab‑delimited SNP table with columns: SNP_INDEX, CHROMOSOME, POSITION, REFCALL for 32 million SNPs")

def opt_cluster_json():
    return typer.Option(get_default_Cluster_json(), "--cluster-json", "-cj", help="JSON file mapping cluster names to lists of accession IDs")

def opt_representative_accessions():
    return typer.Option(get_default_representative_accessions(), "--representative-accessions", "-ra", help="FASTA file containing representative sequences for each cluster")

# detection inputs
def opt_barcodes():
    return typer.Option(None, "--barcodes", "-b", exists=True, dir_okay=False, readable=True,
                        help="tsv file of the fingerprint (2-SNV barcodes). If omitted, all 3KRG fingerprints are tested.")

def opt_ref_vcf():
    return typer.Option(None, "--ref-vcf", "-v", exists=True, dir_okay=False, readable=True,
                        help="VCF used to generate the fingerprint. If omitted, all 3KRG fingerprints are tested.")

def opt_top_n():
    return typer.Option(10, "--top-n", "-n", min=1, help="Number of top-matching rice lines to show in the report.")

def opt_identification_threshold():
    return typer.Option(MATCH_THRESHOLD, "--identification-threshold", "-it", min=1,
                        help="Minimum number of detected barcodes for a line to be considered identified")
def opt_centromere():
    return typer.Option(get_default_centromere_path(), "--centromere", "-c", help="Centromere TSV file")

# QC / variant-calling parameters
def opt_fastp_window_size():
    return typer.Option(4, "--fastp-window-size", help="window size for fastp")

def opt_fastp_quality():
    return typer.Option(20, "--fastp-quality", help="quality cutoff for fastp")

def opt_fastp_min_len():
    return typer.Option(40, "--fastp-min-len", help="minimum length for fastp")

def opt_chopper_quality():
    return typer.Option(10, "--chopper-quality", help="quality cutoff for chopper")

def opt_chopper_min_len():
    return typer.Option(150, "--chopper-min-len", help="minimum length for chopper")

def opt_breadth_depth_cutoff():
    return typer.Option(5, "--breadth-of-coverage-depth-cutoff", help="depth cutoff for calculating the breadth of coverage")

def opt_min_qual():
    return typer.Option(20, "--variant-filtering-min-QUAL", help="minimum QUAL for variant filtering")

def opt_min_dp():
    return typer.Option(5, "--variant-filtering-min-DP", help="minimum DP for variant filtering")


def opt_regions(default):
    return typer.Option(default, "--regions", "-r", exists=True, dir_okay=False, readable=True,
                        help="BED file restricting variant calling (DeepVariant) to these regions")
"""
-----------------------------
Config building blocks
-----------------------------
"""

def fastp_params(window_size, quality, min_len) -> dict:
    return {"fastp": {"window_size": int(window_size), "quality": int(quality), "min_len": int(min_len)}}

def chopper_params(quality, min_len) -> dict:
    return {"chopper": {"quality": int(quality), "min_len": int(min_len)}}

def calling_params(depth_cutoff, checkpoint, min_qual, min_dp, regions) -> dict:
    return {
        "breadth_of_coverage": {"depth_cutoff": int(depth_cutoff)},
        "call_variants": {"CHECKPOINT": str(checkpoint), "regions": str(regions)},
        "variant_filtering": {"min_QUAL": int(min_qual), "min_DP": int(min_dp)},
    }

def fingerprints_to_config(fingerprints: dict | None) -> dict | None:
    """{name: {"barcodes": Path, "ref_vcf": Path}} -> same structure with string paths."""
    if fingerprints is None:
        return None
    return {name: {"barcodes": str(v["barcodes"]), "ref_vcf": str(v["ref_vcf"])} for name, v in fingerprints.items()}

def set_paths(output_dir, fastq_R1=None, fastq_R2=None, fastq=None,
              generation=None, centromere=None, fingerprints=None) -> dict:
    """
    Paths block of the config for the FASTQ pipelines.
    The reference genome and chain file are always the packaged ones.
    generation = (hq_snp_positions, pickle_file, snp_info, cluster_json, representative_accessions)
    """
    hq, pickle_file, snp_info, cluster_json, representatives = generation or (None,) * 5
    return {
        "fastq_R1": str(fastq_R1),
        "fastq_R2": str(fastq_R2),
        "fastq": str(fastq),
        "output_dir": str(output_dir),
        "reference": str(TEST_REFERENCE if is_test_data(*[p for p in (fastq_R1, fastq_R2, fastq) if p]) else REFERENCE),
        "hq_positions": str(hq),
        "chain_file": str(CHAIN_FILE),
        "chromosome_mapping": str(CHROMOSOME_MAPPING),
        "pickle_file": str(pickle_file),
        "SNP_Info": str(snp_info),
        "Cluster_json": str(cluster_json),
        "representative_accessions": str(representatives),
        "centromere": str(centromere),
        "fingerprints": fingerprints_to_config(fingerprints),
    }

"""
-----------------------------
Config writing and Snakemake execution
-----------------------------
"""

def merge_config(template: dict, overrides: dict) -> dict:
    """
    Recursively merge overrides into template.
    Only assign override if the template key is missing or None.
    """
    for k, v in overrides.items():
        if k in template:
            if isinstance(v, dict) and isinstance(template[k], dict):
                merge_config(template[k], v)
            elif template[k] is None:
                template[k] = v
        else:
            template[k] = v
    return template

def load_and_write_config(template_name: str, overrides: dict, output_dir: Path, config_name: str = "config.yaml") -> Path:
    """
    Load template config, apply overrides only if the key is missing or has no value, and write to config_name.
    """
    with resources.files(rispy_cli).joinpath(f"templates/{template_name}").open("r") as f:
        config = yaml.safe_load(f)
    config = merge_config(config, overrides)
    output_dir.mkdir(exist_ok=True)
    config_path = output_dir / config_name
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    typer.echo(f"Config written to {config_path}")
    return config_path

def run_snakemake(snakefile: Path, config_path: Path, cores: int, mem_mb: int, dry_run: bool) -> None:
    """
    Run Snakemake with given parameters, logging stdout/stderr.
    """
    # Create logger
    log_dir = Path(config_path).parent
    log_file = log_dir / "RiSpy.log"
    logger = logging.getLogger("RiSpy")
    logger.setLevel(logging.INFO)
    # Avoid duplicate handlers if function is called multiple times
    if not logger.handlers:
        handler = RotatingFileHandler(
            log_file,
            maxBytes=5_000_000,   # 5 MB
            backupCount=5
        )
        formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    # Build Snakemake command
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--configfile", str(config_path),
        "--cores", str(cores),
        "--resources", f"mem_mb={mem_mb}",
        "--rerun-triggers", "mtime",
    ]
    if dry_run:
        cmd.extend(["--dry-run", "--printshellcmds"])
    logger.info(f"Running Snakemake command: {' '.join(cmd)}")
    # Execute Snakemake with full logging
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        if result.stdout:
            logger.info(f"Snakemake STDOUT:\n{result.stdout}")
        if result.stderr:
            logger.info(f"Snakemake STDERR:\n{result.stderr}")
        logger.info("Snakemake completed successfully.")
        typer.echo(f"{snakefile.name} completed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error("Snakemake failed.")
        logger.error(f"Exit code: {e.returncode}")
        logger.error(f"Command: {' '.join(cmd)}")
        if e.stdout:
            logger.error(f"STDOUT:\n{e.stdout}")
        if e.stderr:
            logger.error(f"STDERR:\n{e.stderr}")
        raise RuntimeError("Snakemake execution failed. See log file for details.")

def write_combined_snakefile(snakefiles: list, output_dir: Path) -> Path:
    """
    For --dry-run only: writes a wrapper Snakefile that loads all Snakefiles as modules,
    so Snakemake sees one DAG and can dry-run every step (separate runs cannot, because
    a dry-run does not create the files the next Snakefile needs).
    """
    lines = []
    for i, sf in enumerate(snakefiles):
        lines += [f"module m{i}:",
                  f"    snakefile: '{SNAKEFILES / sf}'",
                  "    config: config",
                  f"use rule * from m{i} as m{i}_*",
                  ""]
    lines += ["rule all:",
              "    default_target: True",
              "    input: " + ", ".join(f"rules.m{i}_all.input" for i in range(len(snakefiles)))]
    combined = Path(output_dir) / "dry_run.smk"
    combined.write_text("\n".join(lines) + "\n")
    return combined

def run_pipeline(template, overrides, output_dir, snakefiles, cores, mem_mb, dry_run, config_name="config.yaml"):
    """Writes the config, points TMPDIR next to it and runs the Snakefiles in order (or dry-runs them all together)."""
    config_path = load_and_write_config(template, overrides, output_dir, config_name)
    os.environ["TMPDIR"] = str(Path(config_path).parent / "tmp")
    report = Path(output_dir) / "report.html"
    if not dry_run and report.exists():
        report.unlink()  # a report must only exist for a run that finished without errors
    if dry_run:
        run_snakemake(write_combined_snakefile(snakefiles, Path(config_path).parent), config_path, cores, mem_mb, dry_run)
        return
    for sf in snakefiles:
        run_snakemake(SNAKEFILES / sf, config_path, cores, mem_mb, dry_run)

def is_test_data(*paths: Path) -> bool:
    """
    True if the given inputs are RiSpy's test files, i.e. they live in
    a .../rispy_cli/test_data/... folder.
    """
    def _in_test_data(p: Path) -> bool:
        parts = Path(p).resolve().parts
        return any(parts[i] == "rispy_cli" and parts[i + 1] == "test_data"
                   for i in range(len(parts) - 1))
    return all(_in_test_data(p) for p in paths)