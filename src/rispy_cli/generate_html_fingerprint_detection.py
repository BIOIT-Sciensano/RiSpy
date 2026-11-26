import base64
import platform
import re
import shutil
import subprocess
import importlib.resources as resources
from importlib.metadata import version as package_version, PackageNotFoundError
from datetime import datetime
from pathlib import Path
import yaml
from rispy_cli.io import build_fq_by_sample

MATCH_THRESHOLD = 43  # default identification threshold

CITATION = ("Zolfaghari, A., Fraiture, M-A., Vanneste, K., Stuyts, A., Frouin, J., Meunier, A-C., Deforce, D., "
            "Roosens, N. H. C., &amp; D'aes, J. (2026). RiSpy: A feature selection-based fingerprinting framework "
            "for accurate identification of genome-edited rice lines. <i>Briefings in Bioinformatics</i>, 27(4), bbag406.")
DOI = "https://doi.org/10.1093/bib/bbag406"

COMMANDS = {
    "illugen": "Fingerprint generation from short-read Illumina paired-end WGS data (preceded by detection against the RiSpy fingerprint atlas)",
    "ontgen": "Fingerprint generation from long-read ONT WGS data (preceded by detection against the RiSpy fingerprint atlas)",
    "illudetect": "Fingerprint detection from short-read Illumina paired-end WGS data",
    "ontdetect": "Fingerprint detection from long-read ONT WGS data",
}

STYLE = """
        body { font-family: Arial, sans-serif; margin: 30px; background-color: #f9f9f9; }
        h1 { text-align: center; color: #2c3e50; margin: 0; padding-top: 20px; }
        h2 { color: #34495e; margin-top: 30px; text-align: center; }
        table { width: 60%; margin: 20px auto; border-collapse: collapse; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        th, td { padding: 12px 15px; text-align: left; }
        td { word-break: break-word; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #ecf0f1; }
        tr:hover { background-color: #d1ecf1; }
        .metric { font-weight: bold; color: #2c3e50; white-space: nowrap; }
        .header { position: relative; min-height: 80px; }
        .logo { position: absolute; left: 0; top: 0; height: 70px; }
        .subtitle { text-align: center; font-style: italic; color: #7f8c8d; margin-top: 8px; }
        .summary { width: 60%; margin: 20px auto; padding: 15px 20px; background-color: #e8f4fd;
                   border-left: 6px solid #3498db; color: #1b3a57; line-height: 1.5; }
        .warning { width: 60%; margin: 20px auto; padding: 15px 20px; background-color: #fff3cd;
                   border-left: 6px solid #e67e22; color: #7a4b00; line-height: 1.5; }
        .note { width: 60%; margin: 10px auto; color: #555; font-size: 0.9em; line-height: 1.5; }
        .citation { width: 60%; margin: 40px auto 20px; padding: 15px 20px; border-top: 2px solid #3498db;
                    color: #2c3e50; font-size: 0.9em; line-height: 1.5; }
"""

QC_NOTE = ("<b>How to interpret.</b> Mean and median depth: number of reads covering each reference position. "
           "Breadth of coverage: percentage of the reference covered at or above the depth cutoff (see Analysis Settings). "
           "Fold80: evenness of coverage (1 = perfectly even). Low depth or breadth reduces the number of barcodes "
           "that can be detected.")

DETECTION_NOTE = ("<b>How to interpret.</b> Each fingerprint consists of a set of 2-SNV barcodes. "
                  "A barcode is counted as detected only when both of its SNVs are confirmed in the sample "
                  "(alternative allele called, or reference allele supported by at least 5 reads). "
                  "A line is considered identified when at least {threshold} of its barcodes are detected. "
                  "Low depth or breadth of coverage reduces the number of "
                  "detectable barcodes, so check the QC Report before concluding that a line is absent.")

TOP_LINES_NOTE = (" Lines with the same number of detected barcodes share one place in the top list, so more lines "
                  "than requested can be listed.")

# Third-party tools: name -> (executable, conda package, arguments that print the version)
TOOLS = {
    "FastQC": ("fastqc", "fastqc", ["--version"]),
    "fastp": ("fastp", "fastp", ["--version"]),
    "BWA": ("bwa", "bwa", []),
    "minimap2": ("minimap2", "minimap2", ["--version"]),
    "NanoPlot": ("NanoPlot", "nanoplot", ["--version"]),
    "SeqKit": ("seqkit", "seqkit", ["version"]),
    "Chopper": ("chopper", "chopper", ["--version"]),
    "SAMtools": ("samtools", "samtools", ["--version"]),
    "Qualimap": ("qualimap", "qualimap", ["--help"]),
    "BCFtools": ("bcftools", "bcftools", ["--version"]),
    "BEDTools": ("bedtools", "bedtools", ["--version"]),
    "liftOver": ("liftOver", "ucsc-liftover", None),
    "CD-HIT": ("cd-hit-est-2d", "cd-hit", ["-h"]),
}
ILLUMINA_TOOLS = ["FastQC", "fastp", "BWA"]
ONT_TOOLS = ["NanoPlot", "SeqKit", "Chopper", "minimap2"]
COMMON_TOOLS = ["SAMtools", "Qualimap", "BCFtools", "BEDTools", "liftOver"]
GENERATION_TOOLS = ["CD-HIT"]
PYTHON_PACKAGES = ["snakemake", "pandas", "numpy", "scikit-learn", "pysam"]

"""
-----------------------------
Page building blocks
-----------------------------
"""

def safe_version(package: str) -> str:
    try:
        return package_version(package)
    except PackageNotFoundError:
        return "unknown"

def logo_html() -> str:
    """Sciensano logo embedded as base64 (empty string if the file is missing)."""
    try:
        data = (resources.files("rispy_cli.data") / "sciensano_logo.jpg").read_bytes()
    except (FileNotFoundError, OSError):
        return ""
    mime = "image/png" if data[:8] == b"\x89PNG\r\n\x1a\n" else "image/jpeg"
    return f'<img class="logo" src="data:{mime};base64,{base64.b64encode(data).decode()}" alt="Sciensano">'

def write_table(f, title: str, metrics: list, col1: str = "Metric", col2: str = "Value") -> None:
    """
    Writes a single table to an HTML file.
    f: file object
    title: string, table title
    metrics: list of (key, value) tuples
    col1, col2: optional column headers
    """
    f.write(f"<h2>{title}</h2>\n")
    f.write("<table>\n")
    f.write(f"<tr><th>{col1}</th><th>{col2}</th></tr>\n")
    for key, value in metrics:
        f.write(f"<tr><td class='metric'>{key}</td><td>{value}</td></tr>\n")
    f.write("</table>\n\n")

def write_box(f, css_class: str, html: str) -> None:
    """Writes a text box (summary / warning / note)."""
    f.write(f'<div class="{css_class}">{html}</div>\n')

def write_header(f, title: str = "Fingerprint Detection Report") -> None:
    """Writes the page head, logo (top left), tool name, version and generation timestamp."""
    report_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    f.write(f"""<!DOCTYPE html>
        <html lang="en">
        <head>
        <meta charset="UTF-8">
        <title>RiSpy {title}</title>
        <style>{STYLE}</style>
        </head>
        <body>
        <div class="header">{logo_html()}<h1>RiSpy &ndash; {title}</h1></div>
        <p class="subtitle">RiSpy version {safe_version("RiSpy")} &middot; Generated on: {report_datetime}</p>
        """)

def write_footer(f) -> None:
    """Writes the citation box and closes the page."""
    f.write(f"""<div class="citation">
        <b>How to cite.</b> If you use RiSpy in your work, please cite:<br>
        {CITATION} {DOI}
        </div>
        </body></html>""")

"""
-----------------------------
What was run (read from the config itself)
-----------------------------
"""

def _is_set(value) -> bool:
    return value not in (None, "None", "null", "")

def is_ont(run_config: dict) -> bool:
    return _is_set(run_config["paths"].get("fastq"))

def is_generation(run_config: dict) -> bool:
    return _is_set(run_config["paths"].get("hq_positions"))

def command_name(run_config: dict) -> str:
    command = (run_config.get("run_info") or {}).get("command")
    if command in COMMANDS:
        return command
    return ("ont" if is_ont(run_config) else "illu") + ("gen" if is_generation(run_config) else "detect")

def sample_name(run_config: dict) -> str:
    paths = run_config["paths"]
    return next(iter(build_fq_by_sample(fastq_R1=paths.get("fastq_R1"), fastq_R2=paths.get("fastq_R2"), fastq=paths.get("fastq"))))

"""
-----------------------------
Results
-----------------------------
"""

def load_qc_metrics(config_qc: Path | None) -> list:
    """QC metrics as [(key, value), ...]; empty list when config_qc is None."""
    if config_qc is None:
        return []
    with open(config_qc, "r") as f:
        qc = yaml.safe_load(f)
    return [
        ("Mean Depth of Coverage", f"{qc['Mean Depth of Coverage']}x"),
        ("Median Depth of Coverage", f"{qc['Median Depth of Coverage']}x"),
        ("Breadth of Coverage", f"{qc['Breadth of Coverage']} %"),
        ("Number of Primary Mapped Reads", f"{int(qc['Number of Primary Mapped Reads']):,}"),
        ("Fold80", qc["Fold80"]),
    ]

def load_detection(config_fingerprint: Path) -> tuple[list, list]:
    """
    Reads the detection output.
    Returns (fingerprint_metrics, ranking_table):
      - single-fingerprint mode (.yaml): metrics = [("Fingerprint", name), ("Number of Detected Barcodes", "93 / 100")]
      - multi-fingerprint mode (.tsv):   ranking = [(variety, detected, total), ...] sorted by detected
    """
    config_fingerprint = Path(config_fingerprint)
    fingerprint_metrics = []
    ranking_table = []
    if config_fingerprint.suffix == ".yaml":
        with open(config_fingerprint, "r") as f:
            fp = yaml.safe_load(f)
        total = fp.get("Number of Barcodes in Fingerprint", "?")
        fingerprint_metrics.append(("Fingerprint", "Supplied fingerprint"))
        fingerprint_metrics.append(("Number of Detected Barcodes", f"{fp['Number of Detected Barcodes']} / {total}"))
    else:
        with open(config_fingerprint, "r") as f:
            next(f)  # skip header
            for line in f:
                fp, count, total = line.strip().split("\t")
                ranking_table.append((fp, int(count), int(total)))
        ranking_table = sorted(ranking_table, key=lambda x: x[1], reverse=True)
    return fingerprint_metrics, ranking_table

def find_matches(fingerprint_metrics: list, ranking_table: list, threshold: int = MATCH_THRESHOLD) -> list:
    """Lines with at least `threshold` detected barcodes, as [(variety, detected, total), ...], highest first."""
    if ranking_table:
        return [row for row in ranking_table if row[1] >= threshold]
    metrics = dict(fingerprint_metrics)
    value = metrics.get("Number of Detected Barcodes")
    if value is not None:
        count, total = [v.strip() for v in str(value).split("/")]
        if int(count) >= threshold:
            return [(metrics.get("Fingerprint", "the supplied fingerprint"), int(count), total)]
    return []

def identification_threshold(run_config: dict | None) -> int:
    """Identification threshold of this run (config), MATCH_THRESHOLD if not set."""
    value = ((run_config or {}).get("parameters") or {}).get("identification_threshold")
    return int(value) if _is_set(value) else MATCH_THRESHOLD

def top_rows(ranking_table: list, top_n: int) -> list:
    """
    Lines with at least 1 detected barcode whose count is among the top_n highest counts.
    Lines with the same count share one place, so more than top_n rows can be returned.
    """
    top_counts = set(sorted({c for _, c, _ in ranking_table if c > 0}, reverse=True)[:top_n])
    return [row for row in ranking_table if row[1] in top_counts]

def names_with_counts(rows: list) -> str:
    return ", ".join(f"<b>{v}</b> ({c} / {t} barcodes detected)" for v, c, t in rows)

def tie_note(rows: list) -> str:
    """Sentence naming the lines that share the highest number of detected barcodes (empty if no tie)."""
    top = rows[0][1]
    tied = [v for v, c, _ in rows if c == top]
    if len(tied) < 2:
        return ""
    return f" {', '.join(f'<b>{v}</b>' for v in tied)} have the same number of detected barcodes ({top})."

def detection_outcome(fingerprint_metrics: list, ranking_table: list, threshold: int = MATCH_THRESHOLD) -> str:
    """One or two sentences summarising the detection result."""
    matches = find_matches(fingerprint_metrics, ranking_table, threshold)
    if ranking_table:
        if len(matches) == 1:
            return (f"The sample was identified as {names_with_counts(matches)}, the only line in the RiSpy fingerprint atlas "
                    f"that reached the identification threshold of {threshold} detected barcodes.")
        if matches:
            return (f"More than one line reached the identification threshold of {threshold} detected barcodes: "
                    f"{names_with_counts(matches)}.{tie_note(matches)} This can indicate closely related lines "
                    "or a sample containing more than one line.")
        best = ranking_table[0][1]
        if best == 0:
            return ("No barcode of any line in the RiSpy database was detected in the sample, so the sample does not "
                    "match any line in the fingerprint atlas. Check the QC Report: very low depth or breadth of coverage "
                    "can also cause this.")
        closest = [row for row in ranking_table if row[1] == best]
        if len(closest) == 1:
            return (f"No line in the RiSpy database reached the identification threshold of {threshold} detected "
                    f"barcodes. The closest line is {names_with_counts(closest)}.")
        return (f"No line in the RiSpy database reached the identification threshold of {threshold} detected "
                f"barcodes. The closest lines, all with the same number of detected barcodes, are {names_with_counts(closest)}.")
    if fingerprint_metrics:
        value = dict(fingerprint_metrics)["Number of Detected Barcodes"]
        if matches:
            return (f"The supplied fingerprint was detected in the sample: {value} barcodes detected "
                    f"(identification threshold: {threshold}).")
        return (f"The supplied fingerprint was not detected in the sample: {value} barcodes detected, "
                f"below the identification threshold of {threshold}.")

def write_detection_sections(f, fingerprint_metrics: list, ranking_table: list, top_n: int, threshold: int = MATCH_THRESHOLD) -> None:
    """Writes the detection table (single or multi mode) and the interpretation note."""
    if fingerprint_metrics:
        metrics = dict(fingerprint_metrics)
        write_table(f, "Fingerprint Detection Results",
                    [(metrics.get("Fingerprint", "Supplied fingerprint"), metrics["Number of Detected Barcodes"])],
                    col1="Fingerprint", col2="Detected / Total Barcodes")
    rows = top_rows(ranking_table, top_n)
    if rows:
        write_table(f, f"Top {top_n} Matching Lines", [(v, f"{c} / {t}") for v, c, t in rows],
                    col1="Line", col2="Detected / Total Barcodes")
    if fingerprint_metrics or rows:
        write_box(f, "note", DETECTION_NOTE.format(threshold=threshold) + (TOP_LINES_NOTE if rows else ""))

"""
-----------------------------
Run details (overview, settings, outputs, software)
-----------------------------
"""

def overview_rows(run_config: dict) -> list:
    """Which analysis was run, on which data, and with which options."""
    info = run_config.get("run_info") or {}
    paths = run_config["paths"]
    command = command_name(run_config)
    inputs = [p for p in (paths.get("fastq_R1"), paths.get("fastq_R2"), paths.get("fastq")) if _is_set(p)]
    rows = [
        ("Analysis", f"<b>{command}</b>: {COMMANDS[command]}"),
        ("Sample", sample_name(run_config)),
        ("Input data", "<br>".join(inputs)),
        ("Output directory", paths["output_dir"]),
        ("Reference genome", paths["reference"]),
    ]
    fingerprints = paths.get("fingerprints") or {}
    if len(fingerprints) > 1:
        rows.append(("Fingerprints tested", f"All {len(fingerprints)} lines in the RiSpy fingerprint atlas"))
    elif fingerprints:
        fp = next(iter(fingerprints.values()))
        rows.append(("Fingerprint tested", f"{fp['barcodes']}<br>{fp['ref_vcf']}"))
    options = info.get("non_default_parameters") or {}
    rows.append(("Options set by the user",
                 "<br>".join(f"{k} {v}".strip() for k, v in options.items()) or "None (all defaults used)"))
    return rows

def write_overview(f, run_config: dict) -> None:
    write_table(f, "Analysis Overview", overview_rows(run_config), col1="Item", col2="Details")

def settings_rows(run_config: dict) -> list:
    """All analysis parameters from the config, flattened as section.key (top_n left out when one fingerprint is tested)."""
    single = len(run_config["paths"].get("fingerprints") or {}) == 1
    rows = []
    def flatten(d, prefix=""):
        for k, v in d.items():
            if isinstance(v, dict):
                flatten(v, f"{prefix}{k}.")
            elif not (single and k == "top_n"):
                rows.append((f"{prefix}{k}", v))
    flatten(run_config.get("parameters") or {})
    return rows

def tool_version(executable: str, conda_package: str, version_args: list | None) -> str:
    """
    Version of a third-party tool: first from the conda metadata of the environment the
    executable lives in, otherwise by running the tool and taking the first version number it prints.
    """
    path = shutil.which(executable)
    if path is None:
        return "not found"
    for meta in (Path(path).parent.parent / "conda-meta").glob(f"{conda_package}-[0-9]*.json"):
        return meta.name[len(conda_package) + 1:].split("-")[0]
    if version_args is None:
        return "unknown"
    try:
        out = subprocess.run([executable] + version_args, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    m = re.search(r"\d+\.\d+(?:\.\d+)?", out.stdout + out.stderr)
    return m.group(0) if m else "unknown"

def software_rows(run_config: dict) -> list:
    """Versions of the software actually used by this run (Illumina or ONT; generation adds CD-HIT)."""
    names = (ONT_TOOLS if is_ont(run_config) else ILLUMINA_TOOLS) + COMMON_TOOLS + (GENERATION_TOOLS if is_generation(run_config) else [])
    rows = [(name, tool_version(*TOOLS[name])) for name in names]
    checkpoint = str(((run_config.get("parameters") or {}).get("call_variants") or {}).get("CHECKPOINT", ""))
    m = re.search(r"deepvariant/(\d+\.\d+\.\d+)/", checkpoint)
    rows.append(("DeepVariant", m.group(1) if m else "user-supplied checkpoint (see Analysis Settings)"))
    rows.append(("Python", platform.python_version()))
    rows += [(package, safe_version(package)) for package in PYTHON_PACKAGES]
    return rows

def write_run_details(f, run_config: dict) -> None:
    """Settings and software versions (everything needed to repeat the analysis)."""
    write_table(f, "Analysis Settings", settings_rows(run_config), col1="Parameter", col2="Value")
    write_box(f, "note", "To repeat this analysis, run the same command with the options listed in the "
                         "Analysis Overview.")
    write_table(f, "Software Versions", software_rows(run_config), col1="Software", col2="Version")

"""
-----------------------------
Report
-----------------------------
"""

def generate_report_html(html_file: str, config_qc: Path | None, config_fingerprint: Path, top_n: int, run_config: dict | None = None) -> None:
    """
    Generate the fingerprint detection HTML report.
    run_config: the Snakemake config (for overview, settings, software versions and the identification threshold).
    """
    threshold = identification_threshold(run_config)
    qc_metrics = load_qc_metrics(config_qc)
    fingerprint_metrics, ranking_table = load_detection(config_fingerprint)
    with open(html_file, "w") as f:
        write_header(f, "Fingerprint Detection Report")
        write_box(f, "summary", "<b>Result summary.</b> " + detection_outcome(fingerprint_metrics, ranking_table, threshold))
        if run_config:
            write_overview(f, run_config)
        if qc_metrics:
            write_table(f, "QC Report", qc_metrics)
            write_box(f, "note", QC_NOTE)
        write_detection_sections(f, fingerprint_metrics, ranking_table, top_n, threshold)
        if run_config:
            write_run_details(f, run_config)
        write_footer(f)