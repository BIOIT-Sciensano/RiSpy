import yaml
from pathlib import Path
from rispy_cli.generate_html_fingerprint_detection import (
    MATCH_THRESHOLD,
    QC_NOTE,
    write_table,
    write_box,
    write_header,
    write_footer,
    write_overview,
    write_run_details,
    load_qc_metrics,
    load_detection,
    find_matches,
    identification_threshold,
    write_detection_sections,
)

GENERATION_NOTE = ("<b>How to use the fingerprint.</b> The barcode table "
                   "(<i>fingerprint/2snp_fingerprints_annotated_withranks.tsv</i>) and the reference VCF "
                   "(<i>fingerprint/barcodes.vcf</i>) together form the fingerprint of this line. Pass them to "
                   "<code>illudetect</code> or <code>ontdetect</code> with <code>--barcodes</code> and "
                   "<code>--ref-vcf</code> to test new samples for this line.")

def write_match_warning(f, matches: list, threshold: int = MATCH_THRESHOLD) -> None:
    """Writes the 'matching variety detected' warning box."""
    names = ", ".join(f"<b>{v}</b> ({c} / {t} barcodes)" for v, c, t in matches)
    write_box(f, "warning", f"""<b>&#9888; Warning: a line matching the sample was detected in the RiSpy fingerprint atlas.</b><br>
        Fingerprint detection found {threshold} or more barcodes for: {names}.<br>
        The generated fingerprint may not be able to distinguish the line in the sample
        from the matching line(s) in the database.""")

def generation_outcome(fingerprint_metrics: list) -> str:
    """One sentence summarising the generation result."""
    metrics = dict(fingerprint_metrics)
    n_barcodes = str(metrics["Number of generated 2-SNV barcodes"])
    n_snvs = metrics["Count of unique SNVs comprising the fingerprint"]
    if not n_barcodes.isdigit():
        return ("No fingerprint could be generated: no combination of two SNVs distinguishes this line "
                "from all lines in the RiSpy fingerprint atlas.")
    return f"A fingerprint of <b>{n_barcodes}</b> 2-SNV barcodes, built from <b>{n_snvs}</b> unique SNVs, was generated for this line."

def generate_report_html(html_file: str, config_qc: Path | None, config_fingerprint: Path,
                         detection_result: Path | None = None, top_n: int = 10, run_config: dict | None = None) -> None:
    """
    Generate one HTML report with QC, fingerprint detection and fingerprint generation.
    detection_result: output of the detection run before generation (.yaml or .tsv);
                      None when detection was skipped (e.g. test data) -> no detection section.
    top_n: number of top-matching lines to show (multi-fingerprint detection).
    run_config: the Snakemake config (for overview, settings, output files and software versions).
    """
    threshold = identification_threshold(run_config)
    qc_metrics = load_qc_metrics(config_qc)
    if detection_result is not None:
        detection_metrics, ranking_table = load_detection(detection_result)
    else:
        detection_metrics, ranking_table = [], []
    matches = find_matches(detection_metrics, ranking_table, threshold)
    fingerprint_metrics = []
    with open(config_fingerprint, "r") as f:
        fp = yaml.safe_load(f)
        fingerprint_metrics.append(("Number of generated 2-SNV barcodes", fp["Number of generated 2-SNV barcodes"]))
        fingerprint_metrics.append(("Count of unique SNVs comprising the fingerprint", fp["Count of unique SNVs comprising the fingerprint"]))
    with open(html_file, "w") as f:
        write_header(f,
                     "Fingerprint Detection and Generation Report" if detection_result is not None else "Fingerprint Generation Report")
        if matches:
            write_match_warning(f, matches, threshold)
        write_box(f, "summary", "<b>Result summary.</b> " + generation_outcome(fingerprint_metrics))
        if run_config:
            write_overview(f, run_config)
        if qc_metrics:
            write_table(f, "QC Report", qc_metrics)
            write_box(f, "note", QC_NOTE)
        write_detection_sections(f, detection_metrics, ranking_table, top_n, threshold)
        write_table(f, "Fingerprint Generation Results", fingerprint_metrics)
        write_box(f, "note", GENERATION_NOTE)
        if run_config:
            write_run_details(f, run_config)
        write_footer(f)