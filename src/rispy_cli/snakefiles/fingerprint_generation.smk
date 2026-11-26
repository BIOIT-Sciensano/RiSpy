from pathlib import Path
from rispy_cli.fingerprint_generation_pipeline import (
    hq_positions,
    chromosome_name_conversion,
    genotype_detection,
    pickle_to_fasta,
    update_cluster_json,
    feature_selection,
    add_info_to_indices,
    list_remaining_accessions,
    screen_2_snp_combinations,
    add_info_to_snp_combinations,
    unique_SNPs,
    bed_maker,
    vcf_maker
)
from rispy_cli.io import build_fq_by_sample
from rispy_cli.generate_html_fingerprint_generation import generate_report_html

fastq_R1 = config["paths"].get("fastq_R1")
fastq_R2 = config["paths"].get("fastq_R2")
fastq    = config["paths"].get("fastq")

fq_by_sample = build_fq_by_sample(fastq_R1=fastq_R1,fastq_R2=fastq_R2,fastq=fastq)

root = Path(config["paths"]["output_dir"])

HQ_positions = config["paths"]["hq_positions"]

chain_file_old_nip_to_new_nip_ref = config["paths"]["chain_file"]

chromosome_mapping_file = config["paths"]["chromosome_mapping"]

pickle_file = config["paths"]["pickle_file"]

SNP_Info = config["paths"]["SNP_Info"]

Cluster_json = config["paths"]["Cluster_json"]

representative_accessions = config["paths"]["representative_accessions"]

filtered_vcf = config["paths"].get("filtered_vcf", None)

unfiltered_vcf = config["paths"].get("unfiltered_vcf", None)

bam = config["paths"].get("bam", None)

ref = config["paths"].get("reference", None)

report_html_file = root / 'report.html'

rule all:
    input:
        VCF = expand(str(root / '{sample}' / 'fingerprint' / 'barcodes.vcf'),sample=fq_by_sample.keys()),
        HTML = report_html_file
    resources: tmpdir = str(root / 'tmp')

rule hq_positions:
    """
    Making a bed file out of 160K SNPs positions
    """
    input:
        HQ_POSITIONS = HQ_positions
    output:
        HQ_BED = root / '{sample}' / 'fingerprint' / '160K_positions_with_ids.bed'
    threads: config["resources"]["hq_positions"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["hq_positions"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'hq_positions.txt'
    run:
        hq_positions(input.HQ_POSITIONS,output.HQ_BED)

rule liftover:
    """
    Executing liftOver using the annotated BED file, to get the positions in the new reference genome
    """
    input:
        HQ_BED = rules.hq_positions.output.HQ_BED,
        CHAINFILE = chain_file_old_nip_to_new_nip_ref
    output:
        HQ_BED_LIFTOVER = root / '{sample}' / 'fingerprint' / '160K_positions_liftOver.bed',
        UNMAPPED = root / '{sample}' / 'fingerprint' / 'unmapped.bed'
    threads: config["resources"]["liftover"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["liftover"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'liftover.txt'
    shell:
        """
        liftOver {input.HQ_BED} {input.CHAINFILE} {output.HQ_BED_LIFTOVER} {output.UNMAPPED}
        """

rule chromosome_name_conversion:
    """
    Converting the chromosome names in the liftOver output to match the new reference genome
    """
    input:
        CHROM_MAP = chromosome_mapping_file,
        HQ_BED_LIFTOVER = rules.liftover.output.HQ_BED_LIFTOVER
    output:
        HQ_BED_CONVERT_CHR = root / '{sample}' / 'fingerprint' / '160K_positions_liftOver_converted_chr_names.bed'
    threads: config["resources"]["chromosome_name_conversion"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["chromosome_name_conversion"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'chromosome_name_conversion.txt'
    run:
        chromosome_name_conversion(input.CHROM_MAP, input.HQ_BED_LIFTOVER, output.HQ_BED_CONVERT_CHR)

rule hq_variants:
    """
    Checking the variants in the filtered VCF file that are in the 160K positions
    """
    input:
        HQ_BED_CONVERT_CHR = rules.chromosome_name_conversion.output.HQ_BED_CONVERT_CHR,
        VCF_FILE = lambda wildcards: filtered_vcf or root / wildcards.sample / 'deep_variant' / 'Filtering' / 'filtered_variants.vcf.gz',
        No_Filter_VCF = lambda wildcards: unfiltered_vcf or root  / wildcards.sample / 'deep_variant' / 'variants.vcf.gz',
    output:
        HQ_VARIANTS = root / '{sample}' / 'fingerprint' / 'hg_variants_in_160K_positions.txt',
        REF_POSITIONS = root / '{sample}' / 'fingerprint' / 'reference_alleles.bed',
        Low_Confidence_Ref_Alleles = root / '{sample}' / 'fingerprint' / 'Low_Confidence_Ref_Alleles.txt',
    threads: config["resources"]["hq_variants"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["hq_variants"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'hq_variants.txt'
    shell:
        """
        bedtools intersect -a {input.HQ_BED_CONVERT_CHR} -b {input.VCF_FILE} -wa -wb > {output.HQ_VARIANTS};
        bedtools intersect -a {input.HQ_BED_CONVERT_CHR} -b {input.VCF_FILE} -v > {output.REF_POSITIONS};
        bedtools intersect -a {output.REF_POSITIONS} -b {input.No_Filter_VCF} -wa -wb > {output.Low_Confidence_Ref_Alleles};
        """

rule genotype_detection:
    """
    Detecting the genotype of our variety of interest in the 160K positions
    """
    input:
        HQ_VARIANTS = rules.hq_variants.output.HQ_VARIANTS,
        PICKLEFILE = pickle_file,
        UNMAPPED = rules.liftover.output.UNMAPPED,
        Low_Confidence_Ref_Alleles = rules.hq_variants.output.Low_Confidence_Ref_Alleles,
        REF_POSITIONS = rules.hq_variants.output.REF_POSITIONS,
        BAM = lambda wildcards: bam or root / wildcards.sample / 'mapping' / 'mapped_reads_sorted.bam',
        Fasta_path = ref or root / 'reference' / 'reference.fasta',
        Mapping_tsv = root / '{sample}' / 'fingerprint' / '160K_positions_liftOver_converted_chr_names.bed'
    output:
        PICKLE_OUT = root / '{sample}' / 'fingerprint' / 'hq_genotype_included_variety_of_interest.pkl'
    threads: config["resources"]["genotype_detection"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["genotype_detection"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'genotype_detection.txt'
    run:
        genotype_detection(input.HQ_VARIANTS, input.PICKLEFILE, input.UNMAPPED, input.Low_Confidence_Ref_Alleles, input.REF_POSITIONS, input.BAM, input.Mapping_tsv, input.Fasta_path, output.PICKLE_OUT)

rule pickle_to_fasta:
    """
    Generating a fasta file for all cultivars of 3KRG including our cultivar of interest
    """
    input:
        picklefile1 = rules.genotype_detection.output.PICKLE_OUT,
    output:
        fasta = root / '{sample}' / 'fingerprint' / 'VOI.fasta',
    threads: config["resources"]["pickle_to_fasta"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["pickle_to_fasta"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'pickle_to_fasta.txt'
    run:
        pickle_to_fasta(input.picklefile1, output.fasta)

rule clustering_cd_hit:
    """
    Finding the cluster of the variety of interest using cd-hit
    """
    input:
        voi_fasta = rules.pickle_to_fasta.output.fasta,
        rep_fasta = representative_accessions
    output:
        fasta = root / '{sample}' / 'fingerprint' / 'cultivars_90.fasta.clstr'
    params:
        output_base_name = lambda wildcards: root / wildcards.sample / 'fingerprint' / 'cultivars_90.fasta',
        clustering_identity = config["parameters"]["clustering_cd_hit"]["clustering_identity"]
    threads: config["resources"]["clustering_cd_hit"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["clustering_cd_hit"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'clustering_cd_hit.txt'
    shell:
        """
        cd-hit-est-2d \
        -i {input.rep_fasta} \
        -i2 {input.voi_fasta} \
        -o {params.output_base_name} \
        -c {params.clustering_identity} \
        -T {threads} \
        -M {resources.mem_mb}
        """

rule update_cluster_json:
    """
    Updating the cluster json file with the cluster of the variety of interest
    """
    input:
        clusters_cd_hit = rules.clustering_cd_hit.output.fasta,
        Cluster_json = Cluster_json
    output:
        json = root / '{sample}' / 'fingerprint' / 'clusters_updated.json',
    threads: config["resources"]["update_cluster_json"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["update_cluster_json"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'update_cluster_json.txt'
    run:
        update_cluster_json(input.clusters_cd_hit, input.Cluster_json, output.json)

rule feature_selection:
    """
    Selecting the best 2000 SNPs for discriminating our variety of interest from the rest of the varieties
    """
    input:
        picklefile1 = rules.genotype_detection.output.PICKLE_OUT,
        cluster_json = rules.update_cluster_json.output.json
    output:
        feature_scores = root / '{sample}' / 'fingerprint' / 'fs_scores_chi2_2.txt',
        picklefile2 = root / '{sample}' / 'fingerprint' / 'hq_genotype_preprocessed.pkl',
    threads: config["resources"]["feature_selection"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["feature_selection"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'feature_selection.txt'
    run:
        feature_selection(input.picklefile1, input.cluster_json, output.feature_scores, output.picklefile2)

rule add_info_to_indices:
    """
    Adding info (chr, pos, ref) to selected features
    """
    input:
        Top_features = rules.feature_selection.output.feature_scores,
        SNP_INFO = SNP_Info,
    output:
        Features_with_annotation = root / '{sample}' / 'fingerprint' / 'features_with_annotation.tsv',
    threads: config["resources"]["add_info_to_indices"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["add_info_to_indices"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'add_info_to_indices.txt'
    run:
        add_info_to_indices(input.Top_features, input.SNP_INFO, output.Features_with_annotation)

rule list_remaining_accessions:
    """
    Adding the number of varieties with the same genotype as the variety of interest
    """
    input:
        picklefile2 = rules.feature_selection.output.picklefile2,
        Features_with_annotation = rules.add_info_to_indices.output.Features_with_annotation
    output:
        Features_with_annotation_withremainingaccession = root / '{sample}' / 'fingerprint' / 'features_with_annotation_withremainingaccession.tsv',
    threads: config["resources"]["list_remaining_accessions"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["list_remaining_accessions"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'list_remaining_accessions.txt'
    run:
        list_remaining_accessions(input.picklefile2, input.Features_with_annotation, output.Features_with_annotation_withremainingaccession)


rule screen_2_snp_combinations:
    """
    Exhaustive screening of all combinations of 2 snps of the list with selected features
    """
    input:
        picklefile = rules.feature_selection.output.picklefile2,
        Features_with_annotation = rules.add_info_to_indices.output.Features_with_annotation,
        check = rules.list_remaining_accessions.output.Features_with_annotation_withremainingaccession
    output:
        Fingerprints = root / '{sample}' / 'fingerprint' / '2snp_fingerprints.tsv'
    params:
        max_fing = config["parameters"]["screen_2_snp_combinations"]["max_fingerprints"]
    threads: config["resources"]["screen_2_snp_combinations"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["screen_2_snp_combinations"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'screen_2_snp_combinations.txt'
    run:
        screen_2_snp_combinations(input.picklefile, input.Features_with_annotation, params.max_fing, output.Fingerprints)

rule add_info_to_snp_combinations:
    """
    Adding the annotation (chr, pos, genotype) to the snp-combinations
    """
    input:
        Combinations = rules.screen_2_snp_combinations.output.Fingerprints,
        Scores = rules.feature_selection.output.feature_scores,
        SNP_INFO = SNP_Info,
        picklefile = rules.feature_selection.output.picklefile2,
    output:
        Fingerprints_annotated_withranks = root / '{sample}' / 'fingerprint' / '2snp_fingerprints_annotated_withranks.tsv'
    params:
        config = root / 'config_results.yaml'
    threads: config["resources"]["add_info_to_snp_combinations"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["add_info_to_snp_combinations"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'add_info_to_snp_combinations.txt'
    run:
        add_info_to_snp_combinations(input.Combinations, input.Scores, input.SNP_INFO, input.picklefile, params.config, output.Fingerprints_annotated_withranks)

rule unique_SNPs:
    """
    Detecting the unique SNPs in the combinations of 2 SNPs fingerprints
    """
    input:
        Fingerprints = rules.add_info_to_snp_combinations.output.Fingerprints_annotated_withranks,
        Features_with_annotation_withremainingaccession = rules.list_remaining_accessions.output.Features_with_annotation_withremainingaccession,
    output:
        Unique_SNPs = root / '{sample}' / 'fingerprint' / 'unique_snps.tsv'
    params:
        config = root / 'config_results.yaml'
    threads: config["resources"]["unique_SNPs"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["unique_SNPs"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'unique_SNPs.txt'
    run:
        unique_SNPs(input.Fingerprints, input.Features_with_annotation_withremainingaccession, params.config, output.Unique_SNPs)

rule bed_maker:
    """
    Convert unique SNPs TSV to BED for liftOver
    """
    input:
        unqiue_snps = rules.unique_SNPs.output.Unique_SNPs,
    output:
        alt_unique_snps = root / '{sample}' / 'fingerprint' / 'alt_unique_snps.bed'
    threads: config["resources"]["bed_maker"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["bed_maker"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'bed_maker.txt'
    run:
        bed_maker(input.unqiue_snps, output.alt_unique_snps)


rule liftover2:
    """
    Executing liftOver using the BED file, to get the positions in the new reference genome for generating the vcf file which is needed for fingerprint detection
    """
    input:
        FINGER_BED = rules.bed_maker.output.alt_unique_snps,
        CHAINFILE = chain_file_old_nip_to_new_nip_ref
    output:
        FINGER_BED_LIFTOVER = root / '{sample}' / 'fingerprint' / 'fingerprints_liftover.bed',
        UNMAPPED = root / '{sample}' / 'fingerprint' / 'unmapped2.bed'
    threads: config["resources"]["liftover"]["threads"]
    resources:
        tmpdir = str(root / 'tmp'),
        mem_mb = config["resources"]["liftover"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'liftover2.txt'
    shell:
        """
        liftOver {input.FINGER_BED} {input.CHAINFILE} {output.FINGER_BED_LIFTOVER} {output.UNMAPPED}
        """

rule vcf_maker:
    """
    Generating vcf file for the generated fingerprint
    """
    input:
        FINGER_BED_LIFTOVER = rules.liftover2.output.FINGER_BED_LIFTOVER,
        Fasta_path = ref or root / 'reference' / 'reference.fasta'
    output:
        vcf = root / '{sample}' / 'fingerprint' / 'barcodes.vcf',
    threads: config["resources"]["vcf_maker"]["threads"]
    resources:
        tmpdir= str(root / 'tmp'),
        mem_mb= config["resources"]["vcf_maker"]["memory"]
    benchmark: root / '{sample}' / 'fingerprint' / 'benchmarks' / 'vcf_maker.txt'
    run:
        vcf_maker(input.FINGER_BED_LIFTOVER, input.Fasta_path, output.vcf)

rule report_html:
    """
    Generating the final HTML report
    """
    input:
        VCF = expand(str(root / '{sample}' / 'fingerprint' / 'barcodes.vcf'), sample=fq_by_sample.keys()) # To ensure this rule runs at the end of the pipeline
    output:
        HTML = report_html_file
    params:
        config_fing = root / 'config_results.yaml',
        config_qc = next((root / s / "config_results.yaml" for s in fq_by_sample.keys() if (root / s / "config_results.yaml").exists()), None),
        detection_result = next((root / s / 'Fingerprint_Detection' / 'fingerprint_summary.tsv' for s in fq_by_sample.keys() if (root / s / 'Fingerprint_Detection' / 'fingerprint_summary.tsv').exists()),None)
    threads: config["resources"]["report_html"]["threads"]
    resources:
        mem_mb = config["resources"]["report_html"]["memory"]
    run:
        generate_report_html(output.HTML,params.config_qc,params.config_fing,params.detection_result,run_config=config)