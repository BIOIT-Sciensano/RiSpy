# RiSpy

### Fingerprint Generation & Detection Workflow

A workflow for **genetic fingerprint generation/detection** using Whole-Genome Sequencing (WGS) data from **Illumina** short-read and **Oxford Nanopore Technology (ONT)** long-read platforms.    

RiSpy is a fingerprint‑based identification tool for rice lines. It can detect the identity of any rice line using previously generated fingerprints stored in the RiSpy repository,  
and it can also generate new fingerprints for rice lines whose profiles do not yet exist in the fingerprint atlas. RiSpy supports fingerprinting in both gene‑edited (GE) and non‑GE rice lines,  
enabling reliable identification, verification, and tracking across breeding, research, and genome‑editing workflows.

- A **fingerprint** is defined as a collection of barcodes.  
- Each **barcode** consists of **two SNVs (single nucleotide variants)**.  

```text
Fingerprint
│
├── Barcode 1 ──> SNV1 + SNV2
├── Barcode 2 ──> SNV3 + SNV4
├── Barcode 3 ──> SNV5 + SNV6
└── Barcode N ──> SNV(2N-1) + SNV(2N)
```
----

## Installation

### Step 1 — Install RiSpy
You can install RiSpy directly from the repository:

```bash
pip install git+https://github.com/BIOIT-Sciensano/RiSpy.git
````

or, if you have cloned the repo locally:

```bash
pip install .
````

### Step 2 — Install additional dependencies with Conda
We provide an environment.yaml file with the required tools

```bash
conda env create -f environment.yaml
conda activate rispy_env
````
---
## Dependencies
This workflow was developed and tested using **Python 3.12**.   
In addition to the Python dependencies specified in `pyproject.toml`, **RiSpy** requires several bioinformatics tools.  

These are the versions we have **tested and validated** the workflow on:  

  - fastqc = 0.12.1
  - fastp = 1.0.1
  - bwa = 0.7.17
  - samtools = 1.21
  - qualimap = 2.2.2
  - bcftools = 1.21
  - ucsc-liftover = 469
  - bedtools = 2.31.0
  - minimap2 = 2.28
  - nanoplot = 1.48.0
  - seqkit = 2.8.2
  - chopper = 0.8.0
  - cd-hit = 4.6.8
  - deepvariant = 1.10.0

Other versions of these tools **may also work**, but they have not been tested by us.

### Manual installation
If you install these tools manually rather than through Conda,  
please ensure that the corresponding **binaries are available in your `PATH`**,  
otherwise RiSpy will not be able to run them.

**DeepVariant Integration:** RiSpy automatically installs *DeepVariant 1.10.0* the first time you run it.
> The required binaries and Python environment are set up and added to your PATH so  
  tools like make_examples, call_variants, and postprocess_variants are available system‑wide without any manual configuration.  
> If you prefer to use a different DeepVariant version, you can install it yourself and point your environment to that installation. To do so:  
> Generate the DeepVariant binaries and include their installation directory in your PATH environment variable.  
> Please follow the official installation instructions here: [DeepVariant Quick Start](https://github.com/google/deepvariant/blob/r1.9/docs/deepvariant-build-test.md)

---

### Environment & Dependency Validation
To ensure that your environment is correctly configured and that all required dependencies are available, we provide minimal Illumina and ONT test datasets.

### illugen Validation Test
Run the following command to validate the illugen workflow:
```bash
rispy illugen \
  --fastq-R1 test_data/Illumina_fastq/test_R1.fastq.gz \
  --fastq-R2 test_data/Illumina_fastq/test_R2.fastq.gz \
  --pickle-file test_data/test.pkl \
  --output-dir output_dir
```

### ontgen validation test
Run the following command to validate the ontgen workflow:
```bash
rispy ontgen \
  --fastq test_data/ONT_fastq/test_ONT.fastq.gz \
  --pickle-file test_data/test.pkl \
  --output-dir output_dir
```

### illudetect validation test
Run the following command to validate the illudetect workflow:
```bash
rispy illudetect \
  --fastq-R1 test_data/Illumina_fastq/test_R1.fastq.gz \
  --fastq-R2 test_data/Illumina_fastq/test_R2.fastq.gz \
  --barcodes test_data/test_barcodes.tsv \
  --ref-vcf test_data/test.vcf \
  --output-dir output_dir
```

### ontdetect validation test
Run the following command to validate the ontdetect workflow:
```bash
rispy ontdetect \
  --fastq test_data/ONT_fastq/test_ONT.fastq.gz \
  --barcodes test_data/test_barcodes.tsv \
  --ref-vcf test_data/test.vcf \
  --output-dir output_dir
```

## Command-Line Parameters

`RiSpy` workflow provides multiple subcommands.  
Each subcommand has its own parameter set.

---

## illugen

### Description:

```markdown
Analysis of Illumina sequencing data, including read QC, trimming, mapping, variant calling, and fingerprint generation.
NOTE: Before fingerprint generation is executed, RiSpy always performs fingerprint detection first 
      to verify whether there is already a fingerprint in the database that matches the provided dataset; 
      only after this check does the fingerprint‑generation workflow run.
      At the end of the analysis, the tool generates one report with the fingerprint detection and fingerprint generation results. 
      The fingerprint detection results allow users to assess the uniqueness of their line and identify whether 
      similar or matching lines already exist in the reference database. If a matching variety is detected, 
      the generated fingerprint may not be able to distinguish the sample from the matching line(s) present in the database.
```
### Expected outputs:

```markdown
OUTPUT_DIR/report.html - HTML report combining the QC, the fingerprint detection against the RiSpy database and the fingerprint generation results

OUTPUT_DIR/sample/fingerprint/2snp_fingerprints_annotated_withranks.tsv – generated 2-SNV barcodes reported based on the 3KRG genomics coordinates

OUTPUT_DIR/sample/fingerprint/barcodes.vcf - generated ref-vcf for the fingerprint

OUTPUT_DIR/sample/fastqc_before_trimming/ – fastqc results for raw reads

OUTPUT_DIR/sample/fastqc_after_trimming/ – fastqc results for trimmed reads

OUTPUT_DIR/sample/depth/ – mean and median depth of coverage results

OUTPUT_DIR/sample/breadth/ – breadth of coverage results

OUTPUT_DIR/sample/count/ – number of primary reads mapped

OUTPUT_DIR/sample/fold_80/ – fold80 result

OUTPUT_DIR/sample/mapping/ – bam file

OUTPUT_DIR/sample/qualimap/ – bam qc results

OUTPUT_DIR/sample/trimming/ – trimmed reads

OUTPUT_DIR/sample/deep_variant/ – raw and filtered variants results with stats
```

## usage: 

```bash
rispy illugen [OPTIONS]
```
## Options:

```markdown
 *  --fastq-R1                                -f1       PATH     Paired-end Illumina FASTQ file: forward reads (R1) [required]                                                                                                                                                                                    
 *  --fastq-R2                                -f2       PATH     Paired-end Illumina FASTQ file: reverse reads (R2) [required]                                                                                                                                                                                    
 *  --output-dir                              -o        PATH     Output directory [required]                                                                                                                                                                                                                      
    --call-variants-checkpoint-model          -vcm      PATH     Call variants checkpoint model path [default: /src/rispy_cli/deepvariant/1.10.0/checkpoints/wgs/deepvariant.wgs.ckpt]                                 
    --max-barcodes                            -mb       INTEGER  Maximum number of 2-SNV barcodes to generate (default: 100) [default: 100]                                                                                                                                                                       
    --identification-threshold                -it       INTEGER  Minimum number of detected barcodes for a line to be considered identified [default: 43]
    --hq-snp-positions                        -hsp      PATH     Tab‑delimited file listing SNP index, chromosome, and genomic position for 160K high‑quality SNPs [default: /src/rispy_cli/data/3K-HDRA-snp-comm-miss5pc_index_position.tsv]                                                                            
    --pickle-file                             -pf       PATH     Pickle file containing the high‑quality genotype matrix with a 'VARIETY' column and SNP genotype columns [default: /src/rispy_cli/data/hq_genotype.pkl]                                                                                                        
    --snp-info                                -si       PATH     Tab‑delimited SNP table with columns: SNP_INDEX, CHROMOSOME, POSITION, REFCALL for 32 million SNPs [default: /src/rispy_cli/data/3kall_snpposition_map.tsv]                                                                                              
    --cluster-json                            -cj       PATH     JSON file mapping cluster names to lists of accessions IDs [default: /src/rispy_cli/data/clusters.json]                                                 
    --representative-accessions               -ra       PATH     FASTA file containing representative sequences for each cluster [default: /src/rispy_cli/data/representatives.fasta]                                  
    --threads                                 -t        INTEGER  Number of CPU cores [default: 8]                                                                                                                                                                                                                 
    --mem-mb                                  -m        INTEGER  Memory in MB [default: 16000]                                                                                                                                                                                                                    
    --dry-run                                                    Generate the config file and perform a dry-run without executing rules                                                                                                                                                                           
    --fastp-window-size                                 INTEGER  window size for fastp [default: 4]                                                                                                                                                                                                               
    --fastp-quality                                     INTEGER  quality cutoff for fastp [default: 20]                                                                                                                                                                                                           
    --fastp-min-len                                     INTEGER  minimum length for fastp [default: 40]                                                                                                                                                                                                           
    --variant-filtering-min-QUAL                        INTEGER  minimum QUAL for variant filtering [default: 20]                                                                                                                                                                                                 
    --variant-filtering-min-DP                          INTEGER  minimum DP for variant filtering [default: 5]    
    --regions                                 -r        PATH     BED file restricting variant calling (DeepVariant) to these regions [default: /src/rispy_cli/data/Generation_SNVs.bed]
    --breadth-of-coverage-depth-cutoff                  INTEGER  depth cutoff for calculating the breadth of coverage [default: 5]                                                                                                                                                                                
    --fastqc-raw-threads                                INTEGER  threads for fastqc_raw [default: 4]                                                                                                                                                                                                              
    --fastqc-raw-memory                                 INTEGER  memory (MB) for fastqc_raw [default: 8000]                                                                                                                                                                                                       
    --trim-reads-threads                                INTEGER  threads for trim_reads [default: 8]                                                                                                                                                                                                              
    --trim-reads-memory                                 INTEGER  memory (MB) for trim_reads [default: 16000]                                                                                                                                                                                                      
    --fastqc-trimmed-threads                            INTEGER  threads for fastqc_trimmed [default: 4]                                                                                                                                                                                                          
    --fastqc-trimmed-memory                             INTEGER  memory (MB) for fastqc_trimmed [default: 8000]                                                                                                                                                                                                   
    --consensus-fasta-index-threads                     INTEGER  threads for consensus_fasta_index [default: 2]                                                                                                                                                                                                   
    --consensus-fasta-index-memory                      INTEGER  memory (MB) for consensus_fasta_index [default: 4000]                                                                                                                                                                                            
    --alignment-threads                                 INTEGER  threads for alignment [default: 8]                                                                                                                                                                                                               
    --alignment-memory                                  INTEGER  memory (MB) for alignment [default: 16000]                                                                                                                                                                                                       
    --qualimap-threads                                  INTEGER  threads for qualimap [default: 8]                                                                                                                                                                                                                
    --qualimap-memory                                   INTEGER  memory (MB) for qualimap [default: 16000]                                                                                                                                                                                                        
    --samtools-depth-threads                            INTEGER  threads for samtools_depth [default: 4]                                                                                                                                                                                                          
    --samtools-depth-memory                             INTEGER  memory (MB) for samtools_depth [default: 2000]                                                                                                                                                                                                   
    --calculate-mean-and-median-depth-memory            INTEGER  memory (MB) for calculate_mean_and_median_depth [default: 16000]                                                                                                                                                                                 
    --breadth-of-coverage-memory                        INTEGER  memory (MB) for breadth_of_coverage [default: 2000]                                                                                                                                                                                              
    --count-bam-threads                                 INTEGER  threads for count_bam [default: 8]                                                                                                                                                                                                               
    --count-bam-memory                                  INTEGER  memory (MB) for count_bam [default: 2000]                                                                                                                                                                                                        
    --fold-80-memory                                    INTEGER  memory (MB) for fold_80 [default: 16000]                                                                                                                                                                                                         
    --make-examples-memory                              INTEGER  memory (MB) for make-examples [default: 16000]                                                                                                                                                                                                   
    --call-variants-memory                              INTEGER  memory (MB) for call-variants [default: 16000]                                                                                                                                                                                                   
    --postprocess-variants-threads                      INTEGER  threads for postprocess-variants [default: 8]                                                                                                                                                                                                    
    --postprocess-variants-memory                       INTEGER  memory (MB) for postprocess-variants [default: 16000]                                                                                                                                                                                            
    --variant-filtering-memory                          INTEGER  memory (MB) for variant_filtering [default: 2000]                                                                                                                                                                                                
    --hq-positions-memory                               INTEGER  memory (MB) for hq_positions [default: 16000]                                                                                                                                                                                                    
    --liftover-memory                                   INTEGER  memory (MB) for liftover [default: 2000]                                                                                                                                                                                                         
    --chromosome-name-conversion-memory                 INTEGER  memory (MB) for chromosome_name_conversion [default: 2000]                                                                                                                                                                                       
    --hq-variants-memory                                INTEGER  memory (MB) for hq_variants [default: 16000]                                                                                                                                                                                                     
    --genotype-detection-memory                         INTEGER  memory (MB) for genotype_detection [default: 16000]                                                                                                                                                                                              
    --pickle-to-fasta-memory                            INTEGER  memory (MB) for pickle_to_fasta [default: 16000]                                                                                                                                                                                                 
    --clustering-cd-hit-threads                         INTEGER  threads for clustering_cd_hit [default: 8]                                                                                                                                                                                                       
    --clustering-cd-hit-memory                          INTEGER  memory (MB) for clustering_cd_hit [default: 16000]                                                                                                                                                                                               
    --update-cluster-json-memory                        INTEGER  memory (MB) for update_cluster_json [default: 16000]                                                                                                                                                                                             
    --feature-selection-memory                          INTEGER  memory (MB) for feature_selection [default: 16000]                                                                                                                                                                                               
    --add-info-to-indices-memory                        INTEGER  memory (MB) for add_info_to_indices [default: 2000]                                                                                                                                                                                              
    --list-remaining-accessions-memory                  INTEGER  memory (MB) for list_remaining_accessions [default: 16000]                                                                                                                                                                                       
    --screen-2-snp-combinations-memory                  INTEGER  memory (MB) for screen_2_snp_combinations [default: 16000]                                                                                                                                                                                       
    --add-info-to-snp-combinations-memory               INTEGER  memory (MB) for add_info_to_snp_combinations [default: 2000]                                                                                                                                                                                     
    --unique-snps-memory                                INTEGER  memory (MB) for unique_SNPs [default: 2000]                                                                                                                                                                                                      
    --bed-maker-memory                                  INTEGER  memory (MB) for bed_maker [default: 2000]                                                                                                                                                                                                        
    --vcf-maker-memory                                  INTEGER  memory (MB) for vcf_maker [default: 2000]                                                                                                                                                                                                        
    --report-html-memory                                INTEGER  memory (MB) for report_html [default: 2000]                                                                                                                                                                                                      
    --help                                                       Show this message and exit.                                                   
```

---

## ontgen

### Description:

```markdown
Analysis of ONT sequencing data, including read QC, trimming, mapping, variant calling, and fingerprint generation.
NOTE: Before fingerprint generation is executed, RiSpy always performs fingerprint detection first 
      to verify whether there is already a fingerprint in the database that matches the provided dataset; 
      only after this check does the fingerprint‑generation workflow run.
      At the end of the analysis, the tool generates one report with the fingerprint detection and fingerprint generation results. 
      The fingerprint detection results allow users to assess the uniqueness of their line and identify whether 
      similar or matching lines already exist in the reference database. If a matching variety is detected, 
      the generated fingerprint may not be able to distinguish the sample from the matching line(s) present in the database. 
```
### Expected outputs:

```markdown
OUTPUT_DIR/report.html - HTML report combining the QC, the fingerprint detection against the RiSpy database and the fingerprint generation results

OUTPUT_DIR/sample/fingerprint/2snp_fingerprints_annotated_withranks.tsv – generated 2-SNV barcodes reported based on the 3KRG genomics coordinates

OUTPUT_DIR/sample/fingerprint/barcodes.vcf - generated ref-vcf for the fingerprint

OUTPUT_DIR/sample/nanoplot_raw/ – nanoplot results for raw reads

OUTPUT_DIR/sample/nanoplot_after_filtering/ – nanoplot results for trimmed reads

OUTPUT_DIR/sample/depth/ – mean and median depth of coverage results

OUTPUT_DIR/sample/breadth/ – breadth of coverage results

OUTPUT_DIR/sample/count/ – number of primary reads mapped

OUTPUT_DIR/sample/fold_80/ – fold80 result

OUTPUT_DIR/sample/mapping/ – bam file

OUTPUT_DIR/sample/qualimap/ – bam qc results

OUTPUT_DIR/sample/chopper/ – trimmed reads

OUTPUT_DIR/sample/deep_variant/ – raw and filtered variants results with stats
```

## usage: 

```bash
rispy ontgen [OPTIONS]
```
## Options:

```markdown
 *  --fastq                                   -f        PATH     ONT FASTQ file [required]                                                                                                                                                                                                                        
 *  --output-dir                              -o        PATH     Output directory [required]                                                                                                                                                                                                                      
    --call-variants-checkpoint-model          -vcm      PATH     Call variants checkpoint model path [default: /src/rispy_cli/deepvariant/1.10.0/checkpoints/ont/deepvariant.ont.ckpt]                                 
    --max-barcodes                            -mb       INTEGER  Maximum number of 2-SNV barcodes to generate (default: 100) [default: 100]                                                                                                                                                                       
    --identification-threshold                -it       INTEGER  Minimum number of detected barcodes for a line to be considered identified [default: 43]
    --hq-snp-positions                        -hsp      PATH     Tab‑delimited file listing SNP index, chromosome, and genomic position for 160K high‑quality SNPs [default: /src/rispy_cli/data/3K-HDRA-snp-comm-miss5pc_index_position.tsv]                                                                            
    --pickle-file                             -pf       PATH     Pickle file containing the high‑quality genotype matrix with a 'VARIETY' column and SNP genotype columns [default: /src/rispy_cli/data/hq_genotype.pkl]                                                                                                        
    --snp-info                                -si       PATH     Tab‑delimited SNP table with columns: SNP_INDEX, CHROMOSOME, POSITION, REFCALL for 32 million SNPs [default: /src/rispy_cli/data/3kall_snpposition_map.tsv]                                                                                              
    --cluster-json                            -cj       PATH     JSON file mapping cluster names to lists of accessions IDs [default: /src/rispy_cli/data/clusters.json]                                                 
    --representative-accessions               -ra       PATH     FASTA file containing representative sequences for each cluster [default: /src/rispy_cli/data/representatives.fasta]                                  
    --threads                                 -t        INTEGER  Number of CPU cores [default: 8]                                                                                                                                                                                                                 
    --mem-mb                                  -m        INTEGER  Memory in MB [default: 16000]                                                                                                                                                                                                                    
    --dry-run                                                    Generate the config file and perform a dry-run without executing rules                                                                                                                                                                           
    --chopper-quality                                   INTEGER  quality cutoff for chopper [default: 10]                                                                                                                                                                                                         
    --chopper-min-len                                   INTEGER  minimum length for chopper [default: 150]                                                                                                                                                                                                        
    --breadth-of-coverage-depth-cutoff                  INTEGER  depth cutoff for calculating the breadth of coverage [default: 5]                                                                                                                                                                                
    --variant-filtering-min-QUAL                        INTEGER  minimum QUAL for variant filtering [default: 20]  
    --variant-filtering-min-DP                          INTEGER  minimum DP for variant filtering [default: 5]                                                                                                                                                                                                    
    --regions                                 -r        PATH     BED file restricting variant calling (DeepVariant) to these regions [default: /src/rispy_cli/data/Generation_SNVs.bed]
    --nanoplot-raw-threads                              INTEGER  threads for nanoplot_raw [default: 4]                                                                                                                                                                                                            
    --nanoplot-raw-memory                               INTEGER  memory (MB) for nanoplot_raw [default: 2000]                                                                                                                                                                                                     
    --extract-headers-threads                           INTEGER  threads for extract_headers [default: 2]                                                                                                                                                                                                         
    --extract-headers-memory                            INTEGER  memory (MB) for extract_headers [default: 2000]                                                                                                                                                                                                  
    --filter-reads-threads                              INTEGER  threads for filter_reads [default: 4]                                                                                                                                                                                                            
    --filter-reads-memory                               INTEGER  memory (MB) for filter_reads [default: 2000]                                                                                                                                                                                                     
    --trim-reads-threads                                INTEGER  threads for trim_reads [default: 4]                                                                                                                                                                                                              
    --trim-reads-memory                                 INTEGER  memory (MB) for trim_reads [default: 2000]                                                                                                                                                                                                       
    --nanoplot-after-filtering-threads                  INTEGER  threads for nanoplot_after_filtering [default: 4]                                                                                                                                                                                                
    --nanoplot-after-filtering-memory                   INTEGER  memory (MB) for nanoplot_after_filtering [default: 2000]                                                                                                                                                                                         
    --index-reference-threads                           INTEGER  threads for index_reference [default: 4]                                                                                                                                                                                                         
    --index-reference-memory                            INTEGER  memory (MB) for index_reference [default: 5000]                                                                                                                                                                                                  
    --alignment-threads                                 INTEGER  threads for alignment [default: 8]                                                                                                                                                                                                               
    --alignment-memory                                  INTEGER  memory (MB) for alignment [default: 16000]                                                                                                                                                                                                       
    --qualimap-threads                                  INTEGER  threads for qualimap [default: 8]                                                                                                                                                                                                                
    --qualimap-memory                                   INTEGER  memory (MB) for qualimap [default: 16000]                                                                                                                                                                                                        
    --samtools-depth-threads                            INTEGER  threads for samtools_depth [default: 2]                                                                                                                                                                                                          
    --samtools-depth-memory                             INTEGER  memory (MB) for samtools_depth [default: 2000]                                                                                                                                                                                                   
    --calculate-mean-and-median-depth-memory            INTEGER  memory (MB) for calculate_mean_and_median_depth [default: 16000]                                                                                                                                                                                 
    --breadth-of-coverage-memory                        INTEGER  memory (MB) for breadth_of_coverage [default: 2000]                                                                                                                                                                                              
    --count-bam-threads                                 INTEGER  threads for count_bam [default: 8]                                                                                                                                                                                                               
    --count-bam-memory                                  INTEGER  memory (MB) for count_bam [default: 2000]                                                                                                                                                                                                        
    --fold-80-memory                                    INTEGER  memory (MB) for fold_80 [default: 16000]                                                                                                                                                                                                         
    --make-examples-memory                              INTEGER  memory (MB) for make-examples [default: 16000]                                                                                                                                                                                                   
    --call-variants-memory                              INTEGER  memory (MB) for call-variants [default: 16000]                                                                                                                                                                                                   
    --postprocess-variants-threads                      INTEGER  threads for postprocess-variants [default: 8]                                                                                                                                                                                                    
    --postprocess-variants-memory                       INTEGER  memory (MB) for postprocess-variants [default: 16000]                                                                                                                                                                                            
    --variant-filtering-memory                          INTEGER  memory (MB) for variant_filtering [default: 2000]                                                                                                                                                                                                
    --hq-positions-memory                               INTEGER  memory (MB) for hq_positions [default: 16000]                                                                                                                                                                                                    
    --liftover-memory                                   INTEGER  memory (MB) for liftover [default: 2000]                                                                                                                                                                                                         
    --chromosome-name-conversion-memory                 INTEGER  memory (MB) for chromosome_name_conversion [default: 2000]                                                                                                                                                                                       
    --hq-variants-memory                                INTEGER  memory (MB) for hq_variants [default: 16000]                                                                                                                                                                                                     
    --genotype-detection-memory                         INTEGER  memory (MB) for genotype_detection [default: 16000]                                                                                                                                                                                              
    --pickle-to-fasta-memory                            INTEGER  memory (MB) for pickle_to_fasta [default: 16000]                                                                                                                                                                                                 
    --clustering-cd-hit-threads                         INTEGER  threads for clustering_cd_hit [default: 8]                                                                                                                                                                                                       
    --clustering-cd-hit-memory                          INTEGER  memory (MB) for clustering_cd_hit [default: 16000]                                                                                                                                                                                               
    --update-cluster-json-memory                        INTEGER  memory (MB) for update_cluster_json [default: 16000]                                                                                                                                                                                             
    --feature-selection-memory                          INTEGER  memory (MB) for feature_selection [default: 16000]                                                                                                                                                                                               
    --add-info-to-indices-memory                        INTEGER  memory (MB) for add_info_to_indices [default: 2000]                                                                                                                                                                                              
    --list-remaining-accessions-memory                  INTEGER  memory (MB) for list_remaining_accessions [default: 16000]                                                                                                                                                                                       
    --screen-2-snp-combinations-memory                  INTEGER  memory (MB) for screen_2_snp_combinations [default: 16000]                                                                                                                                                                                       
    --add-info-to-snp-combinations-memory               INTEGER  memory (MB) for add_info_to_snp_combinations [default: 2000]                                                                                                                                                                                     
    --unique-snps-memory                                INTEGER  memory (MB) for unique_SNPs [default: 2000]                                                                                                                                                                                                      
    --bed-maker-memory                                  INTEGER  memory (MB) for bed_maker [default: 2000]                                                                                                                                                                                                        
    --vcf-maker-memory                                  INTEGER  memory (MB) for vcf_maker [default: 2000]                                                                                                                                                                                                        
    --report-html-memory                                INTEGER  memory (MB) for report_html [default: 2000]                                                                                                                                                                                                      
    --help                                                       Show this message and exit.                            
```

---

## illudetect

### Description:

```markdown
Analysis of Illumina sequencing data, including read QC, trimming, mapping, variant calling, and fingerprint detection.
```
### Expected outputs:

```markdown
OUTPUT_DIR/report.html - HTML report with summary of results

OUTPUT_DIR/sample/Fingerprint_Detection – fingerprint detection results

OUTPUT_DIR/sample/fastqc_before_trimming/ – fastqc results for raw reads

OUTPUT_DIR/sample/fastqc_after_trimming/ – fastqc results for trimmed reads

OUTPUT_DIR/sample/depth/ – mean and median depth of coverage results

OUTPUT_DIR/sample/breadth/ – breadth of coverage results

OUTPUT_DIR/sample/count/ – number of primary reads mapped

OUTPUT_DIR/sample/fold_80/ – fold80 result

OUTPUT_DIR/sample/mapping/ – bam file

OUTPUT_DIR/sample/qualimap/ – bam qc results

OUTPUT_DIR/sample/trimming/ – trimmed reads

OUTPUT_DIR/sample/deep_variant/ – raw and filtered variants results with stats
```

## usage: 

```bash
rispy illudetect [OPTIONS]
```
## Options:

```markdown
 *  --fastq-R1                                -f1       PATH     Paired-end Illumina FASTQ file: forward reads (R1) [required]                                                                                                                                                                                    
 *  --fastq-R2                                -f2       PATH     Paired-end Illumina FASTQ file: reverse reads (R2) [required]
 *  --output-dir                              -o        PATH     Output directory [required]
    --barcodes                                -b        PATH     tsv file of the fingerprint (2-SNV barcodes). If omitted, all 3KRG fingerprints are tested.                                                                                                                                                      
    --ref-vcf                                 -v        PATH     VCF used to generate the fingerprint. If omitted, all 3KRG fingerprints are tested.                                                                                                                                                              
    --top-n                                   -n        INTEGER  Number of top-matching rice lines to show in the report. Must be used if no fingerprint is passed. [default: 10]
    --identification-threshold                -it       INTEGER  Minimum number of detected barcodes for a line to be considered identified [default: 43]
    --call-variants-checkpoint-model          -vcm      PATH     Call variants checkpoint model path [default: /src/rispy_cli/deepvariant/1.10.0/checkpoints/wgs/deepvariant.wgs.ckpt]                                 
    --dry-run                                                    Generate the config file and perform a dry-run without executing rules                                                                                                                                                                           
    --centromere                              -c        PATH     Centromere TSV file [default: /src/rispy_cli/data/Centromere.tsv]                                                                                     
    --threads                                 -t        INTEGER  Number of CPU cores [default: 8]                                                                                                                                                                                                                 
    --mem-mb                                  -m        INTEGER  Memory in MB [default: 16000]                                                                                                                                                                                                                    
    --fastp-window-size                                 INTEGER  window size for fastp [default: 4]                                                                                                                                                                                                               
    --fastp-quality                                     INTEGER  quality cutoff for fastp [default: 20]                                                                                                                                                                                                           
    --fastp-min-len                                     INTEGER  minimum length for fastp [default: 40]                                                                                                                                                                                                           
    --breadth-of-coverage-depth-cutoff                  INTEGER  depth cutoff for calculating the breadth of coverage [default: 5]
    --variant-filtering-min-QUAL                        INTEGER  minimum QUAL for variant filtering [default: 20]                                                                                                                                                                                                 
    --variant-filtering-min-DP                          INTEGER  minimum DP for variant filtering [default: 5]                                                                                                                                                                                                    
    --regions                                 -r        PATH     BED file restricting variant calling (DeepVariant) to these regions [default: /src/rispy_cli/data/Fingerprint_SNVs.bed]
    --fastqc-raw-threads                                INTEGER  threads for fastqc_raw [default: 4]                                                                                                                                                                                                              
    --fastqc-raw-memory                                 INTEGER  memory (MB) for fastqc_raw [default: 8000]                                                                                                                                                                                                       
    --trim-reads-threads                                INTEGER  threads for trim_reads [default: 8]                                                                                                                                                                                                              
    --trim-reads-memory                                 INTEGER  memory (MB) for trim_reads [default: 16000]                                                                                                                                                                                                      
    --fastqc-trimmed-threads                            INTEGER  threads for fastqc_trimmed [default: 4]                                                                                                                                                                                                          
    --fastqc-trimmed-memory                             INTEGER  memory (MB) for fastqc_trimmed [default: 8000]                                                                                                                                                                                                   
    --consensus-fasta-index-threads                     INTEGER  threads for consensus_fasta_index [default: 2]                                                                                                                                                                                                   
    --consensus-fasta-index-memory                      INTEGER  memory (MB) for consensus_fasta_index [default: 4000]                                                                                                                                                                                            
    --alignment-threads                                 INTEGER  threads for alignment [default: 8]                                                                                                                                                                                                               
    --alignment-memory                                  INTEGER  memory (MB) for alignment [default: 16000]                                                                                                                                                                                                       
    --qualimap-threads                                  INTEGER  threads for qualimap [default: 8]                                                                                                                                                                                                                
    --qualimap-memory                                   INTEGER  memory (MB) for qualimap [default: 16000]                                                                                                                                                                                                        
    --samtools-depth-threads                            INTEGER  threads for samtools_depth [default: 4]                                                                                                                                                                                                          
    --samtools-depth-memory                             INTEGER  memory (MB) for samtools_depth [default: 2000]                                                                                                                                                                                                   
    --calculate-mean-and-median-depth-memory            INTEGER  memory (MB) for calculate_mean_and_median_depth [default: 16000]                                                                                                                                                                                 
    --breadth-of-coverage-memory                        INTEGER  memory (MB) for breadth_of_coverage [default: 2000]                                                                                                                                                                                              
    --count-bam-threads                                 INTEGER  threads for count_bam [default: 8]                                                                                                                                                                                                               
    --count-bam-memory                                  INTEGER  memory (MB) for count_bam [default: 2000]                                                                                                                                                                                                        
    --fold-80-memory                                    INTEGER  memory (MB) for fold_80 [default: 16000]                                                                                                                                                                                                         
    --make-examples-memory                              INTEGER  memory (MB) for make-examples [default: 16000]                                                                                                                                                                                                   
    --call-variants-memory                              INTEGER  memory (MB) for call-variants [default: 16000]                                                                                                                                                                                                   
    --postprocess-variants-threads                      INTEGER  threads for postprocess-variants [default: 8]                                                                                                                                                                                                    
    --postprocess-variants-memory                       INTEGER  memory (MB) for postprocess-variants [default: 16000]                                                                                                                                                                                            
    --variant-filtering-memory                          INTEGER  memory (MB) for variant_filtering [default: 2000]                                                                                                                                                                                                
    --unique-SNPs-memory                                INTEGER  memory (MB) for unique_SNPs [default: 4000]                                                                                                                                                                                                      
    --bed-maker-memory                                  INTEGER  memory (MB) for bed_maker [default: 2000]                                                                                                                                                                                                        
    --liftover-memory                                   INTEGER  memory (MB) for liftover [default: 4000]                                                                                                                                                                                                         
    --chromosome-name-conversion-memory                 INTEGER  memory (MB) for chromosome_name_conversion [default: 4000]                                                                                                                                                                                       
    --filter-fingerprints-memory                        INTEGER  memory (MB) for filter_fingerprints [default: 4000]                                                                                                                                                                                              
    --validation-fingerprints-alt-memory                INTEGER  memory (MB) for validation_fingerprints_alt [default: 4000]                                                                                                                                                                                      
    --validation-fingerprints-ref-memory                INTEGER  memory (MB) for validation_fingerprints_ref [default: 4000]                                                                                                                                                                                      
    --barcodes-positions-new-ref-memory                 INTEGER  memory (MB) for barcodes_positions_new_ref [default: 4000]                                                                                                                                                                                       
    --chromosome-arms-memory                            INTEGER  memory (MB) for chromosome_arms [default: 4000]                                                                                                                                                                                                  
    --aggregate-fingerprint-matches-memory              INTEGER  memory (MB) aggregate_fingerprint_matches [default: 4000]                                                                                                                                                                                        
    --report-html-memory                                INTEGER  memory (MB) for report_html [default: 2000]                                                                                                                                                                                                      
    --help                                                       Show this message and exit.               
```

---

## ontdetect

### Description:

```markdown
Analysis of ONT sequencing data, including read QC, trimming, mapping, variant calling, and fingerprint detection.
```
### Expected outputs:

```markdown
OUTPUT_DIR/report.html - HTML report with summary of results

OUTPUT_DIR/sample/Fingerprint_Detection – fingerprint detection results

OUTPUT_DIR/sample/nanoplot_raw/ – nanoplot results for raw reads

OUTPUT_DIR/sample/nanoplot_after_filtering/ – nanoplot results for trimmed reads

OUTPUT_DIR/sample/depth/ – mean and median depth of coverage results

OUTPUT_DIR/sample/breadth/ – breadth of coverage results

OUTPUT_DIR/sample/count/ – number of primary reads mapped

OUTPUT_DIR/sample/fold_80/ – fold80 result

OUTPUT_DIR/sample/mapping/ – bam file

OUTPUT_DIR/sample/qualimap/ – bam qc results

OUTPUT_DIR/sample/chopper/ – trimmed reads

OUTPUT_DIR/sample/deep_variant/ – raw and filtered variants results with stats
```

## usage:

```bash
rispy ontdetect [OPTIONS]
```
## Options:

```markdown
 *  --fastq                                   -f        PATH     ONT FASTQ file [required]                                                                                                                                                                                                                        
 *  --output-dir                              -o        PATH     Output directory [required]                                                                                                                                                                                                                      
    --barcodes                                -b        PATH     tsv file of the fingerprint (2-SNV barcodes). If omitted, all 3KRG fingerprints are tested.                                                                                                                                                      
    --ref-vcf                                 -v        PATH     VCF used to generate the fingerprint. If omitted, all 3KRG fingerprints are tested.                                                                                                                                                              
    --top-n                                   -n        INTEGER  Number of top-matching rice lines to show in the report. Must be used if no fingerprint is passed. [default: 10]                                                                                                                                    
    --identification-threshold                -it       INTEGER  Minimum number of detected barcodes for a line to be considered identified [default: 43]
    --call-variants-checkpoint-model          -vcm      PATH     Call variants checkpoint model path [default: /src/rispy_cli/deepvariant/1.10.0/checkpoints/ont/deepvariant.ont.ckpt]                                 
    --dry-run                                                    Generate the config file and perform a dry-run without executing rules                                                                                                                                                                           
    --centromere                              -c        PATH     Centromere TSV file [default: /src/rispy_cli/data/Centromere.tsv]                                                                                     
    --threads                                 -t        INTEGER  Number of CPU cores [default: 8]                                                                                                                                                                                                                 
    --mem-mb                                  -m        INTEGER  Memory in MB [default: 16000]                                                                                                                                                                                                                    
    --chopper-quality                                   INTEGER  quality cutoff for chopper [default: 10]                                                                                                                                                                                                         
    --chopper-min-len                                   INTEGER  minimum length for chopper [default: 150]  
    --breadth-of-coverage-depth-cutoff                  INTEGER  depth cutoff for calculating the breadth of coverage [default: 5]                                                                                                                                                                                
    --variant-filtering-min-QUAL                        INTEGER  minimum QUAL for variant filtering [default: 20]                                                                                                                                                                                                 
    --variant-filtering-min-DP                          INTEGER  minimum DP for variant filtering [default: 5]                                                                                                                                                                                                    
    --regions                                 -r        PATH     BED file restricting variant calling (DeepVariant) to these regions [default: /src/rispy_cli/data/Fingerprint_SNVs.bed]
    --nanoplot-raw-threads                              INTEGER  threads for nanoplot_raw [default: 4]                                                                                                                                                                                                            
    --nanoplot-raw-memory                               INTEGER  memory (MB) for nanoplot_raw [default: 2000]                                                                                                                                                                                                     
    --extract-headers-threads                           INTEGER  threads for extract_headers [default: 2]                                                                                                                                                                                                         
    --extract-headers-memory                            INTEGER  memory (MB) for extract_headers [default: 2000]                                                                                                                                                                                                  
    --filter-reads-threads                              INTEGER  threads for filter_reads [default: 4]                                                                                                                                                                                                            
    --filter-reads-memory                               INTEGER  memory (MB) for filter_reads [default: 2000]                                                                                                                                                                                                     
    --trim-reads-threads                                INTEGER  threads for trim_reads [default: 4]                                                                                                                                                                                                              
    --trim-reads-memory                                 INTEGER  memory (MB) for trim_reads [default: 2000]                                                                                                                                                                                                       
    --nanoplot-after-filtering-threads                  INTEGER  threads for nanoplot_after_filtering [default: 4]                                                                                                                                                                                                
    --nanoplot-after-filtering-memory                   INTEGER  memory (MB) for nanoplot_after_filtering [default: 2000]                                                                                                                                                                                         
    --index-reference-threads                           INTEGER  threads for index_reference [default: 4]                                                                                                                                                                                                         
    --index-reference-memory                            INTEGER  memory (MB) for index_reference [default: 5000]                                                                                                                                                                                                  
    --alignment-threads                                 INTEGER  threads for alignment [default: 8]                                                                                                                                                                                                               
    --alignment-memory                                  INTEGER  memory (MB) for alignment [default: 16000]                                                                                                                                                                                                       
    --qualimap-threads                                  INTEGER  threads for qualimap [default: 8]                                                                                                                                                                                                                
    --qualimap-memory                                   INTEGER  memory (MB) for qualimap [default: 16000]                                                                                                                                                                                                        
    --samtools-depth-threads                            INTEGER  threads for samtools_depth [default: 2]                                                                                                                                                                                                          
    --samtools-depth-memory                             INTEGER  memory (MB) for samtools_depth [default: 2000]                                                                                                                                                                                                   
    --calculate-mean-and-median-depth-memory            INTEGER  memory (MB) for calculate_mean_and_median_depth [default: 16000]                                                                                                                                                                                 
    --breadth-of-coverage-memory                        INTEGER  memory (MB) for breadth_of_coverage [default: 2000]                                                                                                                                                                                              
    --count-bam-threads                                 INTEGER  threads for count_bam [default: 8]                                                                                                                                                                                                               
    --count-bam-memory                                  INTEGER  memory (MB) for count_bam [default: 2000]                                                                                                                                                                                                        
    --fold-80-memory                                    INTEGER  memory (MB) for fold_80 [default: 16000]                                                                                                                                                                                                         
    --make-examples-memory                              INTEGER  memory (MB) for make-examples [default: 16000]                                                                                                                                                                                                   
    --call-variants-memory                              INTEGER  memory (MB) for call-variants [default: 16000]                                                                                                                                                                                                   
    --postprocess-variants-threads                      INTEGER  threads for postprocess-variants [default: 8]                                                                                                                                                                                                    
    --postprocess-variants-memory                       INTEGER  memory (MB) for postprocess-variants [default: 16000]                                                                                                                                                                                            
    --variant-filtering-memory                          INTEGER  memory (MB) for variant_filtering [default: 2000]                                                                                                                                                                                                
    --unique-SNPs-memory                                INTEGER  memory (MB) for unique_SNPs [default: 4000]                                                                                                                                                                                                      
    --bed-maker-memory                                  INTEGER  memory (MB) for bed_maker [default: 2000]                                                                                                                                                                                                        
    --liftover-memory                                   INTEGER  memory (MB) for liftover [default: 8000]                                                                                                                                                                                                         
    --chromosome-name-conversion-memory                 INTEGER  memory (MB) for chromosome_name_conversion [default: 4000]                                                                                                                                                                                       
    --filter-fingerprints-memory                        INTEGER  memory (MB) for filter_fingerprints [default: 8000]                                                                                                                                                                                              
    --validation-fingerprints-alt-memory                INTEGER  memory (MB) for validation_fingerprints_alt [default: 4000]                                                                                                                                                                                      
    --validation-fingerprints-ref-memory                INTEGER  memory (MB) for validation_fingerprints_ref [default: 4000]                                                                                                                                                                                      
    --barcodes-positions-new-ref-memory                 INTEGER  memory (MB) for barcodes_positions_new_ref [default: 4000]                                                                                                                                                                                       
    --chromosome-arms-memory                            INTEGER  memory (MB) for chromosome_arms [default: 4000]                                                                                                                                                                                                  
    --aggregate-fingerprint-matches-memory              INTEGER  memory (MB) aggregate_fingerprint_matches [default: 4000]                                                                                                                                                                                        
    --report-html-memory                                INTEGER  memory (MB) for report_html [default: 2000]                                                                                                                                                                                                      
    --help                                                       Show this message and exit.                     
```

---

> **Note:** Variant calling is restricted by default: `illudetect`/`ontdetect` use the fingerprint SNV positions (`/src/rispy_cli/data/Fingerprint_SNVs.bed`),  
> `illugen`/`ontgen` use the 160K high-quality SNP positions (`/src/rispy_cli/data/Generation_SNVs.bed`).

> **Note:** The default number of threads is set conservatively and may not be optimal for large FASTQ datasets.  
> For faster analysis, we recommend increasing the thread count for compute-intensive steps   
> particularly **alignment** and **variant calling** up to the maximum supported by your system.

---

## Help
All commands support:  
```bash
rispy <command> --help
```
## Fingerprint Data
The fingerprint atlas for the 3KRG lines — including Kitaake, which is not present in the original 3KRG database — are stored in
```markdown
/src/rispy_cli/fingerprints/  
```
and are indexed by their accession IDs.  

To look up the correct accession ID, refer to:
```markdown
/src/rispy_cli/data/3kall_variety_map.tsv 
```
and check the VARIETY_INDEX column.

## Automatic Data Download
The first time you run 
```bash
rispy --help
```
or any other rispy command, the file
```markdown
data/3kall_variety_map.tsv 
``` 
along with all other required reference data will be automatically downloaded and placed in the appropriate directories.

## Citation
If you use this tool, please consider citing our [publication](https://doi.org/10.1093/bib/bbag406).  
Zolfaghari, A., Fraiture, M‑A., Vanneste, K., Stuyts, A., Frouin, J., Meunier, A‑C., Deforce, D., Roosens, N. H. C., & D’aes, J. (2026).  
RiSpy: A feature selection‑based fingerprinting framework for accurate identification of genome‑edited rice lines.  
Briefings in Bioinformatics, 27(4), bbag406.