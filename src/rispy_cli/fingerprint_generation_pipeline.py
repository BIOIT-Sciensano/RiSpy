from pathlib import Path
import yaml
import pandas as pd
import numpy as np
import pysam
import os
import pickle
import re
import json
import csv
from sklearn.preprocessing import OrdinalEncoder
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import SelectKBest, chi2

def _check_file_exists(*paths: str) -> None:
    """
    Raise FileNotFoundError if any required file does not exist.
    """
    for path in paths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required file not found: {path}")


def hq_positions(input: str, output: str) -> None:
    """
    Convert a high‑quality SNP table into a BED file containing
    chromosome, start, end, and SNP identifier fields.

    Parameters
    ----------
    input : str
        Path to a TSV file containing high‑quality SNP positions.

    output : str
        Path to the BED file to write.

    Output
    ------
    BED file with columns:
        chrom, start, end, id
    """
    _check_file_exists(input)
    df = pd.read_csv(input, sep='\t', skiprows=1)
    # Create BED format columns
    df['chrom'] = 'chr' + df.iloc[:, 1].astype(str)
    df['start'] = df.iloc[:, 2].astype(int) - 1  # POS (0-based)
    df['end'] = df.iloc[:, 2].astype(int)
    df['id'] = 'ID_' + df.iloc[:, 0].astype(str)
    bed_df = df[['chrom', 'start', 'end', 'id']]
    bed_df.to_csv(output, sep='\t', header=False, index=False)

def chromosome_name_conversion(chromosome_mapping_file: str, HQ_BED_LIFTOVER: str, output: str) -> None:
    """
    Replace chromosome names in a lifted‑over BED file using a user‑provided
    chromosome mapping table.

    Parameters
    ----------
    chromosome_mapping_file : str
        Path to a text file defining chromosome name mappings.

    HQ_BED_LIFTOVER : str
        Path to the BED file produced after liftover.

    output : str
        Path to the output BED file with updated chromosome names.

    Output
    ------
    BED file identical to the input except that:
        - The chromosome field is replaced according to the mapping
          provided in `chromosome_mapping_file`.
        - All other fields are preserved exactly as in the lifted‑over BED
    """
    _check_file_exists(chromosome_mapping_file)
    chromosome_mapping = {}
    with open(chromosome_mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chromosome_mapping[parts[0]] = parts[1]
    with open(HQ_BED_LIFTOVER, 'r') as f, open(output, 'w') as out:
        for line in f:
            parts = line.strip().split()
            if parts[0] in chromosome_mapping:
                parts[0] = chromosome_mapping[parts[0]]  # Replace chromosome name
            out.write('\t'.join(parts) + '\n')

def genotype_detection(HQ_VARIANTS: str, PICKLEFILE: str, UNMAPPED: str, Low_Confidence_Ref_Alleles: str, REF_POSITIONS: str, BAM: str, mapping_tsv: str, fasta_path: str, output: str) -> None:
    """
    Process high‑quality variants, low‑coverage reference positions, and
    liftover‑unmapped SNPs to generate a new genotype row and update an
    existing genotype matrix.

    Parameters
    ----------
    HQ_VARIANTS : str
        Path to a file containing high‑quality variant calls, including
        SNP IDs and ALT alleles.

    PICKLEFILE : str
        Path to a pickle file containing the existing genotype DataFrame.

    UNMAPPED : str
        Path to a file listing SNPs that failed liftover.

    Low_Confidence_Ref_Alleles : str
        Path to a file listing SNPs whose reference allele confidence
        is low.

    REF_POSITIONS : str
        Path to a file containing reference SNP positions to be checked
        for sequencing depth.

    BAM : str
        Path to the BAM file used to compute per‑position coverage.

    Output
    ------
    Updated pickle file containing:
        - A new genotype row constructed from ALT alleles of high‑quality
          variants.
        - Removal of SNP columns that either failed liftover or had
          insufficient reference‑allele coverage (<5x).

    """
    _check_file_exists(HQ_VARIANTS, PICKLEFILE, UNMAPPED, Low_Confidence_Ref_Alleles, REF_POSITIONS, BAM)
    id_to_alt_allele = {}
    with open(HQ_VARIANTS, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            id_field = row[3]
            alt_allele = row[8]
            # Extract the numeric part of the ID (after the underscore)
            id_numeric = id_field.split('_')[1]
            id_to_alt_allele[id_numeric] = alt_allele
    id_to_low_quality_ref_allele = []
    with open(Low_Confidence_Ref_Alleles, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            id_field = row[3]
            # Extract the numeric part of the ID (after the underscore)
            id_numeric = str(id_field.split('_')[1])
            id_to_low_quality_ref_allele.append(id_numeric)
    with pysam.AlignmentFile(BAM, "rb") as bam:
        with open(REF_POSITIONS) as f:
            for line in f:
                chrom, start, end, snp_id = line.strip().split("\t")
                pos = int(end)  # BED end --> 1‑based POS
                ref_len = bam.get_reference_length(chrom)
                s = pos - 1
                e = pos
                if s < 0 or e > ref_len:
                    continue
                # per-base coverage arrays
                a, c, g, t = bam.count_coverage(
                    chrom,
                    start=pos - 1,
                    end=pos,
                    quality_threshold=0
                )
                # depth at this exact position
                depth = a[0] + c[0] + g[0] + t[0]
                numeric_id = str(snp_id.split("_")[1])
                if depth < 5 and numeric_id not in id_to_low_quality_ref_allele:
                    id_to_low_quality_ref_allele.append(numeric_id)
    # retrieve object from pickle
    with open(PICKLEFILE, 'rb') as f:
        X = pickle.load(f)
    # Add the new variety genotype to the DataFrame
    mapping = pd.read_csv(mapping_tsv, sep="\t", header=None, names=["chrom", "start", "end", "id"])
    # Convert ID_369 --> 369
    mapping["col"] = mapping["id"].str.replace("ID_", "").astype(int)
    # Load FASTA
    ref = pysam.FastaFile(fasta_path)
    new_row = {}
    for col in X.columns:
        if col == "VARIETY":
            continue
        col_int = int(col)
        # Check if this SNP exists in the mapping
        match = mapping.loc[mapping["col"] == col_int]
        if match.empty:
            # Column not found --> assign N
            new_row[col] = "N"
            continue
        row = match.iloc[0]
        chrom = row["chrom"]
        pos = row["end"]
        # Extract base from FASTA
        base = ref.fetch(chrom, pos - 1, pos).upper()
        if base == "":
            base = "N"
        new_row[col] = base
    # Add VARIETY value
    max_id = pd.to_numeric(X["VARIETY"]).max()
    new_row["VARIETY"] = int(max_id) + 1
    # Convert to pandas Series in correct column order
    base_row = pd.Series(new_row)[X.columns]
    allele_series = pd.Series(id_to_alt_allele)
    new_row = allele_series.combine_first(base_row)
    # Append the row to the DataFrame
    X.loc[len(X)] = new_row
    # Removing the columns that the SNP index was not able to liftover to the new reference genome
    df = pd.read_csv(UNMAPPED, sep='\t', comment='#', header=None, names=['chrom', 'start', 'end', 'id'])
    ids_to_remove = df['id'].str.split('_').str[1].astype(str)
    # Convert reference alleles list to series
    id_to_low_quality_ref_allele_series = pd.Series(id_to_low_quality_ref_allele, dtype=str)
    # Removing low confidence reference alleles
    ids_to_remove_final = pd.concat([ids_to_remove, id_to_low_quality_ref_allele_series], ignore_index=True).astype(str)
    # Make sure column names are strings for safe comparison
    cols_to_drop = [col for col in X.columns if str(col) in set(ids_to_remove_final)]
    X = X.drop(columns=cols_to_drop)
    # Save the modified DataFrame back to a pickle file
    X.to_pickle(output)

def pickle_to_fasta(input: str, output: str) -> None:
    """
    Extract the genotype of the cultivar of interest in a pickle‑based genotype
    matrix, sanitize all non‑ATCG alleles by converting them to 'N', and
    write the resulting sequence as a single‑entry FASTA file.

    Parameters
    ----------
    input : str
        Path to a pickle file containing a DataFrame with a 'VARIETY'
        column followed by SNP genotype columns.

    output : str
        Path to the FASTA file to be written, containing only the
        cultivar of interest.

    Output
    ------
    FASTA file with:
        - A single record whose header is the VARIETY identifier of the
          cultivar of interest.
        - A sequence constructed by uppercasing all genotype values and
          replacing any non‑ATCG character with 'N'.
    """
    _check_file_exists(input)
    with open(input, "rb") as f:
        df = pickle.load(f)
    last = df.iloc[-1]
    cultivar = str(last["VARIETY"])
    # Convert to string, uppercase, and replace non-ATCG with N
    bases = (
        last.drop(labels=["VARIETY"])
        .astype(str)
        .str.upper()
        .apply(lambda x: x if re.fullmatch(r"[ATCG]", x) else "N")
    )
    seq = "".join(bases.values)
    with open(output, "w") as out:
        out.write(f">{cultivar}\n")
        out.write(seq + "\n")

def update_cluster_json(clusters_cd_hit: str, Cluster_json: str, output: str) -> None:
    """
    Identify the CD‑HIT cluster to which a target cultivar belongs, update
    an existing JSON cluster mapping accordingly, and write the modified
    cluster structure to a new JSON file.

    Parameters
    ----------
    clusters_cd_hit : str
        Path to a CD‑HIT `.clstr` file produced by cd-hit-est-2d.
        The file is scanned to locate the cluster that
        contains the target cultivar identifier.

    Cluster_json : str
        Path to a JSON file mapping cluster names
        to lists of cultivar identifiers as strings.

    output : str
        Path to the updated JSON file that will store the modified
        cluster mapping after inserting the target cultivar.

    Output
    ------
    JSON file in which:
        - If the target cultivar appears in an existing CD‑HIT cluster,
          it is appended to the corresponding entry in the JSON mapping.
        - If the target cultivar does not appear in any cluster, a new
          cluster entry is created with the next available cluster index.
    """
    _check_file_exists(clusters_cd_hit, Cluster_json)
    # Determine next cultivar ID dynamically from existing JSON
    with open(Cluster_json, "r") as f:
        clusters = json.load(f)
    all_ids = []
    for members in clusters.values():
        for cid in members:
            all_ids.append(int(cid))
    target_id = str(max(all_ids) + 1)
    cluster_of_target = None
    current_cluster = None
    with open(clusters_cd_hit, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = line.split()[1]
                continue
            m = re.search(r">(.*?)\.\.\.", line)
            if m:
                seq_id = m.group(1)
                if seq_id == target_id:
                    cluster_of_target = current_cluster
                    break
    if cluster_of_target is not None:
        # Target belongs to an existing cluster
        key = f"Cluster {cluster_of_target}"
        if target_id not in clusters[key]:
            clusters[key].append(target_id)
    else:
        # Target is a new cluster --> create a new cluster ID
        existing_numbers = [int(k.split()[1]) for k in clusters.keys() if k.startswith("Cluster")]
        new_cluster_num = max(existing_numbers) + 1
        new_key = f"Cluster {new_cluster_num}"
        clusters[new_key] = [target_id]
    with open(output, "w") as f:
        json.dump(clusters, f, indent=2)

def feature_selection(picklefile1: str, cluster_json: str, feature_scores: str, picklefile2: str) -> None:
    """
    Restrict the genotype matrix to the representative cultivars of all
    clusters, retain only valid SNP features for the variety of interest,
    resolve missing or ambiguous genotypes, and perform chi‑square–based
    feature selection to identify the top 2,000 discriminatory SNPs.

    Parameters
    ----------
    picklefile1 : str
        Path to a pickle file containing the full genotype DataFrame with
        a 'VARIETY' column and SNP genotype columns.

    cluster_json : str
        Path to a JSON file mapping cluster names to lists of cultivar
        identifiers. The first identifier of each list is treated as the
        representative of that cluster.

    picklefile2 : str
        Path to a pickle file where the filtered and preprocessed
        genotype matrix (restricted to allowed varieties and valid SNPs)
        will be stored.

    feature_scores : str
        Path to a tab‑separated file where the ranked list of selected
        SNP features and their chi‑square scores will be written.

    Output
    ------
    Two files:
        - A pickle file containing the cleaned genotype matrix restricted
          to the representatives of all clusters and the variety of
          interest, with invalid SNPs removed and ambiguous genotypes
          imputed using the genotype of the variety of interest.
        - A TSV file listing the top 2,000 SNP features selected by the
          chi‑square statistic, sorted by decreasing score and then by
          feature name.

    """
    _check_file_exists(picklefile1, cluster_json)
    # Determine next cultivar ID dynamically from the pickle
    with open(picklefile1, 'rb') as f:
        X = pickle.load(f)
    X["VARIETY"] = pd.to_numeric(X["VARIETY"])
    max_id = X["VARIETY"].max()
    variety_id = int(max_id)
    X["VARIETY"] = X["VARIETY"].astype(int)
    # load clusters json
    with open(cluster_json, "r") as f:
        clusters = json.load(f)
    # determine representatives of all clusters
    representatives = {}
    cluster_of_interest = None
    members_of_interest = None
    for cl, ids in clusters.items():
        ids_int = [int(i) for i in ids]
        representatives[cl] = ids_int[0]
        if variety_id in ids_int:
            cluster_of_interest = cl
            members_of_interest = ids_int
    if cluster_of_interest is None:
        raise ValueError(f"Variety {variety_id} not found in any cluster")
    allowed_varieties = set()
    for cl, rep in representatives.items():
        if cl == cluster_of_interest:
            allowed_varieties.add(variety_id)  # keep only the target variety from its cluster
        else:
            allowed_varieties.add(rep)  # keep only the representative from other clusters
    X = X[X["VARIETY"].isin(allowed_varieties)]
    idx = X[X["VARIETY"] == variety_id].index[0]
    # Save the identifier column before dropping it
    variety_series = X["VARIETY"].copy()
    # remove the identifier column because this is not a feature
    X.drop(['VARIETY'], axis=1, inplace=True)
    # Define valid genotype values
    valid_genotypes = {'A', 'C', 'G', 'T'}
    # Identify columns where the value at index 'idx' is not valid
    invalid_columns = [col for col in X.columns if X.at[idx, col] not in valid_genotypes]
    # Drop invalid columns from the DataFrame
    X.drop(columns=invalid_columns, inplace=True)
    # preprocess missing/ambiguous data: for all ? or degenerate bases, set value to that of the variety of interest
    for (colname, column) in X.items():
        ref_val = X.at[idx, colname]
        column.where(column.isin(['A', 'C', 'G', 'T']), other=ref_val, inplace=True)
    # Restore the VARIETY column using the original index
    X["VARIETY"] = variety_series
    # store features as pickle
    with open(picklefile2, 'wb') as f:
        pickle.dump(X, f)
    X.drop(['VARIETY'], axis=1, inplace=True)
    # save column names for later
    feature_names = X.columns.values
    # generate array for the label of the observations
    y = ["N"] * len(X)
    y = pd.Series("N", index=X.index)
    y[idx] = "Y"  # set the value for the variety of interest to yes
    # transform input data
    oe = OrdinalEncoder()
    oe.fit(X)
    X_trans = oe.transform(X)
    le = LabelEncoder()
    le.fit(y)
    y_trans = le.transform(y)
    if X_trans.shape[1] == 0:
        print("No features available in test data")
        # Create an empty dataframe with the same columns
        ns_df_sorted = pd.DataFrame(columns=['Feat_names', 'Scores'])
        # Write empty file with header
        ns_df_sorted.to_csv(feature_scores, sep="\t", index=False)
    else:
        fs = SelectKBest(score_func=chi2, k=2000)
        fs.fit(X_trans, y_trans)
        # New dataframe with the selected features
        names = feature_names[fs.get_support()]
        scores = fs.scores_[fs.get_support()]
        names_scores = list(zip(names, scores))
        ns_df = pd.DataFrame(data=names_scores, columns=['Feat_names', 'Scores'])
        # Sort the dataframe
        ns_df_sorted = ns_df.sort_values(['Scores', 'Feat_names'], ascending=[False, True])
        # write to file
        ns_df_sorted.to_csv(feature_scores, sep="\t")

def add_info_to_indices(Top_features: str, SNP_INFO: str, Features_with_annotation: str) -> None:
    """
    Merge the list of top‑ranked SNP features with an external SNP
    annotation table and write the combined feature‑annotation mapping
    to a tab‑separated output file.

    Parameters
    ----------
    Top_features : str
        Path to a TSV file containing the ranked SNP features produced
        by the feature‑selection step. The file is expected to contain
        feature names as the index column.

    SNP_INFO : str
        Path to a TSV file containing SNP annotation metadata, indexed
        by SNP identifier. These annotations are merged onto the feature
        list using the shared index.

    Features_with_annotation : str
        Path to the output TSV file where the merged feature list and
        corresponding annotation fields will be written.

    Output
    ------
    TSV file in which:
        - Each selected SNP feature is retained as an index entry.
        - Annotation fields from the SNP metadata table are appended
          using a left join on the feature index.
        - The resulting annotated feature table is written with tab
          separation and index preservation.
    """
    _check_file_exists(Top_features, SNP_INFO)
    top_features = pd.read_csv(Top_features, sep='\t', index_col=1)
    snp_info = pd.read_csv(SNP_INFO, sep='\t', index_col=0)
    merged = top_features.merge(snp_info, how='left', left_index=True, right_index=True)
    merged.to_csv(Features_with_annotation, sep='\t', index=True)

def list_remaining_accessions(picklefile2: str, Features_with_annotation: str, Features_with_annotation_withremainingaccession: str) -> None:
    """
    Compute per‑SNP genotype matching statistics
    across all remaining accessions, and append these metrics to the
    annotated feature table.

    Parameters
    ----------
    picklefile2 : str
        Path to a pickle file containing the preprocessed genotype
        DataFrame restricted to the representatives of all clusters and
        the variety of interest. The DataFrame must include a 'VARIETY'
        column followed by SNP genotype columns.

    Features_with_annotation : str
        Path to a TSV file containing the selected SNP features together
        with their annotation fields.

    Features_with_annotation_withremainingaccession : str
        Path to the output TSV file where the annotated feature table
        will be extended with two additional fields: the number of
        accessions matching the genotype of the variety of interest for
        each SNP, and the set of distinct genotypes observed at that SNP.

    Output
    ------
    TSV file in which:
        - Each selected SNP feature is annotated with the number of
          accessions whose genotype matches that of the variety of
          interest.
        - A list of all distinct genotypes observed at each SNP is
          included for downstream inspection.
    """
    _check_file_exists(picklefile2, Features_with_annotation)
    # retrieve object from pickle
    with open(picklefile2, 'rb') as f:
        data = pickle.load(f)
    features = pd.read_csv(Features_with_annotation, sep='\t')
    # indices of 2000 best ranking snps, in order
    features = features.astype({'Feat_names': 'str'})
    snps = features['Feat_names'].values.tolist()
    # test how many accessions have the same genotype as the variety of interest for each SNP
    remainingaccessions = []
    genotypes = []
    # Determine next cultivar ID dynamically from the pickle
    data["VARIETY"] = pd.to_numeric(data["VARIETY"])
    max_id = data["VARIETY"].max()
    variety_id = int(max_id)
    idx = data[data["VARIETY"] == variety_id].index[0]
    data = data[data["VARIETY"] != 238].copy()
    data.drop(['VARIETY'], axis=1, inplace=True)
    snps_df = data[snps]
    # Row of interest
    row = snps_df.loc[idx].to_numpy()
    # Convert to NumPy
    values = snps_df.to_numpy()
    genotypes = [np.unique(values[:, i]) for i in range(values.shape[1])]
    # remainingaccessions: count matching genotypes per SNP
    remainingaccessions = np.sum(values == row, axis=0).tolist()
    features['remainingaccessions'] = remainingaccessions
    features['genotypes'] = genotypes
    features.to_csv(Features_with_annotation_withremainingaccession, sep='\t', index=True)

def screen_2_snp_combinations(picklefile: str, Features_with_annotation: str, max_fing:int, output: str) -> None:
    """
    Generate 2‑SNV barcodes for a target cultivar by scanning all
    cross‑chromosomal SNP pairs, identifying pairs that uniquely match
    the genotype of the variety of interest, and writing the resulting
    fingerprint list to an output file.

    Parameters
    ----------
    picklefile : str
        Path to a pickle file containing the full genotype matrix with a
        'VARIETY' column followed by SNP genotype columns. The variety of
        interest must be present in this matrix.

    Features_with_annotation : str
        Path to a TSV file containing the selected SNP features together
        with their annotation fields, including a 'Feat_names' column and
        a 'CHROMOSOME' column used to enforce cross‑chromosomal pairing.

    output : str
        Path to the output file where the list of accepted 2‑SNV
        barcodes will be written, one barcode per line.

    Output
    ------
    TSV file in which:
        - Each line contains a tab‑separated SNP pair that uniquely
          identifies the variety of interest across all accessions.
        - Only SNP pairs located on different chromosomes are considered.
        - SNP usage is balanced by gradually increasing the maximum
          number of times each SNP may appear in accepted fingerprints.
        - If no valid fingerprints are detected, a diagnostic message is
          written instead.

    """
    _check_file_exists(picklefile, Features_with_annotation)
    # Determine next cultivar ID dynamically from the pickle
    with open(picklefile, 'rb') as f:
        data = pickle.load(f)
    data["VARIETY"] = pd.to_numeric(data["VARIETY"])
    max_id = data["VARIETY"].max()
    variety_id = int(max_id)
    data["VARIETY"] = data["VARIETY"].astype(int)
    data = data[data["VARIETY"] != 238].copy()
    idx = data[data["VARIETY"] == variety_id].index[0]
    data.drop(['VARIETY'], axis=1, inplace=True)
    # Load features
    features = pd.read_csv(Features_with_annotation, sep='\t')
    features['Feat_names'] = features['Feat_names'].astype(str)
    snps = features['Feat_names'].values
    n_snps = len(snps)
    # Pre-map SNP --> chromosome
    chroms = features['CHROMOSOME'].values
    # Convert genotype matrix to NumPy for speed
    data_np = data[snps].to_numpy()
    rows = data.index.to_numpy()
    # Precompute genotype of variety of interest
    voi_genotype = data.loc[idx, snps].to_numpy()
    # Precompute matches for each SNP (boolean matrix)
    valid_match_matrix = (data_np == voi_genotype)
    # Precompute chromosome inequality mask for all SNP pairs
    chrom_diff = chroms[:, None] != chroms[None, :]
    # Get all SNP pairs where chromosomes differ
    pair_indices = np.transpose(np.nonzero(np.triu(chrom_diff, k=1)))
    # Track how many times each SNP is used
    usage_count = np.zeros(n_snps, dtype=int)
    fingerprints = []
    max_fingerprints = max_fing
    # We will increase the allowed usage limit gradually
    limit = 1
    # Continue until we reach 100 fingerprints or exhaust all limits
    while len(fingerprints) < max_fingerprints and limit <= n_snps:
        for i, j in pair_indices:
            # Skip if either SNP already reached the usage limit
            if usage_count[i] >= limit or usage_count[j] >= limit:
                continue
            matches = valid_match_matrix[:, i] & valid_match_matrix[:, j]
            if np.sum(matches) == 1:
                # Accept this fingerprint
                pair = "\t".join(sorted([snps[i], snps[j]]))
                if pair not in fingerprints:
                    fingerprints.append(pair)
                # Update usage counts
                usage_count[i] += 1
                usage_count[j] += 1
                if len(fingerprints) == max_fingerprints:
                    break
        # Increase the allowed usage limit
        limit += 1
    with open(output, 'w') as outfile:
        if fingerprints:
            outfile.write("\n".join(fingerprints) + "\n")
        else:
            outfile.write("#No Fingerprints was detected\n")

def add_info_to_snp_combinations(Combinations: str, Scores: str, SNP_INFO: str, picklefile: str, config: str, output: str) -> None:
    """
    Annotate 2‑SNV barcode combinations with SNP metadata, genotype
    information for the variety of interest, and feature‑selection scores,
    producing a fully annotated fingerprint table for downstream analysis.

    Parameters
    ----------
    Combinations : str
        Path to a TSV file containing the list of SNP index pairs
        representing candidate 2‑SNV barcodes. Each row must contain
        two integer indices ('index1', 'index2').

    Scores : str
        Path to a TSV file containing feature‑selection scores (chi‑square statistics)
        indexed by SNP identifier.

    SNP_INFO : str
        Path to a TSV file containing SNP annotation metadata indexed by
        SNP identifier.

    picklefile : str
        Path to a pickle file containing the full genotype matrix with a
        'VARIETY' column followed by SNP genotype columns.

    output : str
        Specifying the path where the final annotated fingerprint table
        will be written.

    Output
    ------
    TSV file in which:
        - If the fingerprint file contains the message '#No Fingerprints
          was detected', this message is propagated directly to the output.
        - Otherwise, each SNP pair is annotated with:
              • SNP metadata for both SNPs,
              • the genotype of the variety of interest at both SNPs,
              • feature‑selection scores for both SNPs.
    """
    _check_file_exists(Combinations, Scores, SNP_INFO, picklefile)
    # Read the first line of the fingerprint file to check if there are any fingerprints
    with open(Combinations, 'r') as fing:
        first_line = fing.readline().strip()
    # Check if the file contains the specific message
    if first_line.startswith('#No'):
        # Write the message to the output file
        with open(output, 'w') as out_f:
            out_f.write('#No Fingerprints was detected\n')
        results = {"Number of generated 2-SNV barcodes": "No barcode could be generated for this cultivar"}
        config_path = Path(config)
        # Load existing YAML if file exists, otherwise start empty
        if config_path.exists():
            cfg = yaml.safe_load(config_path.read_text()) or {}
        else:
            cfg = {}
        # Remove old keys if they exist
        for key in results:
            cfg.pop(key, None)
        # Add new values
        cfg.update(results)
        # Write updated YAML (creates file if missing)
        config_path.write_text(yaml.dump(cfg, sort_keys=False))
    else:
        combinations = pd.read_csv(Combinations, sep='\t', names=['index1', 'index2'])
        scores = pd.read_csv(Scores, sep='\t', index_col=0)
        snp_info = pd.read_csv(SNP_INFO, sep='\t', index_col=0)
        # Determine next cultivar ID dynamically from the pickle
        with open(picklefile, 'rb') as f:
            data = pickle.load(f)
        data["VARIETY"] = pd.to_numeric(data["VARIETY"])
        max_id = data["VARIETY"].max()
        variety_id = int(max_id)
        data["VARIETY"] = data["VARIETY"].astype(int)
        idx = data[data["VARIETY"] == variety_id].index[0]
        data.drop(['VARIETY'], axis=1, inplace=True)
        genotype = data.loc[idx,]
        genotype = genotype.reset_index()
        genotype = genotype.astype({'index': 'int'})
        genotype = genotype.rename(columns={idx: 'GENOTYPE'})
        merged_1 = combinations.merge(snp_info, how='left', left_on='index1', right_index=True)
        merged_genotype_1 = merged_1.merge(genotype, how='left', left_on='index1', right_on='index')
        merged_genotype_1 = merged_genotype_1.drop(columns='index')
        merged_2 = merged_genotype_1.merge(scores, how='left', left_on='index1', right_on='Feat_names')
        merged_3 = merged_2.merge(snp_info, how='left', left_on='index2', right_index=True)
        merged_genotype_2 = merged_3.merge(genotype, how='left', left_on='index2', right_on='index')
        merged_genotype_2 = merged_genotype_2.drop(columns='index')
        merged_4 = merged_genotype_2.merge(scores, how='left', left_on='index2', right_on='Feat_names')
        merged_4.to_csv(output, sep='\t', index=True)
        results = {"Number of generated 2-SNV barcodes": str(len(merged_4))}
        config_path = Path(config)
        # Load existing YAML if file exists, otherwise start empty
        if config_path.exists():
            cfg = yaml.safe_load(config_path.read_text()) or {}
        else:
            cfg = {}
        # Remove old keys if they exist
        for key in results:
            cfg.pop(key, None)
        # Add new values
        cfg.update(results)
        # Write updated YAML (creates file if missing)
        config_path.write_text(yaml.dump(cfg, sort_keys=False))


def unique_SNPs(Fingerprints: str, Features_with_annotation_withremainingaccession: str, config: str, output: str) -> None:
    """
    Summarize all SNPs participating in 2‑SNV barcodes by extracting
    per‑SNP metadata, genotype information, and feature‑selection scores.

    Parameters
    ----------
    Fingerprints : str
        Path to a TSV file containing annotated 2‑SNV barcode pairs.
        Each row must include SNP indices ('index1', 'index2') together
        with annotation fields and feature‑selection scores for both SNPs.

    Features_with_annotation_withremainingaccession : str
        Path to a TSV file containing SNP‑level annotation fields and the
        number of accessions matching the genotype of the variety of
        interest for each SNP.

    output : str
        Path to the final TSV file where the non‑redundant, per‑SNP
        summary table will be written.

    Output
    ------
    TSV file in which:
        - If the fingerprint file contains the message '#No Fingerprints
          was detected', this message is propagated directly to the output.
        - Otherwise, all SNPs appearing in any fingerprint pair are
          extracted, deduplicated, and sorted by decreasing feature‑
          selection score.
        - For each SNP, the output includes:
              • chromosome, position, reference allele,
              • genotype of the variety of interest,
              • feature‑selection score,
              • the set of genotypes observed across accessions.
    """
    _check_file_exists(Fingerprints, Features_with_annotation_withremainingaccession)
    # Read the first line of the fingerprint file to check if there are any fingerprints
    with open(Fingerprints, 'r') as fing:
        first_line = fing.readline().strip()
    # Check if the file contains the specific message
    if first_line.startswith('#No'):
        # Write the message to the output file
        with open(output, 'w') as out_f:
            out_f.write('#No Fingerprints was detected\n')
        results = {"Count of unique SNVs comprising the fingerprint": "No SNV could be generated for this cultivar"}
        config_path = Path(config)
        # Load existing YAML if file exists, otherwise start empty
        if config_path.exists():
            cfg = yaml.safe_load(config_path.read_text()) or {}
        else:
            cfg = {}
        # Remove old keys if they exist
        for key in results:
            cfg.pop(key, None)
        # Add new values
        cfg.update(results)
        # Write updated YAML (creates file if missing)
        config_path.write_text(yaml.dump(cfg, sort_keys=False))
    else:
        df = pd.read_csv(Fingerprints, sep='\t')
        df_index1 = df[['index1', 'CHROMOSOME_x', 'POSITION_x', 'REFCALL_x', 'GENOTYPE_x', 'Scores_x']].copy()
        df_index1.columns = ['snp_index', 'chromosome', 'position', 'refcall', 'genotype', 'score']
        df_index2 = df[['index2', 'CHROMOSOME_y', 'POSITION_y', 'REFCALL_y', 'GENOTYPE_y', 'Scores_y']].copy()
        df_index2.columns = ['snp_index', 'chromosome', 'position', 'refcall', 'genotype', 'score']
        df_combined = pd.concat([df_index1, df_index2], ignore_index=True)
        df_unique = df_combined.drop_duplicates(subset='snp_index')
        df_unique = df_unique.sort_values(by='score', ascending=False).reset_index(drop=True)
        features = pd.read_csv(Features_with_annotation_withremainingaccession, sep='\t')
        features = features[['Feat_names', 'genotypes']]
        df_final = df_unique.merge(features, how='left', left_on='snp_index', right_on='Feat_names')
        df_final = df_final.drop(columns='Feat_names')
        df_final.to_csv(output, sep='\t', index=False)
        results = {"Count of unique SNVs comprising the fingerprint": str(len(df_final))}
        config_path = Path(config)
        # Load existing YAML if file exists, otherwise start empty
        if config_path.exists():
            cfg = yaml.safe_load(config_path.read_text()) or {}
        else:
            cfg = {}
        # Remove old keys if they exist
        for key in results:
            cfg.pop(key, None)
        # Add new values
        cfg.update(results)
        # Write updated YAML (creates file if missing)
        config_path.write_text(yaml.dump(cfg, sort_keys=False))

def bed_maker(input: str, output: str) -> None:
    """
    Convert a SNP annotation table into BED‑like file suitable
    for downstream liftover and variant processing.

    Parameters
    ----------
    input : str
        Path to a tab‑separated file containing SNP records with at least
        the fields 'chromosome', 'position', 'refcall', 'genotype', and
        'snp_index'. Each row represents a single SNP.

    output : str
        Path to the BED‑formatted file to be written.

    Output
    ------
    BED‑like file with one line per SNP.
    All coordinates follow UCSC BED conventions
    to ensure compatibility with liftover tools.
    """
    _check_file_exists(input)
    with open(input) as tsv, open(output, "w") as bed:
        reader = csv.DictReader(tsv, delimiter="\t")
        for row in reader:
            chrom_num = row["chromosome"]
            pos = int(row["position"])
            ref = row["refcall"]
            alt = row["genotype"]
            # liftover requires chr1, chr2, ...
            chrom = f"chr{chrom_num}"
            # BED is 0-based start, 1-based end
            start = pos - 1
            end = pos
            snp_id = row["snp_index"]
            bed.write(f"{chrom}\t{start}\t{end}\t{snp_id}\t{ref}\t{alt}\n")

def vcf_maker(input: str, fasta_path: str, output: str) -> None:
    """
    Convert a BED‑like SNP file produced after liftover into a sorted
    VCFv4.2 file by mapping UCSC chromosome labels to reference contig
    identifiers, reconstructing 1‑based variant positions, and emitting a
    minimal VCF with contig metadata and variant records.

    Parameters
    ----------
    input : str
        Path to a BED‑like file containing SNPs after liftover.
        Each row must contain:
            chrN, start, end, snp_index, ref, alt
        where coordinates follow UCSC BED conventions (0‑based start,
        1‑based end).

    output : str
        Path to the VCF file to be written. The output will contain a
        VCFv4.2 header, contig declarations for all chromosomes, and one
        variant record per SNP.

    Output
    ------
    VCFv4.2 file in which:
        - UCSC chromosome labels (e.g., 'chr1') are mapped to their
          corresponding reference contig identifiers (e.g., 'NC_089035.1').
        - Variant positions are reconstructed using the BED end coordinate
          as the 1‑based POS field.
        - All records are sorted by contig order and genomic position.
        - Each variant line contains:
              CHROM, POS, ID, REF, ALT, QUAL='.', FILTER='.', INFO='.'
        - The resulting file is fully compatible with RiSpy for line identification and downstream analysis.
    """
    _check_file_exists(input, fasta_path)
    records = []
    REF = pysam.FastaFile(fasta_path)
    CHR_MAP = {}
    CONTIGS = []
    # Take only the first 12 contigs in FASTA order
    selected_contigs = REF.references[:12]
    for i, chrom in enumerate(selected_contigs, start=1):
        length = REF.get_reference_length(chrom)
        CHR_MAP[str(i)] = chrom
        CONTIGS.append((chrom, length))
    CONTIG_ORDER = {c[0]: i for i, c in enumerate(CONTIGS)}
    with open(input) as bed:
        reader = csv.reader(bed, delimiter="\t")
        for row in reader:
            chrom_chr, start, end, snp_id, ref, alt = row
            chrom_num = chrom_chr.replace("chr", "")
            chrom = CHR_MAP.get(chrom_num)
            pos = int(end)
            records.append((chrom, pos, snp_id, ref, alt))
    # Sort by contig order + position
    records.sort(key=lambda x: (CONTIG_ORDER.get(x[0], 999), x[1]))
    updated_records = []
    for chrom, pos, snp_id, ref_old, alt in records:
        try:
            base = REF.fetch(chrom, pos - 1, pos).upper()
            if base not in ("A", "C", "G", "T"):
                base = "N"
        except Exception:
            base = "N"
        # Append updated tuple with corrected REF
        updated_records.append((chrom, pos, snp_id, base, alt))
    with open(output, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        for contig, length in CONTIGS:
            vcf.write(f"##contig=<ID={contig},length={length}>\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, snp_id, ref, alt in updated_records:
            if ref == alt:
                continue  # skip non-variant sites
            vcf.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t.\t.\t.\n")