#!/usr/bin/env python3

import os
import re
import pandas as pd
import numpy as np
import logging
import time
import tempfile
from pathlib import Path
from phylodm import PhyloDM
from Bio import SeqIO
from Bio import Phylo
from typing import Dict, Set, List, Tuple

# Set up logging
logger = logging.getLogger("lasvdedup.duplicates")

def setup_logging(level=logging.INFO):
    """Configure logging for the duplicate detection module."""
    logger.setLevel(level)

    # Create console handler if none exists
    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

def to_distance_matrix(tree: PhyloDM) -> Tuple[List, np.ndarray]:
    """Create a distance matrix (NumPy array) from clades/branches in tree.

    A cell (i,j) in the array is the path distance between labels[i]
    and labels[j], as calculated from the phylogenetic tree.

    Returns a tuple of (labels, distance_matrix) where labels is a list of
    terminal clades and distance_matrix is a NumPy array of distances.
    """
    logger.info("Calculating distance matrix from phylogenetic tree")
    start_time = time.time()

    dm = tree.dm(norm=False)
    labels = [l.replace(" ", "_") for l in tree.taxa()]

    elapsed = time.time() - start_time
    logger.info(f"Distance matrix calculation completed in {elapsed:.2f} seconds")
    return (labels, dm)

def write_distance_matrix(dist_matrix: np.ndarray, labels: List, output_path: str) -> None:
    """
    Write distance matrix to a file.

    Args:
        dist_matrix: Distance matrix as numpy array
        labels: List of terminal labels (strings) or objects with name attribute
        output_path: Path to write the distance matrix
    """
    logger.info(f"Writing distance matrix to {output_path}")
    n_labels = len(labels)

    with open(output_path, "w") as mldist_file:
        mldist_file.write(f"{n_labels}\n")
        for i, term in enumerate(labels):
            # Handle both string labels and objects with name attribute
            label = term if isinstance(term, str) else term.name
            line = [label]
            for j in range(n_labels):
                line.append(f"{dist_matrix[i, j]:.6f}")
            mldist_file.write("\t".join(line) + "\n")
    logger.info("Distance matrix written successfully")

def load_read_counts(table_path: Path, reads_column:str) -> Dict[str, Dict[str, int]]:
    """
    Extract contig ID to read count mapping from contigs table.

    Args:
        table_path: Path to the contigs table

    Returns:
        Dictionary mapping contig IDs to read counts
    """
    logger.info(f"Loading read counts from {table_path}")
    try:
        contigs_df = pd.read_csv(str(table_path), sep="\t")
        logger.info(f"Loaded table with {len(contigs_df)} rows and {len(contigs_df.columns)} columns")

        if reads_column in contigs_df.columns:
            logger.info(f"Using specified reads column: {reads_column}")
            contigs_df["reads"] = contigs_df[reads_column]
        elif "(samtools Post-dedup) reads mapped (R1+R2)" in contigs_df.columns:
            logger.info("Using '(samtools Post-dedup) reads mapped (R1+R2)' column")
            contigs_df["reads"] = contigs_df["(samtools Post-dedup) reads mapped (R1+R2)"]
        elif "(samtools Raw) reads mapped (R1+R2)" in contigs_df.columns:
            logger.info("Using '(samtools Raw) reads mapped (R1+R2)' column")
            contigs_df["reads"] = contigs_df["(samtools Raw) reads mapped (R1+R2)"]
        else:
            logger.error(f"Reads column {reads_column} not found in table")
            raise ValueError(f"Reads column {reads_column} not found in table: {table_path}")

        contigs_df.set_index("index", inplace=True)
        contigs_df = contigs_df[["reads"]]
        contig_to_reads = contigs_df.to_dict('index')
        logger.info(f"Processed read counts for {len(contig_to_reads)} contigs")
        return contig_to_reads

    except Exception as e:
        logger.error(f"Error loading read counts: {e}", exc_info=True)
        raise

def group_sequences_by_sample(tips: List, sample_regex: str) -> Dict[str, List[str]]:
    """
    Group sequences by sample ID using regex.

    Args:
        tips: List of tree tips
        sample_regex: Regular expression to extract sample identifiers

    Returns:
        Dictionary mapping sample IDs to lists of sequence names
    """
    logger.info(f"Grouping sequences by sample using regex: {sample_regex}")
    sample_to_seqs = {}
    pattern = re.compile(sample_regex)
    no_match_count = 0

    for t in tips:
        match = pattern.search(t)
        if match:
            sample_id = match.group(0)
            if sample_id not in sample_to_seqs:
                sample_to_seqs[sample_id] = []
            sample_to_seqs[sample_id].append(t)
        else:
            no_match_count += 1

    if no_match_count > 0:
        logger.warning(f"No sample ID match for {no_match_count} sequences")

    logger.info(f"Found {len(sample_to_seqs)} unique sample IDs")
    # Log stats on sequences per sample
    seqs_per_sample = [len(seqs) for seqs in sample_to_seqs.values()]
    if seqs_per_sample:
        logger.info(f"Sequences per sample - Min: {min(seqs_per_sample)}, Max: {max(seqs_per_sample)}, "
                   f"Avg: {sum(seqs_per_sample)/len(seqs_per_sample):.2f}")

    return sample_to_seqs

def find_duplicates(
    sample_to_seqs: Dict[str, List[str]],
    tips: List,
    dist_matrix: np.ndarray,
    contig_to_reads: Dict[str, Dict[str, int]],
    threshold: float,
) -> Dict[str, str]:
    """
    Process each sample group to identify duplicates.

    Args:
        sample_to_seqs: Dictionary mapping sample IDs to lists of sequence names
        dist_matrix: Distance matrix as numpy array
        contig_to_reads: Dictionary mapping contig IDs to read counts
        threshold: Distance threshold to identify duplicates

    Returns:
        Dictionary mapping each tip name to its classification (good, bad, or comorbidity)
    """
    logger.info(f"Finding duplicates with distance threshold: {threshold}")
    classifications = {}  # Track classifications for each tip
    tips_lookup = {t: i for i, t in enumerate(tips)}

    # Count statistics
    good_count = 0
    bad_count = 0
    comorbidity_count = 0

    total_samples = len(sample_to_seqs)
    for i, (sample, seq_names) in enumerate(sample_to_seqs.items()):
        logger.debug(f"Processing sample {i+1}/{total_samples}: {sample} with {len(seq_names)} sequences")

        # Skip if only one sequence for this sample
        if len(seq_names) <= 1:
            for seq in seq_names:
                classifications[seq] = "good"
                good_count += 1
            continue

        distances = get_distances(seq_names, tips_lookup, dist_matrix)

        # Consensus sequences are all similar
        if all(d <= threshold for d in distances):
            logger.debug(f"Sample {sample}: All sequences are similar (below threshold)")
            reads = {seq: get_reads(seq, contig_to_reads) for seq in seq_names}
            # Keep sequence with highest read count
            max_reads_seq = max(reads, key=reads.get)
            classifications[max_reads_seq] = "good"
            good_count += 1
            for seq in set(seq_names) - {max_reads_seq}:
                classifications[seq] = "bad"
                bad_count += 1
        # Consensus sequences are distinct
        else:
            logger.debug(f"Sample {sample}: Sequences are distinct (potential comorbidity)")
            for seq in seq_names:
                classifications[seq] = "comorbidity"
                comorbidity_count += 1

    logger.info(f"Classification results: {good_count} good, {bad_count} bad, {comorbidity_count} comorbidity")
    return classifications

def get_distances(names: List , tips_lookup:Dict, matrix:np.ndarray) -> Dict[str, Dict[str, float]]:
    """
    Get distances between all pairs of sequences in the tree.

    Args:
        names: List of sequence names
        tips_lookup: Dictionary of tree tips with location in matrix
        matrix: Distance matrix as numpy array

    Returns:
        Dictionary mapping pairs of sequence names to distances
    """
    distances = []
    for i, name1 in enumerate(names):
        for j, name2 in enumerate(names):
            if i >= j:
                continue
            idx1 = tips_lookup[name1]
            idx2 = tips_lookup[name2]
            dist = matrix[idx1, idx2]
            distances.append(dist)
    return distances

def get_reads (seq_name: str, contig_to_reads: Dict[str, Dict[str, int]]) -> int:
    """
    Get the number of reads for a given sequence.

    Args:
        seq_name: Sequence name
        contig_to_reads: Dictionary mapping contig IDs to read counts
    """
    if seq_name in contig_to_reads:
        return contig_to_reads[seq_name]["reads"]
    elif seq_name.replace('_R_', '').split(".", 1)[0] in contig_to_reads:
        return contig_to_reads[seq_name.replace('_R_', '').split(".", 1)[0]]["reads"]

    logger.warning(f"Read count not found for sequence: {seq_name}")
    raise ValueError(f"Read count not found for sequence: {seq_name}")


def write_results(
    classifications: Dict[str, str],
    sequences: Dict[str, SeqIO.SeqRecord],
    species: str,
    segment: str,
    prefix: str,
    sample_regex: str
) -> None:
    """
    Write sequences to organized FASTA files and create a classifications summary.

    Args:
        classifications: Dictionary mapping tip names to classifications
        sequences: Dictionary mapping sequence names to SeqIO records
        species: Species name for the output file
        segment: Segment name for the output file
        prefix: Output file prefix (parent directory)
        sample_regex: Regular expression to extract sample identifiers
    """
    logger.info(f"Writing results to {prefix}")
    # Create parent directory if it doesn't exist
    os.makedirs(prefix, exist_ok=True)

    # Compile the sample regex pattern
    pattern = re.compile(sample_regex)

    # Write each sequence to its own file in the appropriate directory
    for seq_name, classification in classifications.items():
        # Extract sample ID from sequence name
        match = pattern.search(seq_name)
        if not match:
            logger.warning(f"Could not extract sample ID from {seq_name}, skipping...")
            continue
        sample_id = match.group(0)

        # Determine effective classification for directory structure
        # (comorbidity sequences go in "good" directory)
        file_classification = "good" if classification in ["good", "comorbidity"] else "bad"

        # Create directory structure
        output_dir = os.path.join(prefix, sample_id, file_classification)
        os.makedirs(output_dir, exist_ok=True)

        # Create filename
        basename = os.path.basename(seq_name)
        output_filename = f"{basename}_{species}_{segment}.fasta"
        output_path = os.path.join(output_dir, output_filename)

        # Write the sequence to a FASTA file
        if seq_name in sequences:
            with open(output_path, "w") as fasta_file:
                SeqIO.write(sequences[seq_name], fasta_file, "fasta")
        else:
            logger.warning(f"Sequence {seq_name} not found in the provided FASTA file, skipping...")

    # Write classification summary file
    with open(f"{prefix}/classifications.tsv", "w") as class_file:
        class_file.write("tip name\tclassification\n")
        for tip_name in sorted(classifications.keys()):
            class_file.write(f"{tip_name}\t{classifications[tip_name]}\n")

def root_tree_at_midpoint(tree_path: Path) -> PhyloDM:
    """
    Root a phylogenetic tree at midpoint and return a PhyloDM object.

    Args:
        tree_path: Path to the unrooted tree file in Newick format

    Returns:
        PhyloDM object containing the rooted tree
    """
    logger.info(f"Loading phylogenetic tree from {tree_path}")

    # Read the tree with Bio.Phylo
    bio_tree = Phylo.read(str(tree_path), "newick")

    # Root the tree at midpoint
    logger.info("Rooting tree at midpoint")
    bio_tree.root_at_midpoint()

    # Write the rooted tree to a temporary file
    with tempfile.NamedTemporaryFile(suffix='.treefile', delete=False) as tmp_tree_file:
        Phylo.write(bio_tree, tmp_tree_file.name, "newick")
        rooted_tree_path = tmp_tree_file.name

    # Load the rooted tree with PhyloDM
    logger.info("Loading rooted tree into PhyloDM")
    tree_obj = PhyloDM.load_from_newick_path(rooted_tree_path)

    # Clean up the temporary file
    try:
        os.unlink(rooted_tree_path)
    except:
        logger.warning(f"Could not remove temporary tree file: {rooted_tree_path}")

    logger.info("Phylogenetic tree loaded and rooted successfully")
    return tree_obj

def determine_duplicates(
        tree: Path,
        sequences: Path,
        prefix: str,
        table: Path,
        sample_regex: str,
        reads_column: str,
        species: str,
        segment: str,
        threshold: float,
        log_level: str
) -> None:
    """
    Main function to detect duplicate sequences in a phylogenetic tree.

    Args:
        tree: Path to the phylogenetic tree file
        sequences: Path to the sequences file (FASTA format)
        prefix: Output file prefix
        table: Path to the table with contig information
        sample_regex: Regular expression to extract sample identifiers
        reads_column: Name of the column containing read counts
        species: Species name for the output files
        segment: Segment name for the output files
        threshold: Distance threshold to identify duplicates
    """
    # Setup logging
    setup_logging(log_level)

    # Log start of the process with parameters
    logger.info("Starting duplicate detection process")
    logger.info(f"Parameters: tree={tree}, sequences={sequences}, prefix={prefix}, table={table}")
    logger.info(f"  sample_regex={sample_regex}, reads_column={reads_column}, threshold={threshold}")
    logger.info(f"  species={species}, segment={segment}")

    start_time = time.time()

    # Create output directory if it doesn't exist
    os.makedirs(prefix, exist_ok=True)

    # Load and root the tree using the modular function
    tree_obj = root_tree_at_midpoint(tree)

    # Load sequences
    logger.info(f"Loading sequences from {sequences}")
    seq_records = SeqIO.to_dict(SeqIO.parse(str(sequences), "fasta"))

    # Calculate distances from tree
    tips, dist_matrix = to_distance_matrix(tree_obj)

    # Write distance matrix to file
    write_distance_matrix(dist_matrix, tips, f"{prefix}/distance_matrix.mldist")

    # Load read counts from contigs table
    contig_to_reads = load_read_counts(table, reads_column)

    # Group sequences by sample
    sample_to_seqs = group_sequences_by_sample(tips, sample_regex)

    # Find duplicates
    classifications = find_duplicates(
        sample_to_seqs, tips, dist_matrix, contig_to_reads, threshold
    )

    # Write results to output files
    write_results(classifications, seq_records, species, segment, prefix, sample_regex)

    # Log completion time
    elapsed = time.time() - start_time
    logger.info(f"Duplicate detection completed in {elapsed:.2f} seconds")
