#!/usr/bin/env python3

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import Phylo
from typing import Dict, Set, List, Tuple


def to_distance_matrix(tree: Phylo.BaseTree.Tree) -> Tuple[List, np.ndarray]:
    """Create a distance matrix (NumPy array) from clades/branches in tree.

    A cell (i,j) in the array is the path distance between terminals[i]
    and terminals[j], as calculated from the phylogenetic tree.

    Returns a tuple of (terminals, distance_matrix) where terminals is a list of
    terminal clades and distance_matrix is a NumPy array of distances.
    """
    terminals = tree.get_terminals()
    n_terms = len(terminals)

    # Initialize matrix with zeros
    distmat = np.zeros((n_terms, n_terms))

    # Fill matrix with pairwise distances
    for i, t1 in enumerate(terminals):
        for j, t2 in enumerate(terminals):
            if i < j:  # Calculate each pair only once
                dist = tree.distance(t1, t2)
                distmat[i, j] = dist
                distmat[j, i] = dist

    return (terminals, distmat)

def write_distance_matrix(dist_matrix: np.ndarray, terminals: List, output_path: str) -> None:
    """
    Write distance matrix to a file.

    Args:
        dist_matrix: Distance matrix as numpy array
        terminals: List of terminal objects
        output_path: Path to write the distance matrix
    """
    n_terminals = len(terminals)

    with open(output_path, "w") as mldist_file:
        mldist_file.write(f"{n_terminals}\n")
        for i, term in enumerate(terminals):
            line = [term.name]
            for j in range(n_terminals):
                line.append(f"{dist_matrix[i, j]:.6f}")
            mldist_file.write("\t".join(line) + "\n")


def load_read_counts(table_path: Path, reads_column:str) -> Dict[str, Dict[str, int]]:
    """
    Extract contig ID to read count mapping from contigs table.

    Args:
        table_path: Path to the contigs table

    Returns:
        Dictionary mapping contig IDs to read counts
    """
    contigs_df = pd.read_csv(str(table_path), sep="\t")

    if reads_column in contigs_df.columns:
        contigs_df["reads"] = contigs_df[reads_column]
    elif "(samtools Post-dedup) reads mapped (R1+R2)" in contigs_df.columns:
        contigs_df["reads"] = contigs_df["(samtools Post-dedup) reads mapped (R1+R2)"]
    elif "(samtools Raw) reads mapped (R1+R2)" in contigs_df.columns:
        contigs_df["reads"] = contigs_df["(samtools Raw) reads mapped (R1+R2)"]
    else:
        raise ValueError(f"Reads column {reads_column} not found in table: {table_path}")

    contigs_df.set_index("index", inplace=True)

    contigs_df = contigs_df[["reads"]]

    contig_to_reads = contigs_df.to_dict('index')
    return contig_to_reads


def group_sequences_by_sample(tips: List, sample_regex: str) -> Dict[str, List[str]]:
    """
    Group sequences by sample ID using regex.

    Args:
        tips: List of tree tips
        sample_regex: Regular expression to extract sample identifiers

    Returns:
        Dictionary mapping sample IDs to lists of sequence names
    """
    sample_to_seqs = {}
    pattern = re.compile(sample_regex)

    for t in tips:
        match = pattern.search(t.name)
        if match:
            sample_id = match.group(0)
            if sample_id not in sample_to_seqs:
                sample_to_seqs[sample_id] = []
            sample_to_seqs[sample_id].append(t.name)

    return sample_to_seqs


def find_duplicates(
    sample_to_seqs: Dict[str, List[str]],
    tips: List,
    dist_matrix: np.ndarray,
    contig_to_reads: Dict[str, int],
    threshold: float
) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Process each sample group to identify duplicates.

    Args:
        sample_to_seqs: Dictionary mapping sample IDs to lists of sequence names
        dist_matrix: Distance matrix as numpy array
        contig_to_reads: Dictionary mapping contig IDs to read counts
        threshold: Distance threshold to identify duplicates

    Returns:
        Tuple containing:
        - Set of good sequences to keep
        - Set of bad sequences to discard
        - Set of sequences with comorbidities
    """
    good_seqs = set()
    bad_seqs = set()
    comorbidities = set()

    tips_lookup = {t.name: i for i, t in enumerate(tips)}

    for _, seq_names in sample_to_seqs.items():
        # Skip if only one sequence for this sample
        if len(seq_names) <= 1:
            good_seqs.update(seq_names)
            continue

        distances = get_distances(seq_names, tips_lookup, dist_matrix)

        # Consensus sequences are all similar
        if all(d <= threshold for d in distances):  # Fixed: Check each distance individually
            reads = {seq: get_reads(seq, contig_to_reads) for seq in seq_names}
            # Keep sequence with highest read count
            max_reads_seq = max(reads, key=reads.get)
            good_seqs.add(max_reads_seq)
            bad_seqs.update(set(seq_names) - {max_reads_seq})
        # Consensus sequences are distinct
        else:
            good_seqs.update(seq_names)
            comorbidities.update(seq_names)

    return good_seqs, bad_seqs, comorbidities

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
    elif seq_name.replace('_R_','') in contig_to_reads:
        return contig_to_reads[seq_name.replace('_R_','')]["reads"]

    raise ValueError(f"Read count not found for sequence: {seq_name}")


def write_results(good_seqs: Set[str], bad_seqs: Set[str], comorbidities: Set[str], prefix: str) -> None:
    """
    Write good and bad sequences to output files.

    Args:
        good_seqs: Set of good sequences to keep
        bad_seqs: Set of bad sequences to discard
        prefix: Output file prefix
    """
    with open(f"{prefix}.good_seqs.txt", "w") as good_file:
        for seq in sorted(good_seqs):
            good_file.write(f"{seq}\n")

    with open(f"{prefix}.bad_seqs.txt", "w") as bad_file:
        for seq in sorted(bad_seqs):
            bad_file.write(f"{seq}\n")

    with open(f"{prefix}.comorbidities.txt", "w") as comorbidities_file:
        for seq in sorted(comorbidities):
            comorbidities_file.write(f"{seq}\n")


def determine_duplicates(
        tree: Path,
        alignment: Path,
        prefix: str,
        table: Path,
        sample_regex: str,
        reads_column: str,
        threshold: float
) -> None:
    """
    Main function to detect duplicate sequences in a phylogenetic tree.

    Args:
        tree: Path to the phylogenetic tree file
        alignment: Path to the alignment file (not used directly)
        prefix: Output file prefix
        table: Path to the table with contig information
        sample_regex: Regular expression to extract sample identifiers
        threshold: Distance threshold to identify duplicates
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(prefix), exist_ok=True)

    if not tree.exists():
        raise FileNotFoundError(f"Tree file not found: {tree}")
    if not table.exists():
        raise FileNotFoundError(f"Contigs table not found: {table}")

    # Load tree
    tree = Phylo.read(str(tree), "newick")

    # Calculate distances from tree
    tips, dist_matrix  = to_distance_matrix(tree)

    # Write distance matrix to file
    write_distance_matrix(dist_matrix, tips, f"{prefix}.mldist")

    # Load read counts from contigs table
    contig_to_reads = load_read_counts(table, reads_column)

    # Group sequences by sample
    sample_to_seqs = group_sequences_by_sample(tips, sample_regex)

    # Find duplicates
    good_seqs, bad_seqs, comorbidities = find_duplicates(
        sample_to_seqs, tips, dist_matrix, contig_to_reads, threshold
    )

    # Write results to output files
    write_results(good_seqs, bad_seqs, comorbidities, prefix)

    ## Visualise tree in baltic