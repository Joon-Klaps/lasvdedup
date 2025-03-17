import os
import re
import numpy as np
import pandas as pd
import logging
from pathlib import Path
from Bio import SeqIO
from typing import Dict, List

# Import the Classification classes
from .classification import Classification

logger = logging.getLogger("lasvdedup.duplicates")

def sort_table(table_path: Path, length_column:str,
            selection_columns:List[str] = None, expected_length=float) -> pd.DataFrame:
    """
    Extract contig ID to read count, coverage, and sequence length mapping from contigs table.

    Args:
        table_path: Path to the contigs table
        length_column: Name of the column containing length information
        selection_columns: List of columns used for selection
        expected_length: Target length to filter contigs

    Returns:
        DataFrame with contig IDs as index, containing reads, coverage and qlen columns
    """
    logger.info("Loading read counts and coverage from %s", table_path)

    contigs_df = pd.read_csv(str(table_path), sep="\t")
    logger.info("Loaded table with %d rows and %d columns",
                len(contigs_df), len(contigs_df.columns))

    # Ensure selection_columns is a list and not None
    if selection_columns is None:
        selection_columns = []

    # Check if all filtering columns are present in the table
    missing_columns = [col for col in selection_columns + [length_column] if col not in contigs_df.columns]
    if missing_columns:
        logger.error("Selection column(s) %s not found in table: %s", missing_columns, table_path)
        raise ValueError(f"Selection column(s) {missing_columns} not found in table: {table_path}")

    # Calculate distance to expected length
    contigs_df["distance_to_expectation"] = abs(contigs_df[length_column] - expected_length)

    # Filter contigs by selection columns
    rank_columns = ["distance_to_expectation"] + selection_columns

    # Set first column (distance) as ascending=True and all others as ascending=False
    ascending_values = [True] + [False] * len(selection_columns)
    contigs_ranked = contigs_df.sort_values(by=rank_columns, ascending=ascending_values)

    # Have their order written to a column (ie rank)
    contigs_ranked["rank"] = range(1, len(contigs_ranked) + 1)

    logger.info("Sorted contigs by distance to expected length and selection columns")

    # Make sure 'index' column exists before setting it as index
    if 'index' in contigs_ranked.columns:
        contigs_ranked.set_index("index", inplace=True)

    # Select columns for result, be careful not to select columns that don't exist
    result_columns = ['rank'] + rank_columns
    result_df = contigs_ranked[result_columns]

    return result_df.to_dict('index')

def write_distance_matrix(dist_matrix: np.ndarray, labels: List, output_path: str) -> None:
    """Write distance matrix to a file."""
    logger.info("Writing distance matrix to %s", output_path)
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

def write_results(
    classifications: Dict[str, Classification],
    sequences: Dict[str, SeqIO.SeqRecord],
    species: str,
    segment: str,
    prefix: str,
    sample_regex: str
) -> None:
    """
    Write sequences to organized FASTA files and create a detailed classifications summary.

    Args:
        classifications: Dictionary mapping tip names to Classification objects
        sequences: Dictionary mapping sequence names to SeqIO records
        species: Species name for the output file
        segment: Segment name for the output file
        prefix: Output file prefix (parent directory)
        sample_regex: Regular expression to extract sample identifiers
    """
    logger.info("Writing results to %s", prefix)
    # Create parent directory if it doesn't exist
    os.makedirs(prefix, exist_ok=True)

    # Compile the sample regex pattern
    pattern = re.compile(sample_regex)

    # Write each sequence to its own file in the appropriate directory
    for seq_name, classification in classifications.items():
        # Get sample ID from classification or extract from sequence name
        sample_id = classification.sample_id
        if sample_id == "Unknown":
            match = pattern.search(seq_name)
            if not match:
                logger.warning("Could not extract sample ID from %s, skipping...", seq_name)
                continue
            sample_id = match.group(0)

        # Determine directory based on classification type
        file_classification = classification.file_classification

        # Create directory structure
        output_dir = os.path.join(prefix, sample_id, file_classification)
        os.makedirs(output_dir, exist_ok=True)

        # Create filename
        basename = os.path.basename(seq_name)
        output_filename = f"{basename}-{species}-{segment}.fasta"
        output_path = os.path.join(output_dir, output_filename)

        # Write the sequence to a FASTA file
        if seq_name in sequences:
            with open(output_path, "w") as fasta_file:
                SeqIO.write(sequences[seq_name], fasta_file, "fasta")
        else:
            logger.warning("Sequence %s not found in the provided FASTA file, skipping...", seq_name)

    # Write detailed classification summary file
    parent_dir = os.path.dirname(prefix)
    with open(f"{parent_dir}/{species}-{segment}-classifications.tsv", "w") as class_file:
        # Write header
        class_file.write(Classification.header_line() + "\n")

        # Write each classification
        for seq_name in sorted(classifications.keys()):
            class_file.write(classifications[seq_name].to_line() + "\n")

    logger.info("Wrote detailed classification summary to %s/%s-%s-classifications.tsv",
                parent_dir, species, segment)


    with open(f"{parent_dir}/{species}-{segment}.class.figtree.ann", "w") as class_file:
        # Write header
        class_file.write(Classification.annotation_line() + "\n")

        # Write each classification
        for seq_name in sorted(classifications.keys()):
            class_file.write(classifications[seq_name].to_annotation_line() + "\n")

    logger.info("Wrote detailed classification summary to %s/%s-%s.class.figtree.ann",
                parent_dir, species, segment)
