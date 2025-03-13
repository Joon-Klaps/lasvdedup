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

def load_read_counts(table_path: Path, reads_column:str) -> Dict[str, Dict[str, int]]:
    """
    Extract contig ID to read count mapping from contigs table.

    Args:
        table_path: Path to the contigs table

    Returns:
        Dictionary mapping contig IDs to read counts
    """
    logger.info("Loading read counts from %s", table_path)
    try:
        contigs_df = pd.read_csv(str(table_path), sep="\t")
        logger.info("Loaded table with %d rows and %d columns",
                    len(contigs_df), len(contigs_df.columns))

        if reads_column in contigs_df.columns:
            logger.info("Using specified reads column: %s", reads_column)
            contigs_df["reads"] = contigs_df[reads_column]
        elif "(samtools Post-dedup) reads mapped (R1+R2)" in contigs_df.columns:
            logger.info("Using '(samtools Post-dedup) reads mapped (R1+R2)' column")
            contigs_df["reads"] = contigs_df["(samtools Post-dedup) reads mapped (R1+R2)"]
        elif "(samtools Raw) reads mapped (R1+R2)" in contigs_df.columns:
            logger.info("Using '(samtools Raw) reads mapped (R1+R2)' column")
            contigs_df["reads"] = contigs_df["(samtools Raw) reads mapped (R1+R2)"]
        else:
            logger.error("Reads column %s not found in table", reads_column)
            raise ValueError(f"Reads column {reads_column} not found in table: {table_path}")

        contigs_df.set_index("index", inplace=True)
        contigs_df = contigs_df[["reads"]]
        contig_to_reads = contigs_df.to_dict('index')
        logger.info("Processed read counts for %d contigs", len(contig_to_reads))
        return contig_to_reads

    except Exception as e:
        logger.error("Error loading read counts: %s", e, exc_info=True)
        raise

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
