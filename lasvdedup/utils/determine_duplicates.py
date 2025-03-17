#!/usr/bin/env python3

import os
import logging
import time
from Bio import SeqIO
from pathlib import Path
from typing import Dict


# Import sub utility functions
from .tree_utils import root_tree_at_midpoint
from .distance_matrix import to_distance_matrix
from .io_utils import sort_table, write_results, write_distance_matrix
from .sequence_grouping import group_sequences_by_sample, find_duplicates

# Set up logging
logger = logging.getLogger("lasvdedup.duplicates")

def setup_logging(level=logging.INFO, filepath=None):
    """Configure logging for the duplicate detection module."""
    logger.setLevel(level)

    # Create console handler if none exists
    if not logger.handlers:
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        fileHandler = logging.FileHandler(filepath)
        fileHandler.setFormatter(formatter)
        logger.addHandler(fileHandler)

def determine_duplicates(
            config: Dict,  # Pass full config object
            tree: Path = None,
            sequences: Path =None,
            prefix: Path = None,
            table: Path = None,
            segment: str = None,
        ):
    """
    Main function to detect duplicate sequences in a phylogenetic tree.

    Args:
        config: Path to config file or config dictionary (optional)
    """

    # Create output directory if it doesn't exist
    prefix = Path(prefix) or config["prefix"]
    os.makedirs(prefix, exist_ok=True)

    # Extract segment for threshold determination
    if (segment := segment or config.get("segment")) is None:
        raise ValueError("Segment not provided in config or CLI arguments")

    if (species := config.get("species")) is None:
        raise ValueError("Species not provided in config or CLI arguments")

    # Setup logging
    log_level = config.get("log_level", "INFO")
    setup_logging(log_level, f"{prefix}/{species}-{segment}.log")

    # Log start of the process with parameters
    logger.info("Starting duplicate detection process")

    if (thresholds := config.get("DEDUPLICATE", {}).get("THRESHOLDS", {}).get(segment)) is None:
        logger.error(f"No thresholds found for segment {segment}")

    logger.info("Using thresholds for segment %s: %s", segment, thresholds)

    start_time = time.time()

    tree_path = Path(tree) or config["tree"]
    # Load and root the tree using the modular function
    tree = root_tree_at_midpoint(tree_path)

    # Calculate distances from tree
    tips, dist_matrix = to_distance_matrix(tree["phylodm"])

    # Write distance matrix to file
    write_distance_matrix(dist_matrix, tips, f"{prefix}/distance_matrix.mldist")

    # Load read counts and coverage from contigs table
    if (length_column := config.get("length_column")) is None:
        logger.error("Length column not found in config or CLI arguments")
    selection_columns = config.get("selection_column") or []

    # Create dictionary of all contigs with their rank {sample: {rank: x, ...}}
    contig_table = Path(table) or config["table"]
    contigs_ranked = sort_table(contig_table, length_column,
        selection_columns, expected_length = thresholds["target_length"])

    sample_regex = config.get("sample_regex", r'(LVE\d+)_.*')

    # Group sequences by sample
    sample_to_seqs = group_sequences_by_sample(tips, sample_regex)

    # Find duplicates with the determined thresholds - pass segment parameter
    classifications = find_duplicates(
        sample_to_seqs, tips, dist_matrix, contigs_ranked,
        tree["biophylo"], segment, thresholds
        )

    # Write results to output files
    sequence_path = Path(sequences) or config["sequences"]
    logger.info("Loading sequences from %s", sequence_path)
    seq_records = SeqIO.to_dict(SeqIO.parse(str(sequence_path), "fasta"))
    write_results(classifications, seq_records, species, segment, prefix, sample_regex)

    # Log completion time
    elapsed = time.time() - start_time
    logger.info("Duplicate detection completed in %.2f seconds", elapsed)

    # Return classifications for testing purposes
    return classifications
