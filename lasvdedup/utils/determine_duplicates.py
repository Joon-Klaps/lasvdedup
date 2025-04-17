#!/usr/bin/env python3

import os
import logging
import time
import yaml
from Bio import SeqIO
from pathlib import Path
from typing import Dict, Union
import pprint

# Import sub utility functions
from .tree_utils import root_tree_at_midpoint
from .distance_matrix import to_distance_matrix
from .io_utils import sort_table, write_results, write_distance_matrix
from .sequence_grouping import group_sequences_by_sample, find_duplicates
from .config_setup import build_config

# Set up logging
logger = logging.getLogger("lasvdedup.duplicates")

def setup_logging(level=logging.INFO, filepath=None):
    """Configure logging for the duplicate detection module."""
    logger.setLevel(level)

    # Prevent propagation to prevent duplicate log messages
    logger.propagate = False

    # Clear existing handlers to avoid duplicates
    if logger.handlers:
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

    # Create console handler
    formatter = logging.Formatter('%(asctime)s - %(filename)s:%(funcName)s - %(levelname)s - %(message)s')
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Add file handler if filepath provided
    if filepath:
        log_dir = os.path.dirname(filepath)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)

        fileHandler = logging.FileHandler(filepath)
        fileHandler.setLevel(level)
        fileHandler.setFormatter(formatter)
        logger.addHandler(fileHandler)

def determine_duplicates(
            config: Union[Dict, str, Path],
            tree: Union[str, Path] = None,
            sequences: Union[str, Path] = None,
            prefix: Union[str, Path] = None,
            table: Union[str, Path] = None,
            segment: str = None,
        ):
    """
    Main function to detect duplicate sequences in a phylogenetic tree.

    Args:
        config: Config dictionary or path to config file
        tree: Override tree path from config
        sequences: Override sequences path from config
        prefix: Override output prefix from config
        table: Override table path from config
        segment: Override segment from config
    """
    args = type('Args', (), {})()

    args.config = config
    args.tree = Path(tree) if tree else None
    args.sequences = Path(sequences) if sequences else None
    args.prefix = Path(prefix) if prefix else None
    args.table = Path(table) if table else None
    args.segment = segment
    args.command = None
    config = build_config(args)

    # Create output directory if it doesn't exist
    prefix = Path(config["PREFIX"])
    os.makedirs(prefix, exist_ok=True)

    # Extract segment for threshold determination
    if (segment := config.get("segment")) is None:
        raise ValueError("Segment not provided in config or CLI arguments")

    if (species := config.get("SPECIES")) is None:
        raise ValueError("Species not provided in config or CLI arguments")

    # Setup logging
    setup_logging(level="INFO",
                 filepath=str(prefix / f"{species}-{segment}.log"))

    # Log start of the process with parameters
    logger.info("Starting duplicate detection process.")

    if (thresholds := config.get("DEDUPLICATE", {}).get("THRESHOLDS", {}).get(segment)) is None:
        logger.error(f"No thresholds found for segment {segment}")

    thresholds = {k:float(v) for k,v in thresholds.items()}

    logger.info("Using thresholds for segment %s: %s", segment, thresholds)

    start_time = time.time()

    # Load and root the tree using the modular function
    tree_path = config["TREE"]
    tree = root_tree_at_midpoint(tree_path)

    # Calculate distances from tree
    tips, dist_matrix = to_distance_matrix(tree["PHYLODM"])

    # Write distance matrix to file
    write_distance_matrix(dist_matrix, tips, str(prefix / "distance_matrix.mldist"))

    if (length_column := config.get("LENGTH_COLUMN") or config.get("DEDUPLICATE", {}).get("LENGTH_COLUMN")) is None:
        logger.error("LENGTH_COLUMN not found in config or CLI arguments")
    if not (selection_columns := config.get("DEDUPLICATE", {}).get("SELECTION_COLUMNS") or config.get("SELECTION_COLUMNS", [])):
        logger.warning("No selection columns found in config or CLI arguments")

    # Create dictionary of all contigs with their rank {sample: {rank: x, ...}}
    contig_table = config["CONTIGS_TABLE"]
    contigs_ranked = sort_table(contig_table, length_column,
        selection_columns, expected_length=thresholds["TARGET_LENGTH"])

    if (sample_regex := config.get("SAMPLE_REGEX") or config.get("DEDUPLICATE", {}).get("SAMPLE_REGEX")) is None:
        logger.warning("SAMPLE_REGEX not found in config or CLI arguments - defaulting to '(LVE\\d+)'")
        sample_regex = r'(LVE\d+)'

    # Group sequences by sample
    sample_to_seqs = group_sequences_by_sample(tips, sample_regex)

    # Find duplicates with the determined thresholds - pass segment parameter
    classifications = find_duplicates(
        sample_to_seqs, tips, dist_matrix, contigs_ranked,
        tree["BIOPHYLO"], segment, thresholds
        )

    # Write results to output files
    sequence_path = config["SEQUENCES"]
    logger.info("Loading sequences from %s", sequence_path)
    seq_records = SeqIO.to_dict(SeqIO.parse(str(sequence_path), "fasta"))
    write_results(classifications, seq_records, species, segment, prefix, sample_regex)

    # Log completion time
    elapsed = time.time() - start_time
    logger.info("Duplicate detection completed in %.2f seconds", elapsed)

    # Return classifications for testing purposes
    return classifications
