#!/usr/bin/env python3

import os
import logging
import time
from Bio import SeqIO
import yaml

# Import sub utility functions
from .tree_utils import root_tree_at_midpoint
from .distance_matrix import to_distance_matrix
from .io_utils import load_read_counts, write_results, write_distance_matrix
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

def get_segment_thresholds(config, segment):
    """
    Get the appropriate thresholds for the given segment.

    Args:
        config: Configuration dictionary or segment name
        segment: Segment name (e.g. 'L', 'S')

    Returns:
        tuple: (lower_threshold, upper_threshold, clade_size, z_threshold)
    """
    # If config is a dictionary, extract thresholds from it
    logger.debug("Getting segment thresholds for %s", segment)
    logger.debug("Config: %s", config)

    if isinstance(config, dict):
        thresholds = config.get('DEDUPLICATE', {}).get('THRESHOLDS', {})
        segment_thresholds = thresholds.get(segment, None)

        if segment_thresholds:
            return (
                float(segment_thresholds.get('LOWER')),
                float(segment_thresholds.get('UPPER')),
                int(segment_thresholds.get('CLADE_SIZE', 10)),
                float(segment_thresholds.get('Z_THRESHOLD', 2.0))
            )

        # Use default thresholds if segment-specific ones aren't defined
        default = config.get('DEDUPLICATE', {}).get('DEFAULT_THRESHOLD', {})
        return (
            float(default.get('LOWER', config.get('lowerthreshold', 0.02))),
            float(default.get('UPPER', config.get('upperthreshold', 0.05))),
            int(default.get('CLADE_SIZE', config.get('clade_size', 10))),
            float(default.get('Z_THRESHOLD', config.get('z_threshold', 2.0)))
        )

    # If no config dictionary was provided, use the parameters directly
    return (
        float(config.get('lowerthreshold', 0.02)),
        float(config.get('upperthreshold', 0.05)),
        int(config.get('clade_size', 10)),
        float(config.get('z_threshold', 2.0))
    )

def determine_duplicates(config=None, **kwargs):
    """
    Main function to detect duplicate sequences in a phylogenetic tree.

    Args:
        config: Path to config file or config dictionary (optional)
        **kwargs: Additional keyword arguments including:
            tree: Path to the phylogenetic tree file
            sequences: Path to the sequences file (FASTA format)
            prefix: Output file prefix
            table: Path to the table with contig information
            sample_regex: Regular expression to extract sample identifiers
            reads_column: Name of the column containing read counts
            species: Species name for the output files
            segment: Segment name for the output files
            lowerthreshold: Lower distance threshold (if not using config)
            upperthreshold: Upper distance threshold (if not using config)
            log_level: Logging level
    """
    # Load configuration if path is provided
    if config and isinstance(config, str):
        try:
            with open(config, 'r') as f:
                config_data = yaml.safe_load(f)
        except Exception as e:
            logger.warning(f"Failed to load config from {config}: {e}")
            config_data = {}
    elif config and isinstance(config, dict):
        config_data = config
    else:
        config_data = {}

    # Create output directory if it doesn't exist
    prefix = kwargs["prefix"]
    prefix_path = os.path.abspath(prefix)
    os.makedirs(prefix, exist_ok=True)

    # Extract segment for threshold determination
    segment = kwargs.get("segment", None)
    species = kwargs.get("species") or config_data.get('SPECIES', 'LASV')

    # Setup logging
    log_level = kwargs.get("log_level", "INFO")
    setup_logging(log_level, f"{prefix_path}/{species}-{segment}.log")

    # Log start of the process with parameters
    logger.info("Starting duplicate detection process")

    # Extract segment for threshold determination
    segment = kwargs.get("segment", None)
    # Get ALL thresholds based on segment
    lowerthreshold, upperthreshold, clade_size, z_threshold = get_segment_thresholds(
        config_data if config_data else kwargs,
        segment
    )

    logger.info("Using thresholds for segment %s: lower=%f, upper=%f, clade_size=%d, z_threshold=%f",
                segment, lowerthreshold, upperthreshold, clade_size, z_threshold)

    start_time = time.time()

    # Load and root the tree using the modular function
    tree = root_tree_at_midpoint(kwargs["tree"])

    # Calculate distances from tree
    tips, dist_matrix = to_distance_matrix(tree["phylodm"])

    # Write distance matrix to file
    write_distance_matrix(dist_matrix, tips, f"{prefix}/distance_matrix.mldist")

    # Load read counts from contigs table
    reads_column = kwargs.get("reads_column") or config_data.get('DEDUPLICATE', {}).get('READS_COLUMN', 'reads')
    contig_to_reads = load_read_counts(kwargs["table"], reads_column)

    # Get sample regex from config or kwargs with a default pattern
    sample_regex = kwargs.get("sample_regex") or config_data.get('DEDUPLICATE', {}).get('SAMPLE_REGEX', r'(LVE\d+)_.*')

    # Group sequences by sample
    sample_to_seqs = group_sequences_by_sample(tips, sample_regex)

    # Find duplicates with the determined thresholds
    classifications = find_duplicates(
        sample_to_seqs, tips, dist_matrix, contig_to_reads,
        lowerthreshold, upperthreshold, tree["biophylo"],
        clade_size=clade_size, z_threshold=z_threshold
    )

    # Write results to output files
    logger.info("Loading sequences from %s", kwargs['sequences'])
    seq_records = SeqIO.to_dict(SeqIO.parse(str(kwargs["sequences"]), "fasta"))
    write_results(classifications, seq_records, species, segment, prefix_path, sample_regex)

    # Log completion time
    elapsed = time.time() - start_time
    logger.info("Duplicate detection completed in %.2f seconds", elapsed)

    # Return classifications for testing purposes
    return classifications
