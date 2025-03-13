import numpy as np
import logging
import time
from typing import List, Tuple, Dict
from phylodm import PhyloDM

logger = logging.getLogger("lasvdedup.duplicates")

def to_distance_matrix(tree: PhyloDM) -> Tuple[List, np.ndarray]:
    """Create a distance matrix from a phylogenetic tree."""
    logger.info("Calculating distance matrix from phylogenetic tree")
    start_time = time.time()

    dm = tree.dm(norm=False)
    labels = [l.replace(" ", "_") for l in tree.taxa()]

    elapsed = time.time() - start_time
    logger.info("Distance matrix calculation completed in %.2f seconds", elapsed)
    return (labels, dm)

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

def get_outliers(clade_members: List[str], seq_names: List[str],
        tips_lookup: Dict[str, int], dist_matrix: np.ndarray, z_threshold: float = 2.0) -> Dict[str, Dict[str, float]]:
    """Identify sequences that are outliers based on their distances."""
    logger.debug("Checking for outliers among %d sequences", len(seq_names))
    outliers = {}
    mad = 0
    references = set(clade_members) - set(seq_names)
    distances = []

    while mad == 0:
        logger.debug(f"Currently the mad is {mad}, trying to find a reference sequence")
        # Select a random reference sequence from those in seq_names
        if len(references) == 0:
            logger.debug("No reference sequences left to try")
            break

        reference = references.pop()
        logger.debug("Using %s as reference sequence", reference)
        logger.debug("References: %s", references)
        references = references - {reference}

        # Collect distances from reference to other clade members
        refid = tips_lookup[reference]
        distances = [
            dist_matrix[refid, tips_lookup[member]]
            for member in clade_members
            if member != reference
        ]

        for member in seq_names:
            if member == reference:
                continue
            logger.debug(f"Distance from {member} to ref {reference}: {dist_matrix[refid, tips_lookup[member]]}")

        # Calculate statistics
        median = np.median(distances)
        logger.debug("Median distance from %s to clade members: %.4f", reference, median)
        logger.debug("distances to median: %s", np.sort(np.abs(distances - median)))
        mad = np.median(np.abs(distances - median))

    if mad == 0:
        logger.warning("Could not calculate MAD for clade %s, determining variance", clade_members)
        mad = np.std(distances)

    logger.debug("MAD: %.4f for clade %s", mad, clade_members)

    outlier_threshold = z_threshold * mad

    # Identify outliers
    for seq in seq_names:
        dist = dist_matrix[tips_lookup[reference], tips_lookup[seq]]
        logger.debug("Distance from %s to %s: %.4f with median: %.4f and threshold of %.4f", reference, seq, dist, median, outlier_threshold)
        if np.abs(dist - median) > outlier_threshold:
            outliers[seq] = {
                'distance': float(dist),
                'median': float(median),
                'threshold': float(outlier_threshold),
                'reference': reference
            }

    return outliers
