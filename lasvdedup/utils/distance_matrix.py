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

def calculate_dist_bootstrap(distances: np.ndarray, n_bootstraps: int = 1000, quantile: float = 0.85) -> float:
    """Calculate the bootstrap estimate of variation in distances.

    When MAD is zero, this provides an alternative way to establish a threshold.
    """
    n = len(distances)
    samples = np.random.choice(distances, size=(n_bootstraps, n), replace=True)
    mad_boot = np.apply_along_axis(np.mean, 1, samples)
    return np.quantile(mad_boot, quantile)

def get_outliers( clade_members: List[str], seq_names: List[str], tips_lookup: Dict[str, int], dist_matrix: np.ndarray, quantile: float = 85) -> Dict[str, Dict[str, float]]:
    """Identify sequences that are outliers based on their distances."""
    logger.debug("Checking for outliers among %d sequences", len(seq_names))
    outliers = {}
    references = set(clade_members) - set(seq_names)
    reference = references.pop() # Select one random
    refid = tips_lookup[reference]
    distances = [
        dist_matrix[refid, tips_lookup[member]]
        for member in clade_members
        if member != reference
    ]

    outlier_threshold = calculate_dist_bootstrap(distances, quantile=quantile)
    median = np.median(distances)
    logger.debug("Outlier threshold set to %.4f", outlier_threshold)

    # Identify outliers
    for seq in seq_names:
        dist = dist_matrix[tips_lookup[reference], tips_lookup[seq]]
        logger.debug("Distance from %s to %s: %.4f with median: %.4f and threshold of %.4f", reference, seq, dist, median, outlier_threshold)
        if dist > outlier_threshold:
            outliers[seq] = {
                'distance': float(dist),
                'median': float(median),
                'threshold': float(outlier_threshold),
                'reference': reference
            }

    return outliers
