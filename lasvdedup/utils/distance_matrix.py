import numpy as np
import logging
import time
from typing import List, Tuple, Dict
from phylodm import PhyloDM
from Bio import Phylo

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

def get_outliers(clade: Phylo, seq_names: List[str], evolution_threshold: float, z_threshold: float = 2) -> Dict[str, Dict[str, float]]:
    """Identify sequences that are outliers based on their distances."""
    logger.debug("Checking for outliers among %d sequences", len(seq_names))
    outliers = {}
    terminals = clade.get_terminals()
    clade_depths = clade.depths()
    distances = [clade_depths[terminal] for terminal in terminals]

    median = np.median(distances)
    mad = np.median(np.abs(distances - median))
    outlier_threshold = median + z_threshold * mad
    if mad == 0:
        logger.warning("MAD is zero - cannot calculate threshold, using evolutionary threshold %.4f", evolution_threshold)
        # Assume evolutionary threshold +/- 6 months of evolution rate:
        outlier_threshold = median + evolution_threshold

    # Identify outliers
    for seq in seq_names:
        dist = clade_depths[clade.find_any(seq)]
        logger.debug("Distance of %s: %.4f with median: %.4f and threshold of %.4f", seq, dist, median, outlier_threshold)
        if dist > outlier_threshold:
            outliers[seq] = {
                'distance': float(dist),
                'median': float(median),
                'threshold': float(outlier_threshold)
            }

    return outliers
