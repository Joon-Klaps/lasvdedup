import re
import logging
from typing import Dict, List
import numpy as np
import Bio.Phylo
from .classification import Classification, ClassificationType, DecisionCategory

from .distance_matrix import get_distances, get_outliers
from .tree_utils import get_mrca_clade

logger = logging.getLogger("lasvdedup.duplicates")

def group_sequences_by_sample(tips: List, sample_regex: str) -> Dict[str, List[str]]:
    """
    Group sequences by sample ID using regex.

    Args:
        tips: List of tree tips
        sample_regex: Regular expression to extract sample identifiers

    Returns:
        Dictionary mapping sample IDs to lists of sequence names
    """
    logger.info("Grouping sequences by sample using regex: %s", sample_regex)
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
        logger.warning("No sample ID match for %d sequences", no_match_count)

    logger.info("Found %d unique sample IDs", len(sample_to_seqs))
    # Log stats on sequences per sample
    seqs_per_sample = [len(seqs) for seqs in sample_to_seqs.values()]
    if seqs_per_sample:
        logger.info("Sequences per sample - Min: %d, Max: %d, Avg: %.2f",
                   min(seqs_per_sample), max(seqs_per_sample),
                   sum(seqs_per_sample)/len(seqs_per_sample))
    return sample_to_seqs

def is_all_distances_below_threshold(distances, threshold):
    """Check if all pairwise distances are below a given threshold."""
    return all(d <= threshold for d in distances)

def get_contig_data(seq_name: str, contigs_ranked: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    """
    Get the data (reads, coverage, qlen) for a given sequence.

    Args:
        seq_name: Sequence name
        contigs_ranked: Dictionary mapping contig IDs to data values
    """
    # Try exact match first
    if seq_name in contigs_ranked:
        return contigs_ranked[seq_name]

    # Try with _R_ replacement and splitting at first dot
    base_name = seq_name.replace('_R_', '').split(".", 1)[0]
    if base_name in contigs_ranked:
        return contigs_ranked[base_name]

    # Try with just the part before any dots (for cases like seq.with.dots.extra)
    for key in contigs_ranked.keys():
        # Check if the key is a prefix of seq_name or matches before any dots
        if seq_name.startswith(key) and (len(seq_name) == len(key) or seq_name[len(key)] == '.'):
            return contigs_ranked[key]

    logger.warning(f"Data not found for sequence: {seq_name}")
    raise ValueError(f"Data not found for sequence: {seq_name}")

def select_best_sequence(seq_names, seq_data):
    """
    Select the best sequence based on previously identified rank
    """
    if not seq_names:
        raise ValueError("Empty sequence list provided")

    return min(seq_names, key=lambda seq: seq_data[seq]["rank"])

def cluster_sequences(seq_names, tips_lookup, dist_matrix, threshold):
    """
    Group sequences into clusters based on a distance threshold.

    Returns:
        List of lists, where each inner list represents a cluster of sequence names
    """
    clusters = []
    unclustered = set(seq_names)

    while unclustered:
        seed = next(iter(unclustered))
        cluster = {seed}
        unclustered.remove(seed)

        while True:
            # Find all sequences that are within threshold distance of any sequence in the cluster
            to_add = {seq for seq in unclustered
                     if any(dist_matrix[tips_lookup[seq], tips_lookup[clustered_seq]] <= threshold
                           for clustered_seq in cluster)}

            if not to_add:
                break

            cluster.update(to_add)
            unclustered -= to_add

        clusters.append(list(cluster))

    return clusters

def find_duplicates(
    sample_to_seqs: Dict[str, List[str]],
    tips: List,
    dist_matrix: np.ndarray,
    contigs_ranked: Dict[str, Dict[str, float]],
    tree: Bio.Phylo = None,
    segment: str = None,
    thresholds: Dict[str, float] = None
) -> Dict[str, Classification]:
    """
    Process each sample group to identify duplicates based on phylogenetic distances.

    Args:
        sample_to_seqs: Dictionary mapping sample IDs to lists of sequence names
        tips: List of tree tips (sequence names)
        dist_matrix: Distance matrix as numpy array
        contigs_ranked: Dictionary mapping contig IDs to ranks
        tree: Bio.Phylo tree object
        segment: Segment name (L or S) to determine target length
        thresholds: Dictionary of thresholds for PWD, Z-score, clade size

    Returns:
        Dictionary mapping each tip name to its Classification object
    """
    logger.info("Finding duplicates with thresholds: %s", thresholds)

    # Create lookup for efficient index retrieval
    tips_lookup = {t: i for i, t in enumerate(tips)}

    # Initialize results dictionary
    all_classifications = {}

    # Track statistics
    good_count = 0
    bad_count = 0
    coinfection_count = 0

    # Process each sample
    total_samples = len(sample_to_seqs)
    for i, (sample_id, seq_names) in enumerate(sample_to_seqs.items()):
        logger.debug("Processing sample %d/%d: %s with %d sequences",
                i+1, total_samples, sample_id, len(seq_names))

        # Classify sequences in this sample - pass segment parameter
        sample_classifications = classify_sample(
            sample_id, seq_names, tips_lookup, dist_matrix,
            contigs_ranked, tree, segment, thresholds
        )

        # Update statistics
        for seq, classification in sample_classifications.items():
            if classification.is_good:
                good_count += 1
            elif classification.is_bad:
                bad_count += 1
            elif classification.is_coinfection:
                coinfection_count += 1

        # Merge classifications
        all_classifications.update(sample_classifications)

    logger.info("Classification results: %d good, %d bad, %d coinfection",
               good_count, bad_count, coinfection_count)
    return all_classifications

def classify_sample(
    sample_id: str,
    seq_names: List[str],
    tips_lookup: Dict[str, int],
    dist_matrix: np.ndarray,
    contigs_ranked: Dict[str, Dict[str, float]],
    tree: Bio.Phylo = None,
    segment: str = None,
    thresholds: Dict[str, float] = None
):
    """
    Classify sequences from a single sample based on phylogenetic distances.

    Args:


    Returns:
        Dictionary mapping sequence names to Classification objects
    """
    # Initialize results dictionary
    classifications = {}
    group_members = list(seq_names)  # All sequences in this sample

    # Case 1:  Handle single sequence case
    if len(seq_names) <= 1:
        for seq in seq_names:
            classifications[seq] = Classification(
                sequence_name=seq,
                classification_type=ClassificationType.GOOD,
                reason="Single sequence in sample",
                sample_id=sample_id,
                group_members=group_members,
                decision_category=DecisionCategory.SINGLE_SEQUENCE,  # Set decision category
                contig_stats=contigs_ranked[seq]
            )
        return classifications

    # Get pairwise distances and ranks
    distances = get_distances(seq_names, tips_lookup, dist_matrix)
    stats = {seq: get_contig_data(seq, contigs_ranked) for seq in seq_names}

    # Case 2: All distances below lower threshold (simple duplicates)
    if is_all_distances_below_threshold(distances, thresholds["PWD"]):
        logger.debug("Sample %s: All sequences are similar (below lower threshold)", sample_id)
        best_seq = select_best_sequence(seq_names, stats)

        classifications[best_seq] = Classification(
            sequence_name=best_seq,
            classification_type=ClassificationType.GOOD,
            reason=f"Selected as best representative (contig stats: {stats[best_seq]}) from nearly identical sequences (distances < {thresholds['PWD']})",
            sample_id=sample_id,
            group_members=group_members,
            decision_category=DecisionCategory.BELOW_THRESHOLD,  # Set decision category
            contig_stats=stats[best_seq]
        )

        for seq in set(seq_names) - {best_seq}:
            classifications[seq] = Classification(
                sequence_name=seq,
                classification_type=ClassificationType.BAD,
                reason=f"Duplicate of {best_seq} (distance < {thresholds['PWD']}, lower contig stats ({stats[seq]}) than {stats[best_seq]})",
                sample_id=sample_id,
                group_members=group_members,
                decision_category=DecisionCategory.BELOW_THRESHOLD,  # Set decision category
                contig_stats=stats[seq]
            )

        return classifications

    logger.debug("Sample %s: Some distances above upper threshold - checking MRCA", sample_id)

    # Get MRCA clade members
    clade_members = get_mrca_clade(seq_names, tree)
    clade_size = len(clade_members)

    # Case 3: Small MRCA clade (likely false positive)
    if clade_size <= thresholds["CLADE_SIZE"]:
        logger.debug("Sample %s: Small MRCA clade size (%d) - likely false positive", sample_id, clade_size)

        best_seq = select_best_sequence(seq_names, contigs_ranked)
        classifications[best_seq] = Classification(
            sequence_name=best_seq,
            classification_type=ClassificationType.GOOD,
            reason=f"Small MRCA clade size ({clade_size} ≤ {thresholds['CLADE_SIZE']}) indicating likely false positive, selected as representative (highest rank: {stats[best_seq]})",
            sample_id=sample_id,
            group_members=group_members,
            decision_category=DecisionCategory.SMALL_CLADE,  # Set decision category
            contig_stats=stats[best_seq]
        )

        for seq in set(seq_names) - {best_seq}:
            classifications[seq] = Classification(
                sequence_name=seq,
                classification_type=ClassificationType.BAD,
                reason=f"Likely false positive with small MRCA clade size ({clade_size} ≤ {thresholds['CLADE_SIZE']}), {best_seq} selected instead (higher rank)",
                sample_id=sample_id,
                group_members=group_members,
                decision_category=DecisionCategory.SMALL_CLADE,  # Set decision category
                contig_stats=stats[seq]
            )

        return classifications

    # Case 4: Check for outliers - use the provided z_threshold parameter
    logger.debug("Sample %s: Large MRCA clade size (%d) - checking for outliers", sample_id, clade_size)

    outliers = get_outliers(clade_members, seq_names, tips_lookup, dist_matrix, thresholds["Z_THRESHOLD"])

    if outliers:
        outlier_names = ", ".join(outliers.keys())
        logger.debug("Sample %s: Outliers found %s - likely false positive", sample_id, outlier_names)

        good_seqs = [seq for seq in seq_names if seq not in outliers]
        if not good_seqs:
            logger.error("No good sequences found after removing outliers: %s", outlier_names)

        best_seq = select_best_sequence(good_seqs, contigs_ranked)

        classifications[best_seq] = Classification(
            sequence_name=best_seq,
            classification_type=ClassificationType.GOOD,
            reason=f"Outliers detected ({outlier_names}), selected as best non-outlier sequence (highest rank: {stats[best_seq]})",
            sample_id=sample_id,
            group_members=group_members,
            decision_category=DecisionCategory.OUTLIERS_DETECTED,
            contig_stats=stats[best_seq]
        )

        for seq in set(seq_names) - {best_seq}:
            if seq in outliers:
                outlier_info = outliers[seq]
                classifications[seq] = Classification(
                    sequence_name=seq,
                    classification_type=ClassificationType.BAD,
                    reason=f"Identified as phylogenetic outlier (distance - median: {outlier_info['distance'] - outlier_info['median']:.4f}, threshold: {outlier_info['threshold']:.4f}), {best_seq} selected instead (non-outlier with highest rank)",
                    sample_id=sample_id,
                    group_members=group_members,
                    decision_category=DecisionCategory.OUTLIERS_DETECTED,
                    contig_stats=stats[seq]
                )
            else:
                classifications[seq] = Classification(
                    sequence_name=seq,
                    classification_type=ClassificationType.BAD,
                    reason=f"Non-outlier but with lower rank ({stats[seq]}) than {best_seq} ({stats[best_seq]})",
                    sample_id=sample_id,
                    group_members=group_members,
                    decision_category=DecisionCategory.OUTLIERS_DETECTED,
                    contig_stats=stats[seq]
                )

        return classifications

    # Case 5: True coinfection of 2 distinct strains
    logger.debug("Sample %s: No outliers in large clade (%d) - TRUE COINFECTION", sample_id, clade_size)

    for seq in seq_names:
        classifications[seq] = Classification(
            sequence_name=seq,
            classification_type=ClassificationType.COINFECTION,
            reason=f"True coinfection: large MRCA clade size ({clade_size} > {thresholds['CLADE_SIZE']}) and no outliers detected",
            sample_id=sample_id,
            group_members=group_members,
            decision_category=DecisionCategory.TRUE_COINFECTION,  # Set decision category
            contig_stats=stats[seq]
        )

    return classifications
