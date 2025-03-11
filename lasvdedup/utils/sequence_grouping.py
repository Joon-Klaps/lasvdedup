import re
import logging
from typing import Dict, List
import numpy as np
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

def select_best_sequence(seq_names, reads):
    """Select the sequence with the highest read count."""
    return max(reads, key=reads.get)

def get_reads (seq_name: str, contig_to_reads: Dict[str, Dict[str, int]]) -> int:
    """
    Get the number of reads for a given sequence.

    Args:
        seq_name: Sequence name
        contig_to_reads: Dictionary mapping contig IDs to read counts
    """
    if seq_name in contig_to_reads:
        return contig_to_reads[seq_name]["reads"]
    elif seq_name.replace('_R_', '').split(".", 1)[0] in contig_to_reads:
        return contig_to_reads[seq_name.replace('_R_', '').split(".", 1)[0]]["reads"]

    logger.warning(f"Read count not found for sequence: {seq_name}")
    raise ValueError(f"Read count not found for sequence: {seq_name}")


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
    contig_to_reads: Dict[str, Dict[str, int]],
    lowerthreshold: float,
    upperthreshold: float,
    tree=None,
    clade_size: int = 10,
    z_threshold: float = 2.0
) -> Dict[str, Classification]:
    """
    Process each sample group to identify duplicates based on phylogenetic distances.

    Args:
        sample_to_seqs: Dictionary mapping sample IDs to lists of sequence names
        tips: List of tree tips (sequence names)
        dist_matrix: Distance matrix as numpy array
        contig_to_reads: Dictionary mapping contig IDs to read counts
        lowerthreshold: Minimum distance to identify duplicates
        upperthreshold: Maximum distance for intrahost variation
        tree: Bio.phylo tree object for MRCA analysis
        clade_size: Minimum clade size to consider a coinfection
        z_threshold: Z-score threshold for outlier detection

    Returns:
        Dictionary mapping each tip name to its Classification object
    """
    logger.info("Finding duplicates with thresholds: lower=%.4f, upper=%.4f, clade_size=%d, z_threshold=%.2f",
               lowerthreshold, upperthreshold, clade_size, z_threshold)

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

        # Classify sequences in this sample
        sample_classifications = classify_sample(
            sample_id, seq_names, tips_lookup, dist_matrix,
            contig_to_reads, lowerthreshold, upperthreshold, tree, clade_size, z_threshold
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
    sample_id,
    seq_names,
    tips_lookup,
    dist_matrix,
    contig_to_reads,
    lowerthreshold,
    upperthreshold,
    tree=None,
    min_clade_size=10,
    z_threshold=2.0
):
    """
    Classify sequences from a single sample based on phylogenetic distances.

    Args:
        sample_id: Sample identifier
        seq_names: List of sequence names from this sample
        tips_lookup: Dictionary mapping sequence names to matrix indices
        dist_matrix: Distance matrix as numpy array
        contig_to_reads: Dictionary mapping sequence names to read counts
        lowerthreshold: Lower distance threshold for duplicates
        upperthreshold: Upper distance threshold for intrahost variants
        tree: Bio.Phylo tree object for MRCA analysis
        min_clade_size: Minimum clade size to consider a coinfection
        z_threshold: Z-score threshold for outlier detection

    Returns:
        Dictionary mapping sequence names to Classification objects
    """
    # Initialize results dictionary
    classifications = {}
    group_members = list(seq_names)  # All sequences in this sample

    # Handle single sequence case
    if len(seq_names) <= 1:
        for seq in seq_names:
            classifications[seq] = Classification(
                sequence_name=seq,
                classification_type=ClassificationType.GOOD,
                reason="Single sequence in sample",
                sample_id=sample_id,
                group_members=group_members,
                decision_category=DecisionCategory.SINGLE_SEQUENCE,  # Set decision category
                read_count=get_reads(seq, contig_to_reads) if contig_to_reads else None
            )
        return classifications

    # Get pairwise distances and read counts
    distances = get_distances(seq_names, tips_lookup, dist_matrix)
    reads = {seq: get_reads(seq, contig_to_reads) for seq in seq_names}

    # Case 1: All distances below lower threshold (simple duplicates)
    if is_all_distances_below_threshold(distances, lowerthreshold):
        logger.debug("Sample %s: All sequences are similar (below lower threshold)", sample_id)
        best_seq = select_best_sequence(seq_names, reads)

        classifications[best_seq] = Classification(
            sequence_name=best_seq,
            classification_type=ClassificationType.GOOD,
            reason=f"Selected as best representative (highest read count: {reads[best_seq]}) from nearly identical sequences (distances < {lowerthreshold})",
            sample_id=sample_id,
            group_members=group_members,
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,  # Set decision category
            read_count=reads[best_seq]
        )

        for seq in set(seq_names) - {best_seq}:
            classifications[seq] = Classification(
                sequence_name=seq,
                classification_type=ClassificationType.BAD,
                reason=f"Duplicate of {best_seq} (distance < {lowerthreshold}, lower read count ({reads[seq]}) than {best_seq} ({reads[best_seq]}))",
                sample_id=sample_id,
                group_members=group_members,
                decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,  # Set decision category
                read_count=reads[seq]
            )

        return classifications

    # Case 2: All distances below upper threshold (potential intrahost variants)
    if is_all_distances_below_threshold(distances, upperthreshold):
        logger.debug("Sample %s: Potential intrahost variants (below upper threshold)", sample_id)

        # Group sequences into clusters based on lower threshold
        clusters = cluster_sequences(seq_names, tips_lookup, dist_matrix, lowerthreshold)

        # For each cluster, keep the sequence with highest read count
        for cluster in clusters:
            cluster_reads = {seq: reads[seq] for seq in cluster}
            best_seq = select_best_sequence(cluster, cluster_reads)

            classifications[best_seq] = Classification(
                sequence_name=best_seq,
                classification_type=ClassificationType.COINFECTION,
                reason=f"Intrahost variant (distances < {upperthreshold}), selected as cluster representative (highest read count: {reads[best_seq]})",
                sample_id=sample_id,
                group_members=cluster,
                decision_category=DecisionCategory.BELOW_UPPER_THRESHOLD,  # Set decision category
                read_count=reads[best_seq]
            )

            for seq in set(cluster) - {best_seq}:
                classifications[seq] = Classification(
                    sequence_name=seq,
                    classification_type=ClassificationType.BAD,
                    reason=f"Intrahost variant duplicate of {best_seq} (distance < {lowerthreshold}, lower read count ({reads[seq]}) than {[best_seq]}",
                    sample_id=sample_id,
                    group_members=cluster,
                    decision_category=DecisionCategory.BELOW_UPPER_THRESHOLD,  # Set decision category
                    read_count=reads[seq]
                )

        return classifications

    # Case 3: Some distances above upper threshold - check MRCA clade size
    logger.debug("Sample %s: Some distances above upper threshold - checking MRCA", sample_id)

    # Skip MRCA analysis if tree is not provided
    if tree is None:
        logger.warning("Sample %s: No tree provided for MRCA analysis, treating as coinfection", sample_id)
        for seq in seq_names:
            classifications[seq] = Classification(
                sequence_name=seq,
                classification_type=ClassificationType.COINFECTION,
                reason="Tree not provided for MRCA analysis, defaulting to coinfection",
                sample_id=sample_id,
                group_members=group_members,
                decision_category=DecisionCategory.NO_TREE,  # Set decision category
                read_count=reads[seq]
            )
        return classifications

    # Get MRCA clade members
    clade_members = get_mrca_clade(seq_names, tree)
    clade_size = len(clade_members)

    # Case 4: Small MRCA clade (likely false positive)
    if clade_size <= min_clade_size:
        logger.debug("Sample %s: Small MRCA clade size (%d) - likely false positive", sample_id, clade_size)

        best_seq = select_best_sequence(seq_names, reads)
        classifications[best_seq] = Classification(
            sequence_name=best_seq,
            classification_type=ClassificationType.GOOD,
            reason=f"Small MRCA clade size ({clade_size} ≤ {min_clade_size}) indicating likely false positive, selected as representative (highest read count: {reads[best_seq]})",
            sample_id=sample_id,
            group_members=group_members,
            decision_category=DecisionCategory.SMALL_CLADE,  # Set decision category
            read_count=reads[best_seq]
        )

        for seq in set(seq_names) - {best_seq}:
            classifications[seq] = Classification(
                sequence_name=seq,
                classification_type=ClassificationType.BAD,
                reason=f"Likely false positive with small MRCA clade size ({clade_size} ≤ {min_clade_size}), {best_seq} selected instead (higher read count)",
                sample_id=sample_id,
                group_members=group_members,
                decision_category=DecisionCategory.SMALL_CLADE,  # Set decision category
                read_count=reads[seq]
            )

        return classifications

    # Case 5: Check for outliers - use the provided z_threshold parameter
    logger.debug("Sample %s: Large MRCA clade size (%d) - checking for outliers", sample_id, clade_size)

    outliers = get_outliers(clade_members, seq_names, tips_lookup, dist_matrix, z_threshold)

    if outliers:
        logger.debug("Sample %s: Outliers found - likely false positive", sample_id)
        outlier_names = ", ".join(outliers)

        good_seqs = [seq for seq in seq_names if seq not in outliers]
        best_seq = select_best_sequence(good_seqs, reads)

        classifications[best_seq] = Classification(
            sequence_name=best_seq,
            classification_type=ClassificationType.GOOD,
            reason=f"Outliers detected ({outlier_names}), selected as best non-outlier sequence (highest read count: {reads[best_seq]})",
            sample_id=sample_id,
            group_members=group_members,
            decision_category=DecisionCategory.OUTLIERS_DETECTED,  # Set decision category
            read_count=reads[best_seq]
        )

        for seq in set(seq_names) - {best_seq}:
            if seq in outliers:
                classifications[seq] = Classification(
                    sequence_name=seq,
                    classification_type=ClassificationType.BAD,
                    reason=f"Identified as phylogenetic outlier, {best_seq} selected instead (non-outlier with highest read count)",
                    sample_id=sample_id,
                    group_members=group_members,
                    decision_category=DecisionCategory.OUTLIERS_DETECTED,  # Set decision category
                    read_count=reads[seq]
                )
            else:
                classifications[seq] = Classification(
                    sequence_name=seq,
                    classification_type=ClassificationType.BAD,
                    reason=f"Non-outlier but with lower read count ({reads[seq]}) than {best_seq} ({reads[best_seq]})",
                    sample_id=sample_id,
                    group_members=group_members,
                    decision_category=DecisionCategory.OUTLIERS_DETECTED,  # Set decision category
                    read_count=reads[seq]
                )

        return classifications

    # Case 6: True coinfection
    logger.debug("Sample %s: No outliers in large clade (%d) - TRUE COINFECTION", sample_id, clade_size)

    for seq in seq_names:
        classifications[seq] = Classification(
            sequence_name=seq,
            classification_type=ClassificationType.COINFECTION,
            reason=f"True coinfection: large MRCA clade size ({clade_size} > {min_clade_size}) and no outliers detected",
            sample_id=sample_id,
            group_members=group_members,
            decision_category=DecisionCategory.TRUE_COINFECTION,  # Set decision category
            read_count=reads[seq]
        )

    return classifications
