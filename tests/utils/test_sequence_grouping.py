import pytest
import numpy as np
from unittest.mock import MagicMock, patch
import logging

from lasvdedup.utils.sequence_grouping import (
    group_sequences_by_sample,
    is_all_distances_below_threshold,
    select_best_sequence,
    get_reads,
    cluster_sequences,
    find_duplicates,
    classify_sample
)
from lasvdedup.utils.classification import Classification, ClassificationType, DecisionCategory

# Set up logging for tests
logging.basicConfig(level=logging.DEBUG)

class TestGroupSequencesBySample:
    """Tests for group_sequences_by_sample function."""

    def test_empty_input(self):
        """Test with empty list of tips."""
        result = group_sequences_by_sample([], r'sample\d+')
        assert result == {}

    def test_basic_grouping(self):
        """Test basic grouping functionality."""
        # Test tips with various patterns
        tips = [
            'sampleA_1_contig',
            'sampleA_2_contig',
            'sampleB_1_contig',
            'noMatch_contig'
        ]

        # Test with pattern matching sampleA_1, sampleA_2, etc.
        sample_to_seqs = group_sequences_by_sample(tips, r'sample[A-Z]_\d')
        assert len(sample_to_seqs) == 3
        assert 'sampleA_1' in sample_to_seqs
        assert 'sampleA_2' in sample_to_seqs
        assert 'sampleB_1' in sample_to_seqs
        assert len(sample_to_seqs['sampleA_1']) == 1

        # Test with more generic pattern
        sample_to_seqs = group_sequences_by_sample(tips, r'sample[A-Z]')
        assert len(sample_to_seqs) == 2  # sampleA and sampleB
        assert len(sample_to_seqs['sampleA']) == 2  # Two sampleA sequences

    def test_complex_regex(self):
        """Test with more complex regex patterns."""
        tips = [
            'Patient123_Sample456_Run789',
            'Patient123_Sample457_Run790',
            'Patient124_Sample458_Run791',
            'NoMatch_Run792'
        ]

        # Match full pattern with capturing groups
        result = group_sequences_by_sample(tips, r'Patient(\d+)_Sample(\d+)')
        assert len(result) == 3
        assert 'Patient123_Sample456' in result
        assert 'Patient123_Sample457' in result
        assert 'Patient124_Sample458' in result

        # Match using just the patient ID
        result = group_sequences_by_sample(tips, r'Patient\d+')
        assert len(result) == 2
        assert 'Patient123' in result
        assert len(result['Patient123']) == 2

    def test_overlapping_patterns(self):
        """Test with patterns that could match multiple parts of the string."""
        tips = [
            'Sample1_Sample2_Sample3',
            'Sample2_OtherData',
            'PrefixSample1Suffix'
        ]

        # Should match the first occurrence in each string
        result = group_sequences_by_sample(tips, r'Sample\d+')
        assert len(result) == 2
        assert 'Sample1' in result
        assert 'Sample2' in result
        assert len(result['Sample1']) == 2  # Matches in first and third tips


class TestDistanceThreshold:
    """Tests for is_all_distances_below_threshold function."""

    def test_basic_threshold(self):
        """Test basic threshold comparisons."""
        # All below
        distances = [0.01, 0.02, 0.03]
        assert is_all_distances_below_threshold(distances, 0.05)

        # Some above
        distances = [0.01, 0.06, 0.03]
        assert not is_all_distances_below_threshold(distances, 0.05)

        # Edge case - exactly at threshold
        distances = [0.01, 0.05, 0.03]
        assert is_all_distances_below_threshold(distances, 0.05)

        # Empty list
        assert is_all_distances_below_threshold([], 0.05)

    def test_negative_distances(self):
        """Test with negative distance values (which shouldn't occur in real data but test robustness)."""
        distances = [-0.01, -0.02, -0.03]
        assert is_all_distances_below_threshold(distances, 0.05)
        assert not is_all_distances_below_threshold(distances, -0.04)

    def test_mixed_distances(self):
        """Test with mix of positive and negative distances."""
        distances = [-0.01, 0.02, -0.03, 0.04]
        assert is_all_distances_below_threshold(distances, 0.05)
        assert not is_all_distances_below_threshold(distances, 0.03)


class TestSelectBestSequence:
    """Tests for select_best_sequence function."""

    def test_basic_selection(self):
        """Test basic selection of best sequence."""
        seq_names = ['seqA', 'seqB', 'seqC']
        reads = {'seqA': 100, 'seqB': 200, 'seqC': 50}

        best_seq = select_best_sequence(seq_names, reads)
        assert best_seq == 'seqB'  # seqB has the highest count (200)

        # Test with tie (should take the first encountered)
        reads = {'seqA': 200, 'seqB': 200, 'seqC': 50}
        best_seq = select_best_sequence(seq_names, reads)
        assert best_seq in ['seqA', 'seqB']  # Either could be selected in case of tie

    def test_with_string_reads(self):
        """Test with string values instead of integers (robustness test)."""
        seq_names = ['seqA', 'seqB', 'seqC']
        reads = {'seqA': '100', 'seqB': '200', 'seqC': '50'}

        with pytest.raises(TypeError):
            select_best_sequence(seq_names, reads)

    def test_with_empty_sequence_list(self):
        """Test behavior with empty sequence list."""
        with pytest.raises(ValueError):
            select_best_sequence([], {'seqA': 100})

    def test_with_missing_read_counts(self):
        """Test behavior when some sequences don't have read counts."""
        seq_names = ['seqA', 'seqB', 'seqC']
        reads = {'seqA': 100, 'seqC': 50}  # seqB missing

        with pytest.raises(KeyError):
            select_best_sequence(seq_names, reads)

    def test_with_infinity_reads(self):
        """Test with infinite read count (edge case)."""
        seq_names = ['seqA', 'seqB', 'seqC']
        reads = {'seqA': 100, 'seqB': float('inf'), 'seqC': 50}

        best_seq = select_best_sequence(seq_names, reads)
        assert best_seq == 'seqB'


class TestGetReads:
    """Tests for get_reads function."""

    def test_basic_reads(self):
        """Test basic read count retrieval."""
        contig_to_reads = {
            'seq1': {'reads': 100},
            'seq2': {'reads': 200},
            'seq3.suffix': {'reads': 300}
        }

        # Direct match
        assert get_reads('seq1', contig_to_reads) == 100

        # Match after _R_ replacement
        assert get_reads('_R_seq2', contig_to_reads) == 200

        # Match with splitting
        assert get_reads('seq3.suffix.moretext', contig_to_reads) == 300

        # Not found
        with pytest.raises(ValueError):
            get_reads('notfound', contig_to_reads)

    def test_with_complex_name_handling(self):
        """Test more complex sequence name handling scenarios."""
        contig_to_reads = {
            'seq1': {'reads': 100},
            'seq2_with_suffix': {'reads': 200},
            'seq3.part1.part2': {'reads': 300}
        }

        # Test with _R_ replacement and additional parts
        assert get_reads('_R_seq1.additional', contig_to_reads) == 100

        # Test with partial match after split
        assert get_reads('seq3.part1.part2.extra', contig_to_reads) == 300

        # Test with combination of replacements
        assert get_reads('_R_seq2_with_suffix.extra', contig_to_reads) == 200

    def test_edge_cases(self):
        """Test edge cases for the get_reads function."""
        contig_to_reads = {
            'empty': {'reads': 0},
            'negative': {'reads': -10},
            'seq.with.dots': {'reads': 500}
        }

        # Test with zero reads
        assert get_reads('empty', contig_to_reads) == 0

        # Test with negative reads (shouldn't occur in real data)
        assert get_reads('negative', contig_to_reads) == -10

        # Test with dots in both name and lookup
        assert get_reads('seq.with.dots.extra', contig_to_reads) == 500

    def test_different_dictionary_structures(self):
        """Test with different contig_to_reads dictionary structures."""
        # Missing 'reads' key
        contig_to_reads = {'seq1': {}}
        with pytest.raises(KeyError):
            get_reads('seq1', contig_to_reads)

        # Direct integer instead of dict
        contig_to_reads = {'seq1': 100}
        with pytest.raises(TypeError):
            get_reads('seq1', contig_to_reads)


class TestClusterSequences:
    """Tests for cluster_sequences function."""

    def test_basic_clustering(self):
        """Test basic clustering functionality."""
        seq_names = ['A', 'B', 'C', 'D', 'E']
        tips_lookup = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}

        # Matrix where A,B,C form one cluster and D,E another
        dist_matrix = np.array([
            [0.0, 0.01, 0.01, 0.5, 0.5],
            [0.01, 0.0, 0.01, 0.5, 0.5],
            [0.01, 0.01, 0.0, 0.5, 0.5],
            [0.5, 0.5, 0.5, 0.0, 0.01],
            [0.5, 0.5, 0.5, 0.01, 0.0]
        ])

        clusters = cluster_sequences(seq_names, tips_lookup, dist_matrix, 0.02)

        # Should form two clusters: [A,B,C] and [D,E]
        assert len(clusters) == 2
        assert set(clusters[0]) == set(['A', 'B', 'C']) or set(clusters[0]) == set(['D', 'E'])
        assert set(clusters[1]) == set(['D', 'E']) or set(clusters[1]) == set(['A', 'B', 'C'])

    def test_empty_input(self):
        """Test with empty input."""
        result = cluster_sequences([], {}, np.array([]), 0.1)
        assert result == []

    def test_single_sequence(self):
        """Test with single sequence."""
        seq_names = ['seqA']
        tips_lookup = {'seqA': 0}
        dist_matrix = np.array([[0.0]])

        result = cluster_sequences(seq_names, tips_lookup, dist_matrix, 0.1)
        assert result == [['seqA']]

    def test_multiple_distinct_clusters(self):
        """Test with multiple completely distinct clusters."""
        seq_names = ['A', 'B', 'C', 'D', 'E', 'F']
        tips_lookup = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}

        # Create matrix with three distinct clusters: [A,B], [C,D], [E,F]
        dist_matrix = np.ones((6, 6))  # Start with all 1s (above threshold)

        # Set within-cluster distances to 0.01 (below threshold)
        dist_matrix[0, 1] = dist_matrix[1, 0] = 0.01  # A-B
        dist_matrix[2, 3] = dist_matrix[3, 2] = 0.01  # C-D
        dist_matrix[4, 5] = dist_matrix[5, 4] = 0.01  # E-F

        # Set diagonal to 0
        np.fill_diagonal(dist_matrix, 0.0)

        clusters = cluster_sequences(seq_names, tips_lookup, dist_matrix, 0.05)

        # Should form three clusters
        assert len(clusters) == 3
        assert sorted([sorted(c) for c in clusters]) == [['A', 'B'], ['C', 'D'], ['E', 'F']]

    def test_threshold_edge_cases(self):
        """Test behavior at the threshold boundary."""
        seq_names = ['A', 'B', 'C']
        tips_lookup = {'A': 0, 'B': 1, 'C': 2}

        # Create matrix with distances exactly at threshold (0.05)
        dist_matrix = np.array([
            [0.0, 0.05, 0.06],
            [0.05, 0.0, 0.05],
            [0.06, 0.05, 0.0]
        ])

        # With threshold 0.05, all sequences should be in same cluster
        # since the algorithm checks distance <= threshold
        clusters = cluster_sequences(seq_names, tips_lookup, dist_matrix, 0.05)
        assert len(clusters) == 1
        assert sorted(clusters[0]) == ['A', 'B', 'C']

        # With threshold 0.04, should form separate clusters
        clusters = cluster_sequences(seq_names, tips_lookup, dist_matrix, 0.04)
        assert len(clusters) == 3


@pytest.fixture
def mock_dependencies():
    """Create mock dependencies for testing find_duplicates and classify_sample."""
    # Create mocks
    mock_tree = MagicMock()
    mock_dist_matrix = np.array([
        [0.0, 0.01, 0.1],
        [0.01, 0.0, 0.1],
        [0.1, 0.1, 0.0]
    ])
    mock_tips_lookup = {'seq1': 0, 'seq2': 1, 'seq3': 2}
    mock_contig_to_reads = {
        'seq1': {'reads': 100},
        'seq2': {'reads': 200},
        'seq3': {'reads': 150}
    }

    # Mock functions
    with patch('lasvdedup.utils.tree_utils.get_mrca_clade') as mock_mrca:
        mock_mrca.return_value = ['seq1', 'seq2', 'seq3', 'seq4', 'seq5']
        yield {
            'tree': mock_tree,
            'dist_matrix': mock_dist_matrix,
            'tips_lookup': mock_tips_lookup,
            'contig_to_reads': mock_contig_to_reads,
            'mock_mrca': mock_mrca
        }


class TestFindDuplicates:
    """Tests for find_duplicates function."""

    @pytest.fixture
    def setup_find_duplicates(self):
        """Set up common test data for find_duplicates tests."""
        # Create a simple distance matrix for 4 sequences
        dist_matrix = np.array([
            [0.0, 0.01, 0.10, 0.50],
            [0.01, 0.0, 0.11, 0.51],
            [0.10, 0.11, 0.0, 0.52],
            [0.50, 0.51, 0.52, 0.0]
        ])

        contig_to_reads = {
            'seqA': {'reads': 100},
            'seqB': {'reads': 200},
            'seqC': {'reads': 150},
            'seqD': {'reads': 300}
        }

        tips = ['seqA', 'seqB', 'seqC', 'seqD']

        mock_tree = MagicMock()

        return {
            'dist_matrix': dist_matrix,
            'contig_to_reads': contig_to_reads,
            'tips': tips,
            'tree': mock_tree
        }

    def test_basic_find_duplicates(self, mock_dependencies):
        """Test basic functionality of find_duplicates."""
        # Mock dependencies
        sample_to_seqs = {
            'sample1': ['seq1', 'seq2'],
            'sample2': ['seq3']
        }

        # Mock classify_sample to return predetermined classifications
        with patch('lasvdedup.utils.sequence_grouping.classify_sample') as mock_classify:
            mock_classify.side_effect = [
                # For sample1
                {
                    'seq1': Classification(
                        'seq1', ClassificationType.GOOD, 'test reason', 'sample1',
                        ['seq1', 'seq2'], DecisionCategory.BELOW_LOWER_THRESHOLD, 100
                    ),
                    'seq2': Classification(
                        'seq2', ClassificationType.BAD, 'test reason', 'sample1',
                        ['seq1', 'seq2'], DecisionCategory.BELOW_LOWER_THRESHOLD, 200
                    )
                },
                # For sample2
                {
                    'seq3': Classification(
                        'seq3', ClassificationType.GOOD, 'test reason', 'sample2',
                        ['seq3'], DecisionCategory.SINGLE_SEQUENCE, 150
                    )
                }
            ]

            # Call find_duplicates
            result = find_duplicates(
                sample_to_seqs, ['seq1', 'seq2', 'seq3'], mock_dependencies['dist_matrix'],
                mock_dependencies['contig_to_reads'], 0.02, 0.05, mock_dependencies['tree']
            )

        # Check results
        assert len(result) == 3
        assert result['seq1'].is_good
        assert result['seq2'].is_bad
        assert result['seq3'].is_good

        # Verify classify_sample was called twice (once per sample)
        assert mock_classify.call_count == 2

    def test_empty_input(self, setup_find_duplicates):
        """Test with empty input."""
        data = setup_find_duplicates

        # Empty sample_to_seqs
        result = find_duplicates(
            {}, data['tips'], data['dist_matrix'],
            data['contig_to_reads'], 0.02, 0.05, data['tree']
        )
        assert result == {}

    def test_non_overlapping_samples(self, setup_find_duplicates):
        """Test with samples that don't overlap."""
        data = setup_find_duplicates

        # Each sequence in a different sample
        sample_to_seqs = {
            'sample1': ['seqA'],
            'sample2': ['seqB'],
            'sample3': ['seqC'],
            'sample4': ['seqD']
        }

        with patch('lasvdedup.utils.sequence_grouping.classify_sample') as mock_classify:
            # Set up return values for each sample
            mock_classify.side_effect = [
                {'seqA': Classification(
                    'seqA', ClassificationType.GOOD, 'Single sequence',
                    'sample1', ['seqA'], DecisionCategory.SINGLE_SEQUENCE, 100
                )},
                {'seqB': Classification(
                    'seqB', ClassificationType.GOOD, 'Single sequence',
                    'sample2', ['seqB'], DecisionCategory.SINGLE_SEQUENCE, 200
                )},
                {'seqC': Classification(
                    'seqC', ClassificationType.GOOD, 'Single sequence',
                    'sample3', ['seqC'], DecisionCategory.SINGLE_SEQUENCE, 150
                )},
                {'seqD': Classification(
                    'seqD', ClassificationType.GOOD, 'Single sequence',
                    'sample4', ['seqD'], DecisionCategory.SINGLE_SEQUENCE, 300
                )}
            ]

            result = find_duplicates(
                sample_to_seqs, data['tips'], data['dist_matrix'],
                data['contig_to_reads'], 0.02, 0.05, data['tree']
            )

        assert len(result) == 4
        assert all(result[seq].is_good for seq in data['tips'])
        assert all(result[seq].decision_category == DecisionCategory.SINGLE_SEQUENCE for seq in data['tips'])

    def test_integration_with_classify_sample(self, setup_find_duplicates):
        """Test integration between find_duplicates and classify_sample."""
        data = setup_find_duplicates

        # Create a sample with seqA and seqB (which are below lower threshold)
        sample_to_seqs = {'sample1': ['seqA', 'seqB']}

        # Don't mock classify_sample to test actual integration
        result = find_duplicates(
            sample_to_seqs, data['tips'], data['dist_matrix'],
            data['contig_to_reads'], 0.02, 0.05, data['tree']
        )

        # Check results - should classify seqB as good (has more reads) and seqA as bad
        assert len(result) == 2
        assert result['seqA'].is_bad
        assert result['seqB'].is_good
        assert result['seqB'].decision_category == DecisionCategory.BELOW_LOWER_THRESHOLD
        assert result['seqA'].decision_category == DecisionCategory.BELOW_LOWER_THRESHOLD


class TestClassifySample:
    """Tests for classify_sample function."""

    @pytest.fixture
    def setup_classify(self):
        """Set up common test data for classify_sample tests."""
        # Create a simple distance matrix for 6 sequences
        dist_matrix = np.array([
            [0.0, 0.01, 0.03, 0.10, 0.50, 0.51],
            [0.01, 0.0, 0.03, 0.10, 0.50, 0.51],
            [0.03, 0.03, 0.0, 0.10, 0.50, 0.51],
            [0.10, 0.10, 0.10, 0.0, 0.50, 0.51],
            [0.50, 0.50, 0.50, 0.50, 0.0, 0.01],
            [0.51, 0.51, 0.51, 0.51, 0.01, 0.0]
        ])

        tips_lookup = {
            'seq1': 0, 'seq2': 1, 'seq3': 2,
            'seq4': 3, 'seq5': 4, 'seq6': 5
        }

        contig_to_reads = {
            'seq1': {'reads': 100},
            'seq2': {'reads': 200},
            'seq3': {'reads': 150},
            'seq4': {'reads': 300},
            'seq5': {'reads': 250},
            'seq6': {'reads': 350}
        }

        mock_tree = MagicMock()

        return {
            'dist_matrix': dist_matrix,
            'tips_lookup': tips_lookup,
            'contig_to_reads': contig_to_reads,
            'tree': mock_tree
        }

    def test_classify_sample_single_sequence(self, mock_dependencies):
        """Test classification of a single sequence."""
        # Call with just one sequence
        result = classify_sample(
            'sample1',
            ['seq1'],
            mock_dependencies['tips_lookup'],
            mock_dependencies['dist_matrix'],
            mock_dependencies['contig_to_reads'],
            0.02, 0.05, mock_dependencies['tree']
        )

        # Check result
        assert len(result) == 1
        assert 'seq1' in result
        assert result['seq1'].is_good
        assert result['seq1'].decision_category == DecisionCategory.SINGLE_SEQUENCE

    def test_classify_sample_below_lower_threshold(self, mock_dependencies):
        """Test classification when all distances are below lower threshold."""
        # Create mock distances
        with patch('lasvdedup.utils.sequence_grouping.get_distances') as mock_get_distances:
            mock_get_distances.return_value = [0.01]  # Below threshold

            # Call the function with seq1 and seq2
            result = classify_sample(
                'sample1',
                ['seq1', 'seq2'],
                mock_dependencies['tips_lookup'],
                mock_dependencies['dist_matrix'],
                mock_dependencies['contig_to_reads'],
                0.02, 0.05, mock_dependencies['tree']
            )

        # Should classify seq2 as good (highest reads) and seq1 as bad
        assert len(result) == 2
        assert result['seq2'].is_good
        assert result['seq2'].decision_category == DecisionCategory.BELOW_LOWER_THRESHOLD
        assert result['seq1'].is_bad
        assert result['seq1'].decision_category == DecisionCategory.BELOW_LOWER_THRESHOLD

    def test_classify_sample_below_upper_threshold(self, mock_dependencies):
        """Test classification when all distances are below upper threshold."""
        # Create mock distances
        with patch('lasvdedup.utils.sequence_grouping.get_distances') as mock_get_distances:
            # All distances are below upper threshold (0.05) but above lower (0.02)
            mock_get_distances.return_value = [0.03]

            # Mock clustering to return two clusters
            with patch('lasvdedup.utils.sequence_grouping.cluster_sequences') as mock_cluster:
                # Each sequence is its own cluster (crucial for them to be marked as coinfection)
                mock_cluster.return_value = [['seq1'], ['seq2']]

                result = classify_sample(
                    'sample1',
                    ['seq1', 'seq2'],
                    mock_dependencies['tips_lookup'],
                    mock_dependencies['dist_matrix'],
                    mock_dependencies['contig_to_reads'],
                    0.02, 0.05, None  # Don't need tree for this path
                )

        # Should classify each sequence as coinfection (as each is their own cluster)
        assert len(result) == 2
        assert result['seq1'].is_coinfection
        assert result['seq2'].is_coinfection
        assert result['seq1'].decision_category == DecisionCategory.BELOW_UPPER_THRESHOLD
        assert result['seq2'].decision_category == DecisionCategory.BELOW_UPPER_THRESHOLD

    def test_multiple_clusters_below_lower_threshold(self, setup_classify):
        """Test scenario with multiple clusters below lower threshold."""
        data = setup_classify

        # Create a custom distance matrix where all distances are below upper threshold
        # but we have two distinct clusters based on lower threshold
        # Cluster 1: seq1, seq2 (distance < 0.02)
        # Cluster 2: seq5, seq6 (distance < 0.02)
        # But distances between clusters are < 0.05 (upper threshold)
        test_dist_matrix = np.array([
            [0.0, 0.01, 0.04, 0.04, 0.04, 0.04],  # seq1
            [0.01, 0.0, 0.04, 0.04, 0.04, 0.04],  # seq2
            [0.04, 0.04, 0.0, 0.01, 0.04, 0.04],  # seq3
            [0.04, 0.04, 0.01, 0.0, 0.04, 0.04],  # seq4
            [0.04, 0.04, 0.04, 0.04, 0.0, 0.01],  # seq5
            [0.04, 0.04, 0.04, 0.04, 0.01, 0.0]   # seq6
        ])

        # Test only with seq1, seq2, seq5, seq6
        test_sequences = ['seq1', 'seq2', 'seq5', 'seq6']

        # This is the key part - we're ensuring the test only looks at the case where
        # all distances are below upper threshold, so tree-related code won't be called
        with patch('lasvdedup.utils.sequence_grouping.get_distances') as mock_distances:
            # All distances are below upper threshold (0.05)
            mock_distances.return_value = [0.01, 0.04, 0.04, 0.04, 0.01, 0.04]

            result = classify_sample(
                'sample1',
                test_sequences,
                data['tips_lookup'],
                test_dist_matrix,  # Use our customized distance matrix
                data['contig_to_reads'],
                0.02,
                0.05,
                None,
            )

        # Should identify as potential intrahost variants
        assert len(result) == 4

        # Best sequence from each cluster should be marked as coinfection
        good_seqs = [seq for seq, cls in result.items() if cls.is_coinfection]
        assert len(good_seqs) == 2
        assert 'seq2' in good_seqs  # Highest reads in first cluster
        assert 'seq6' in good_seqs  # Highest reads in second cluster

        # Other sequences should be marked as bad
        assert result['seq1'].is_bad
        assert result['seq5'].is_bad

        # Decision category should be BELOW_UPPER_THRESHOLD
        for seq in test_sequences:
            assert result[seq].decision_category == DecisionCategory.BELOW_UPPER_THRESHOLD

    def test_outlier_detection(self, setup_classify):
        """Test outlier detection logic."""
        data = setup_classify

        # Create a scenario with distances above thresholds
        with patch('lasvdedup.utils.sequence_grouping.get_distances') as mock_distances:
            # Make sure we're triggering the path for distances above threshold
            mock_distances.return_value = [0.01, 0.06, 0.06, 0.06, 0.01, 0.06]

            # Need to patch the exact import location as used in the sequence_grouping module
            with patch('lasvdedup.utils.sequence_grouping.get_mrca_clade') as mock_mrca:
                # Set the return value for the MRCA clade
                mock_mrca.return_value = ['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6', 'otherseq']

                # Also patch get_outliers
                with patch('lasvdedup.utils.sequence_grouping.get_outliers') as mock_outliers:
                    # seq4 is an outlier
                    mock_outliers.return_value = ['seq4']

                    result = classify_sample(
                        'sample1',
                        ['seq1', 'seq2', 'seq3', 'seq4'],
                        data['tips_lookup'],
                        data['dist_matrix'],
                        data['contig_to_reads'],
                        0.02, 0.05, data['tree'],
                        min_clade_size=5
                    )

        # Should identify seq4 as outlier and mark it as BAD
        assert result['seq4'].is_bad
        assert result['seq4'].decision_category == DecisionCategory.OUTLIERS_DETECTED

        # Should select sequence with highest read count among non-outliers as GOOD
        non_outliers = ['seq1', 'seq2', 'seq3']
        best_seq = max(non_outliers, key=lambda seq: data['contig_to_reads'][seq]['reads'])
        assert result[best_seq].is_good
        assert result[best_seq].decision_category == DecisionCategory.OUTLIERS_DETECTED

    def test_true_inrahost_coinfection(self, setup_classify):
        """Test true coinfection scenario."""
        data = setup_classify

        dist_matrix = np.array([
            [0.0, 0.06, 0.06],
            [0.06, 0.0, 0.06],
            [0.06, 0.06, 0.0],
        ])

        # Need to patch the exact import location as used in the sequence_grouping module
        with patch('lasvdedup.utils.sequence_grouping.get_mrca_clade') as mock_mrca:
            # Set the return value to a large enough clade
            mock_mrca.return_value = ['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6', 'otherseq1', 'otherseq2']

            # Also patch get_outliers
            with patch('lasvdedup.utils.sequence_grouping.get_outliers') as mock_outliers:
                # No outliers
                mock_outliers.return_value = []

                result = classify_sample(
                    'sample1',
                    ['seq1', 'seq2', 'seq4'],
                    data['tips_lookup'],
                    data['dist_matrix'],
                    data['contig_to_reads'],
                    0.02, 0.09, data['tree'],
                    min_clade_size=5
                )

        # All sequences should be marked as coinfection
        assert all(cls.is_coinfection for cls in result.values())
        assert all(cls.decision_category == DecisionCategory.TRUE_COINFECTION for cls in result.values())

    def test_true_diffstrain_coinfection(self, setup_classify):
        """Test true coinfection scenario."""
        data = setup_classify

        dist_matrix = np.array([
            [0.0, 0.06, 0.06],
            [0.06, 0.0, 0.06],
            [0.06, 0.06, 0.0],
        ])

        # Need to patch the exact import location as used in the sequence_grouping module
        with patch('lasvdedup.utils.sequence_grouping.get_mrca_clade') as mock_mrca:
            # Set the return value to a large enough clade
            mock_mrca.return_value = ['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6', 'otherseq1', 'otherseq2']

            # Also patch get_outliers
            with patch('lasvdedup.utils.sequence_grouping.get_outliers') as mock_outliers:
                # No outliers
                mock_outliers.return_value = []

                result = classify_sample(
                    'sample1',
                    ['seq1', 'seq2', 'seq4'],
                    data['tips_lookup'],
                    data['dist_matrix'],
                    data['contig_to_reads'],
                    0.02, 0.05, data['tree'],
                    min_clade_size=5
                )

        # All sequences should be marked as coinfection
        assert all(cls.is_coinfection for cls in result.values())
        assert all(cls.decision_category == DecisionCategory.TRUE_COINFECTION for cls in result.values())

    def test_small_clade(self, setup_classify):
        """Test small clade scenario."""
        data = setup_classify

        # Mock the get_mrca_clade function to return small clade
        with patch('lasvdedup.utils.tree_utils.get_mrca_clade') as mock_mrca:
            mock_mrca.return_value = ['seq1', 'seq2', 'seq4']  # Only 3 sequences

            result = classify_sample(
                'sample1',
                ['seq1', 'seq2', 'seq4'],  # These have distances above thresholds between them
                data['tips_lookup'],
                data['dist_matrix'],
                data['contig_to_reads'],
                0.02, 0.05, data['tree'],
                min_clade_size=5  # Clade size (3) is below this
            )

        # Should select seq4 as good (highest read count) and mark others as bad
        assert result['seq4'].is_good
        assert result['seq1'].is_bad
        assert result['seq2'].is_bad

        assert all(cls.decision_category == DecisionCategory.SMALL_CLADE for cls in result.values())


# Run the tests
if __name__ == "__main__":
    pytest.main(["-v", __file__])
