import pytest
import numpy as np
from unittest.mock import MagicMock, patch
import logging
import pprint

from lasvdedup.utils.sequence_grouping import (
    group_sequences_by_sample,
    is_all_distances_below_threshold,
    select_best_sequence,
    get_contig_data,  # Changed from get_reads
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
        seq_data = {
            'seqA': {'rank': 2, 'reads': 100},
            'seqB': {'rank': 1, 'reads': 200},
            'seqC': {'rank': 3, 'reads': 50}
        }

        best_seq = select_best_sequence(seq_names, seq_data)
        assert best_seq == 'seqB'  # seqB has the lowest rank (1)

        # Test with tie (should take the first encountered)
        seq_data = {
            'seqA': {'rank': 1, 'reads': 200},
            'seqB': {'rank': 1, 'reads': 200},
            'seqC': {'rank': 3, 'reads': 50}
        }
        best_seq = select_best_sequence(seq_names, seq_data)
        assert best_seq in ['seqA', 'seqB']  # Either could be selected in case of tie

    def test_with_empty_sequence_list(self):
        """Test behavior with empty sequence list."""
        with pytest.raises(ValueError):
            select_best_sequence([], {'seqA': {'rank': 1}})

    def test_with_missing_rank_data(self):
        """Test behavior when some sequences don't have rank data."""
        seq_names = ['seqA', 'seqB', 'seqC']
        seq_data = {
            'seqA': {'rank': 1},
            'seqC': {'rank': 3}  # seqB missing
        }

        with pytest.raises(KeyError):
            select_best_sequence(seq_names, seq_data)

    def test_with_infinity_rank(self):
        """Test with infinite rank value (edge case)."""
        seq_names = ['seqA', 'seqB', 'seqC']
        seq_data = {
            'seqA': {'rank': 100},
            'seqB': {'rank': float('-inf')},  # Lower rank is better
            'seqC': {'rank': 50}
        }

        best_seq = select_best_sequence(seq_names, seq_data)
        assert best_seq == 'seqB'


class TestGetContigData:
    """Tests for get_contig_data function."""

    def test_basic_contig_data(self):
        """Test basic contig data retrieval."""
        contigs_ranked = {
            'seq1': {'reads': 100, 'coverage': 95.5, 'rank': 2},
            'seq2': {'reads': 200, 'coverage': 98.0, 'rank': 1},
            'seq3': {'reads': 300, 'coverage': 97.0, 'rank': 3}
        }

        # Direct match
        assert get_contig_data('seq1', contigs_ranked) == {'reads': 100, 'coverage': 95.5, 'rank': 2}

        # Match after _R_ replacement
        assert get_contig_data('_R_seq2', contigs_ranked)['reads'] == 200

        # Match with splitting
        assert get_contig_data('seq3.suffix', contigs_ranked)['rank'] == 3

        # Not found
        with pytest.raises(ValueError):
            get_contig_data('notfound', contigs_ranked)

    def test_with_complex_name_handling(self):
        """Test more complex sequence name handling scenarios."""
        contigs_ranked = {
            'seq1': {'reads': 100, 'rank': 3},
            'seq2_with_suffix': {'reads': 200, 'rank': 1},
            'seq3.part1': {'reads': 300, 'rank': 2}
        }

        # Test with _R_ replacement
        result = get_contig_data('_R_seq1', contigs_ranked)
        assert result['reads'] == 100
        assert result['rank'] == 3

        # Test with partial match after split
        result = get_contig_data('seq3.part1.extra', contigs_ranked)
        assert result['reads'] == 300
        assert result['rank'] == 2

        # Test with combination of replacements
        result = get_contig_data('_R_seq2_with_suffix', contigs_ranked)
        assert result['reads'] == 200
        assert result['rank'] == 1

    def test_edge_cases(self):
        """Test edge cases for the get_contig_data function."""
        contigs_ranked = {
            'empty': {'reads': 0, 'rank': 3},
            'negative': {'reads': -10, 'rank': 2},
            'seq.with.dots': {'reads': 500, 'rank': 1}
        }

        # Test with zero reads
        result = get_contig_data('empty', contigs_ranked)
        assert result['reads'] == 0
        assert result['rank'] == 3

        # Test with negative reads
        result = get_contig_data('negative', contigs_ranked)
        assert result['reads'] == -10
        assert result['rank'] == 2

        # Test with dots in name
        result = get_contig_data('seq.with.dots.extra', contigs_ranked)
        assert result['reads'] == 500
        assert result['rank'] == 1


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
    mock_contigs_ranked = {
        'seq1': {'reads': 100, 'coverage': 95.0, 'rank': 3},
        'seq2': {'reads': 200, 'coverage': 98.0, 'rank': 1},
        'seq3': {'reads': 150, 'coverage': 96.0, 'rank': 2}
    }

    mock_thresholds = {
        "PWD": 0.02,
        "CLADE_SIZE": 5,
        "Z_THRESHOLD": 3.0
    }

    # Mock functions
    with patch('lasvdedup.utils.tree_utils.get_mrca_clade') as mock_mrca:
        mock_mrca.return_value = ['seq1', 'seq2', 'seq3', 'seq4', 'seq5']
        yield {
            'tree': mock_tree,
            'dist_matrix': mock_dist_matrix,
            'tips_lookup': mock_tips_lookup,
            'contigs_ranked': mock_contigs_ranked,
            'mock_mrca': mock_mrca,
            'thresholds': mock_thresholds
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

        contigs_ranked = {
            'seqA': {'reads': 100, 'coverage': 95.0, 'rank': 3},
            'seqB': {'reads': 200, 'coverage': 98.0, 'rank': 1},
            'seqC': {'reads': 150, 'coverage': 96.0, 'rank': 2},
            'seqD': {'reads': 300, 'coverage': 97.0, 'rank': 4}
        }

        tips = ['seqA', 'seqB', 'seqC', 'seqD']

        mock_tree = MagicMock()

        thresholds = {
            "PWD": 0.02,
            "CLADE_SIZE": 5,
            "Z_THRESHOLD": 3.0
        }

        return {
            'dist_matrix': dist_matrix,
            'contigs_ranked': contigs_ranked,
            'tips': tips,
            'tree': mock_tree,
            'thresholds': thresholds
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
                        ['seq1', 'seq2'], DecisionCategory.BELOW_THRESHOLD,
                        {'reads': 100, 'rank': 3}
                    ),
                    'seq2': Classification(
                        'seq2', ClassificationType.BAD, 'test reason', 'sample1',
                        ['seq1', 'seq2'], DecisionCategory.BELOW_THRESHOLD,
                        {'reads': 200, 'rank': 1}
                    )
                },
                # For sample2
                {
                    'seq3': Classification(
                        'seq3', ClassificationType.GOOD, 'test reason', 'sample2',
                        ['seq3'], DecisionCategory.SINGLE_SEQUENCE,
                        {'reads': 150, 'rank': 2}
                    )
                }
            ]

            # Call find_duplicates with updated parameters
            result = find_duplicates(
                sample_to_seqs, ['seq1', 'seq2', 'seq3'],
                mock_dependencies['dist_matrix'],
                mock_dependencies['contigs_ranked'],
                mock_dependencies['tree'],
                segment="L",
                thresholds=mock_dependencies['thresholds']
            )

        # Check results
        assert len(result) == 3
        assert result['seq1'].is_good
        assert result['seq2'].is_bad
        assert result['seq3'].is_good

        # Verify classify_sample was called correctly
        assert mock_classify.call_count == 2
        # Check that mock_classify was called with the correct parameters
        mock_classify.assert_any_call(
            'sample1', ['seq1', 'seq2'],
            mock_dependencies['tips_lookup'],
            mock_dependencies['dist_matrix'],
            mock_dependencies['contigs_ranked'],
            mock_dependencies['tree'],
            "L", mock_dependencies['thresholds']
        )

    def test_empty_input(self, setup_find_duplicates):
        """Test with empty input."""
        data = setup_find_duplicates

        # Empty sample_to_seqs
        result = find_duplicates(
            {}, data['tips'], data['dist_matrix'],
            data['contigs_ranked'], data['tree'],
            segment="L", thresholds=data['thresholds']
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
                    'sample1', ['seqA'], DecisionCategory.SINGLE_SEQUENCE,
                    {'reads': 100, 'rank': 3}
                )},
                {'seqB': Classification(
                    'seqB', ClassificationType.GOOD, 'Single sequence',
                    'sample2', ['seqB'], DecisionCategory.SINGLE_SEQUENCE,
                    {'reads': 200, 'rank': 1}
                )},
                {'seqC': Classification(
                    'seqC', ClassificationType.GOOD, 'Single sequence',
                    'sample3', ['seqC'], DecisionCategory.SINGLE_SEQUENCE,
                    {'reads': 150, 'rank': 2}
                )},
                {'seqD': Classification(
                    'seqD', ClassificationType.GOOD, 'Single sequence',
                    'sample4', ['seqD'], DecisionCategory.SINGLE_SEQUENCE,
                    {'reads': 300, 'rank': 4}
                )}
            ]

            result = find_duplicates(
                sample_to_seqs, data['tips'], data['dist_matrix'],
                data['contigs_ranked'], data['tree'],
                segment="L", thresholds=data['thresholds']
            )

        assert len(result) == 4
        assert all(result[seq].is_good for seq in data['tips'])
        assert all(result[seq].decision_category == DecisionCategory.SINGLE_SEQUENCE for seq in data['tips'])


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

        contigs_ranked = {
            'seq1': {'reads': 100, 'coverage': 95.0, 'rank': 4},
            'seq2': {'reads': 200, 'coverage': 98.0, 'rank': 1},
            'seq3': {'reads': 150, 'coverage': 96.0, 'rank': 3},
            'seq4': {'reads': 300, 'coverage': 97.0, 'rank': 2},
            'seq5': {'reads': 250, 'coverage': 94.0, 'rank': 6},
            'seq6': {'reads': 350, 'coverage': 99.0, 'rank': 5}
        }

        mock_tree = MagicMock()

        thresholds = {
            "PWD": 0.02,
            "CLADE_SIZE": 5,
            "Z_THRESHOLD": 3.0
        }

        return {
            'dist_matrix': dist_matrix,
            'tips_lookup': tips_lookup,
            'contigs_ranked': contigs_ranked,
            'tree': mock_tree,
            'thresholds': thresholds
        }

    def test_classify_sample_single_sequence(self, mock_dependencies):
        """Test classification of a single sequence."""
        # Call with just one sequence
        result = classify_sample(
            'sample1',
            ['seq1'],
            mock_dependencies['tips_lookup'],
            mock_dependencies['dist_matrix'],
            mock_dependencies['contigs_ranked'],
            mock_dependencies['tree'],
            segment="L",
            thresholds=mock_dependencies['thresholds']
        )

        # Check result
        assert len(result) == 1
        assert 'seq1' in result
        assert result['seq1'].is_good
        assert result['seq1'].decision_category == DecisionCategory.SINGLE_SEQUENCE

    def test_classify_sample_below_pwd_threshold(self, mock_dependencies):
        """Test classification when all distances are below pwd threshold."""
        # Create mock distances
        with patch('lasvdedup.utils.sequence_grouping.get_distances') as mock_get_distances:
            mock_get_distances.return_value = [0.01]  # Below threshold

            # Call the function with seq1 and seq2
            result = classify_sample(
                'sample1',
                ['seq1', 'seq2'],
                mock_dependencies['tips_lookup'],
                mock_dependencies['dist_matrix'],
                mock_dependencies['contigs_ranked'],
                mock_dependencies['tree'],
                segment="L",
                thresholds=mock_dependencies['thresholds']
            )

        # Should classify seq2 as good (highest reads) and seq1 as bad
        assert len(result) == 2
        assert result['seq2'].is_good
        assert result['seq2'].decision_category == DecisionCategory.BELOW_THRESHOLD
        assert result['seq1'].is_bad
        assert result['seq1'].decision_category == DecisionCategory.BELOW_THRESHOLD

    def test_multiple_clusters_below_pwd_threshold(self, setup_classify):
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
                data['contigs_ranked'],
                data['tree'],
                segment="L",
                thresholds=data['thresholds']
            )

        # Should identify as potential intrahost variants
        assert len(result) == 4

        pprint.pprint(result)

        # Best sequence as all others are below threshold
        good_seqs = [seq for seq, cls in result.items() if cls.is_good]
        assert len(good_seqs) == 1
        assert 'seq2' in good_seqs  # Highest rank in first cluster

        # Other sequences should be marked as bad
        assert result['seq1'].is_bad
        assert result['seq5'].is_bad

        # Decision category should be BELOW_THRESHOLD
        for seq in test_sequences:
            assert result[seq].decision_category == DecisionCategory.SMALL_CLADE

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
                    mock_outliers.return_value =  { "seq4" :
                        {'distance': float(0.06),
                        'median': float(0.04),
                        'threshold': float(0.05),
                        'reference': 'seq1',
                    }}

                    result = classify_sample(
                        'sample1',
                        ['seq1', 'seq2', 'seq3', 'seq4'],
                        data['tips_lookup'],
                        data['dist_matrix'],
                        data['contigs_ranked'],
                        data['tree'],
                        segment="L",
                        thresholds=data['thresholds']
                    )

                    print(result)

        # Should identify seq4 as outlier and mark it as BAD
        assert result['seq4'].is_bad
        assert result['seq4'].decision_category == DecisionCategory.OUTLIERS_DETECTED

        # Should select sequence with highest read count among non-outliers as GOOD
        non_outliers = ['seq1', 'seq2', 'seq3']
        best_seq = max(non_outliers, key=lambda seq: data['contigs_ranked'][seq]['reads'])
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
                    data['contigs_ranked'],
                    data['tree'],
                    segment="L",
                    thresholds=data['thresholds']
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
                    data['contigs_ranked'],
                    data['tree'],
                    segment="L",
                    thresholds=data['thresholds']
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
                data['contigs_ranked'],
                data['tree'],
                segment="L",
                thresholds=data['thresholds']
            )

        # Should select seq2 as good (highest rank) and mark others as bad
        assert result['seq2'].is_good
        assert result['seq1'].is_bad
        assert result['seq4'].is_bad

        assert all(cls.decision_category == DecisionCategory.SMALL_CLADE for cls in result.values())


# Run the tests
if __name__ == "__main__":
    pytest.main(["-v", __file__])
