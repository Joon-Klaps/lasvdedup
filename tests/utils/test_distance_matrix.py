import os
import pytest
import tempfile
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock, patch

from lasvdedup.utils.distance_matrix import (
    get_distances,
    get_outliers,
    to_distance_matrix,
    calculate_dist_bootstrap
)

class MockPhyloDM:
    """Mock PhyloDM class for testing."""
    def __init__(self, taxa, distances):
        self._taxa = taxa
        self._distances = distances

    def dm(self, norm=False):
        return self._distances

    def taxa(self):
        return self._taxa


def test_to_distance_matrix():
    """Test conversion of PhyloDM object to distance matrix."""
    # Sample data
    taxa = ['A', 'B', 'C']
    distances = np.array([
        [0.0, 0.1, 0.2],
        [0.1, 0.0, 0.3],
        [0.2, 0.3, 0.0]
    ])

    # Create a mock PhyloDM object using the MockPhyloDM class
    mock_tree = MockPhyloDM(taxa, distances)

    # Call the function
    tips, matrix = to_distance_matrix(mock_tree)

    # Verify results
    assert tips == ['A', 'B', 'C']
    assert matrix.shape == (3, 3)
    assert matrix[0, 0] == 0.0  # Distance to self is zero
    assert matrix[0, 1] == 0.1  # A to B
    assert matrix[0, 2] == 0.2  # A to C
    assert matrix[1, 2] == 0.3  # B to C
    assert matrix[1, 0] == 0.1  # B to A (symmetric)
    assert matrix[2, 0] == 0.2  # C to A (symmetric)
    assert matrix[2, 1] == 0.3  # C to B (symmetric)

def test_get_distances():
    """Test getting pairwise distances between sequences."""
    names = ['A', 'B', 'C']
    tips_lookup = {'A': 0, 'B': 1, 'C': 2}
    matrix = np.array([
        [0.0, 0.1, 0.2],
        [0.1, 0.0, 0.3],
        [0.2, 0.3, 0.0]
    ])

    # Call the function
    distances = get_distances(names, tips_lookup, matrix)

    # There should be 3 pairwise distances
    assert len(distances) == 3
    assert 0.1 in distances  # A-B
    assert 0.2 in distances  # A-C
    assert 0.3 in distances  # B-C

@patch('numpy.random.choice')
def test_calculate_dist_bootstrap(mock_choice):
    """Test bootstrap calculation for distance variance."""
    # Mock the random sampling to return predictable values
    mock_choice.return_value = np.array([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]])

    distances = np.array([1.0, 2.0, 3.0])
    result = calculate_dist_bootstrap(distances, n_bootstraps=2, quantile=0.95)

    # Based on our mock, means should be [2.0, 3.0], and 95th percentile would be 2.95
    assert result == 2.95

def test_get_outliers():
    """Test identifying outlier sequences."""
    # Create sample data
    clade_members = ['A', 'B', 'C', 'D', 'E']
    seq_names = ['A', 'B', 'C']  # We only check these
    tips_lookup = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}

    # Create a distance matrix where C is an outlier
    dist_matrix = np.array([
        [0.0, 0.2, 0.5, 0.2, 0.2],
        [0.2, 0.0, 0.5, 0.1, 0.1],
        [0.5, 0.5, 0.0, 0.5, 0.5],
        [0.2, 0.1, 0.5, 0.0, 0.1],
        [0.2, 0.1, 0.5, 0.1, 0.0]
    ])

    # Call the function
    outliers = get_outliers(clade_members, seq_names, tips_lookup, dist_matrix, quantile=0.85)

    # Verify results
    assert isinstance(outliers, dict), "get_outliers should return a dictionary"
    assert 'C' in outliers, "C should be identified as an outlier"
    assert 'A' not in outliers, "A should not be identified as an outlier"
    assert 'B' not in outliers, "B should not be identified as an outlier"

    # Verify the outlier info contains the required fields
    assert 'distance' in outliers['C'], "Outlier info should contain distance"
    assert 'median' in outliers['C'], "Outlier info should contain median"
    assert 'threshold' in outliers['C'], "Outlier info should contain threshold"
    assert 'reference' in outliers['C'], "Outlier info should contain reference"

    # Test with a lower quantile that includes more outliers
    outliers = get_outliers(clade_members, seq_names, tips_lookup, dist_matrix, quantile=0.60)
    assert len(outliers) > 0, "With lower z-threshold, should find at least one outlier"

    # Test with all identical distances (no outliers)
    uniform_matrix = np.ones((5, 5)) * 0.1
    np.fill_diagonal(uniform_matrix, 0.0)

    # Test the MAD = 0 case using a uniform matrix where bootstrap is needed
    with patch('lasvdedup.utils.distance_matrix.calculate_dist_bootstrap', return_value=0.15):
        outliers = get_outliers(clade_members, seq_names, tips_lookup, uniform_matrix)
        assert outliers == {}, "With uniform distances and threshold=0.15, should find no outliers"

    # Test the MAD = 0 case with a higher bootstrap value that would cause outliers
    with patch('lasvdedup.utils.distance_matrix.calculate_dist_bootstrap', return_value=0.05):
        # Setting up a case where some distances exceed the bootstrap threshold
        test_matrix = uniform_matrix.copy()
        test_matrix[3, 2] = 0.12  # Distance from D (ref) to C is 0.12, above threshold of 0.05
        test_matrix[2, 3] = 0.12  # Make symmetric

        outliers = get_outliers(clade_members, seq_names, tips_lookup, test_matrix)
        assert 'C' in outliers, "With MAD=0 and threshold=0.05, C should be identified as outlier"
