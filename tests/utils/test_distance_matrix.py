import os
import pytest
import tempfile
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock

from lasvdedup.utils.distance_matrix import (
    get_distances,
    get_outliers,
    to_distance_matrix
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
    outliers = get_outliers(clade_members, seq_names, tips_lookup, dist_matrix, z_threshold=2.0)

    # Verify results
    assert outliers == ['C']

    # Test with a lower z_threshold that includes more outliers
    outliers = get_outliers(clade_members, seq_names, tips_lookup, dist_matrix, z_threshold=1.0)
    assert len(outliers) > 0

    # Test with all identical distances (no outliers)
    uniform_matrix = np.ones((5, 5)) * 0.1
    np.fill_diagonal(uniform_matrix, 0.0)
    outliers = get_outliers(clade_members, seq_names, tips_lookup, uniform_matrix)
    assert outliers == []
