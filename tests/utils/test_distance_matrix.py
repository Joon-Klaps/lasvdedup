import os
import pytest
import tempfile
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock, patch

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
    clade = MagicMock()
    seq_names = ['A', 'B', 'C']

    # Setup mock depths for the terminals
    terminals = [MagicMock() for _ in range(5)]
    clade.get_terminals.return_value = terminals

    # Create mock clade.depths() return value where C is an outlier
    clade_depths = {
        terminals[0]: 0.2,  # A
        terminals[1]: 0.2,  # B
        terminals[2]: 0.5,  # C is an outlier
        terminals[3]: 0.2,  # D
        terminals[4]: 0.2,  # E
    }
    clade.depths.return_value = clade_depths

    # Setup mock find_any to return the correct terminal for each seq name
    clade.find_any.side_effect = lambda name: {
        'A': terminals[0],
        'B': terminals[1],
        'C': terminals[2],
    }[name]

    # Call the function with evolution_threshold and z_threshold
    evolution_threshold = 0.02
    z_threshold = 2.0
    outliers = get_outliers(clade, seq_names, evolution_threshold, z_threshold)

    # Verify results
    assert isinstance(outliers, dict), "get_outliers should return a dictionary"
    assert 'C' in outliers, "C should be identified as an outlier"
    assert 'A' not in outliers, "A should not be identified as an outlier"
    assert 'B' not in outliers, "B should not be identified as an outlier"

    # Verify the outlier info contains the required fields
    assert 'distance' in outliers['C'], "Outlier info should contain distance"
    assert 'median' in outliers['C'], "Outlier info should contain median"
    assert 'threshold' in outliers['C'], "Outlier info should contain threshold"

    # Test with a higher evolution-threshold that won't detect outliers
    outliers = get_outliers(clade, seq_names, evolution_threshold=0.4, z_threshold=2)
    assert len(outliers) == 0, "With higher evolution threshold, should find no outliers"

    # Test with MAD = 0 case
    uniform_clade = MagicMock()
    uniform_terminals = [MagicMock() for _ in range(5)]
    uniform_clade.get_terminals.return_value = uniform_terminals

    # All terminals have same depth (MAD will be 0)
    uniform_depths = {terminal: 0.1 for terminal in uniform_terminals}
    uniform_clade.depths.return_value = uniform_depths

    # Make one sequence slightly higher to test threshold
    uniform_depths[uniform_terminals[2]] = 0.15  # Above threshold with evolution_threshold = 0.02

    uniform_clade.find_any.side_effect = lambda name: {
        'A': uniform_terminals[0],
        'B': uniform_terminals[1],
        'C': uniform_terminals[2],
    }[name]

    # With evolution_threshold, C should be detected
    outliers = get_outliers(uniform_clade, seq_names, evolution_threshold=0.02, z_threshold=2.0)
    assert 'C' in outliers, "With MAD=0 and evolution_threshold=0.02, C should be identified as outlier"
