import os
import tempfile
import pytest
import yaml
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import logging

from lasvdedup.utils.determine_duplicates import determine_duplicates
from lasvdedup.utils.classification import Classification, ClassificationType, DecisionCategory


@pytest.fixture
def mock_config():
    """Create a mock config dictionary with proper paths."""
    return {
        "tree": str(Path("tree_path.treefile")),
        "sequences": str(Path("sequences.fasta")),
        "table": str(Path("contigs.tsv")),
        "prefix": str(Path("output_directory")),
        "species": "LASV",
        "segment": "L",
        "sample_regex": r"(\w+)_.*",
        "length_column": "length",
        "selection_column": ["coverage"],
        "log_level": "INFO",
        "DEDUPLICATE": {
            "THRESHOLDS": {
                "L": {
                    "lower_threshold": 0.01,
                    "upper_threshold": 0.03,
                    "clade_size": 15,
                    "z_threshold": 2.5,
                    "target_length": 7000
                },
                "S": {
                    "lower_threshold": 0.03,
                    "upper_threshold": 0.07,
                    "clade_size": 8,
                    "z_threshold": 1.8,
                    "target_length": 3500
                }
            }
        }
    }


@pytest.fixture
def mock_tree():
    """Create a mock tree object."""
    return {'biophylo': MagicMock(), 'phylodm': MagicMock()}


@pytest.fixture
def mock_dist_matrix_data():
    """Create mock distance matrix data."""
    tips = ['A', 'B', 'C', 'D', 'E']
    dist_matrix = np.zeros((5, 5))
    # Add some distance values
    dist_matrix[0, 1] = dist_matrix[1, 0] = 0.02  # A-B distance
    dist_matrix[0, 2] = dist_matrix[2, 0] = 0.05  # A-C distance
    dist_matrix[3, 4] = dist_matrix[4, 3] = 0.01  # D-E distance
    return tips, dist_matrix


@pytest.fixture
def mock_contigs_ranked():
    """Create mock ranked contigs data."""
    return {
        'A': {'length': 1000, 'coverage': 20.5, 'rank': 1, 'distance_to_expectation': 6000},
        'B': {'length': 950, 'coverage': 15.2, 'rank': 2, 'distance_to_expectation': 6050},
        'C': {'length': 1200, 'coverage': 18.7, 'rank': 3, 'distance_to_expectation': 5800},
        'D': {'length': 800, 'coverage': 22.3, 'rank': 4, 'distance_to_expectation': 6200},
        'E': {'length': 1100, 'coverage': 19.8, 'rank': 5, 'distance_to_expectation': 5900}
    }


@pytest.fixture
def mock_sample_groups():
    """Create mock sample group data."""
    return {'sample1': ['A', 'B'], 'sample2': ['C', 'D', 'E']}


@pytest.fixture
def mock_classifications():
    """Create mock classification results."""
    return {
        'A': Classification(
            sequence_name='A',
            classification_type=ClassificationType.GOOD,
            reason='Test reason',
            sample_id='sample1',
            group_members=['A', 'B'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            contig_stats={'length': 1000, 'coverage': 20.5}
        ),
        'B': Classification(
            sequence_name='B',
            classification_type=ClassificationType.BAD,
            reason='Test reason',
            sample_id='sample1',
            group_members=['A', 'B'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            contig_stats={'length': 950, 'coverage': 15.2}
        )
    }


@pytest.fixture
def mock_seq_records():
    """Create mock sequence records."""
    mock_records = {}
    for name in ['A', 'B', 'C', 'D', 'E']:
        record = MagicMock()
        record.id = name
        record.seq = MagicMock()
        record.seq.__len__.return_value = 1000
        mock_records[name] = record
    return mock_records


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.sort_table')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
@patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
@patch('os.makedirs')
@patch('lasvdedup.utils.determine_duplicates.setup_logging')  # Patch the entire setup_logging function
def test_determine_duplicates_with_mocks(
    mock_setup_logging,
    mock_makedirs,
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_sort_table,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_config,
    mock_tree,
    mock_dist_matrix_data,
    mock_contigs_ranked,
    mock_sample_groups,
    mock_classifications,
    mock_seq_records
):
    """Test determine_duplicates with mocked dependencies."""
    # Set up mocks
    mock_root_tree.return_value = mock_tree
    mock_to_distance_matrix.return_value = mock_dist_matrix_data
    mock_sort_table.return_value = mock_contigs_ranked
    mock_group_sequences.return_value = mock_sample_groups
    mock_find_duplicates.return_value = mock_classifications

    # Run the function
    with patch('Bio.SeqIO.to_dict', return_value=mock_seq_records):
        with patch('Bio.SeqIO.parse'):  # Mock the actual parsing
            result = determine_duplicates(config=mock_config)

    # Assertions
    assert result == mock_classifications
    mock_makedirs.assert_called_once()
    mock_root_tree.assert_called_once_with(mock_config["tree"])
    mock_to_distance_matrix.assert_called_once_with(mock_tree['phylodm'])
    mock_sort_table.assert_called_once()
    mock_group_sequences.assert_called_once_with(mock_dist_matrix_data[0], mock_config["sample_regex"])
    mock_find_duplicates.assert_called_once()
    mock_write_distance_matrix.assert_called_once()
    mock_write_results.assert_called_once()


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.sort_table')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
@patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
@patch('os.makedirs')
@patch('lasvdedup.utils.determine_duplicates.setup_logging')  # Patch setup_logging
def test_determine_duplicates_with_config_file_loading(
    mock_setup_logging,
    mock_makedirs,
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_sort_table,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_config,
    mock_tree,
    mock_dist_matrix_data,
    mock_contigs_ranked,
    mock_sample_groups,
    mock_classifications,
    mock_seq_records
):
    """Test determine_duplicates with config file loading."""
    # Set up mocks
    mock_root_tree.return_value = mock_tree
    mock_to_distance_matrix.return_value = mock_dist_matrix_data
    mock_sort_table.return_value = mock_contigs_ranked
    mock_group_sequences.return_value = mock_sample_groups
    mock_find_duplicates.return_value = mock_classifications

    # Create a simulated config file path
    config_path = "dummy_config.yaml"

    # Run the function with config path instead of dict
    with patch('builtins.open', mock_open(read_data=yaml.dump(mock_config))):
        with patch('yaml.safe_load', return_value=mock_config):
            with patch('Bio.SeqIO.to_dict', return_value=mock_seq_records):
                with patch('Bio.SeqIO.parse'):
                    result = determine_duplicates(config=config_path)

    # Assertions
    assert result == mock_classifications
    mock_root_tree.assert_called_once()
    mock_find_duplicates.assert_called_once()
    mock_write_results.assert_called_once()


@patch('os.makedirs')
@patch('lasvdedup.utils.determine_duplicates.setup_logging')  # Patch setup_logging
def test_determine_duplicates_missing_segment(mock_setup_logging, mock_makedirs, mock_config):
    """Test error handling when segment is missing."""
    # Remove segment from config
    config_without_segment = mock_config.copy()
    del config_without_segment["segment"]

    # Expect ValueError
    with pytest.raises(ValueError) as excinfo:
        determine_duplicates(config=config_without_segment)

    assert "Segment not provided" in str(excinfo.value)


@patch('os.makedirs')
@patch('lasvdedup.utils.determine_duplicates.setup_logging')  # Patch setup_logging
def test_determine_duplicates_missing_species(mock_setup_logging, mock_makedirs, mock_config):
    """Test error handling when species is missing."""
    # Remove species from config
    config_without_species = mock_config.copy()
    del config_without_species["species"]

    # Expect ValueError
    with pytest.raises(ValueError) as excinfo:
        determine_duplicates(config=config_without_species)

    assert "Species not provided" in str(excinfo.value)


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.sort_table')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
@patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
@patch('os.makedirs')
@patch('lasvdedup.utils.determine_duplicates.setup_logging')  # Patch setup_logging
def test_determine_duplicates_with_explicit_paths(
    mock_setup_logging,
    mock_makedirs,
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_sort_table,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_config,
    mock_tree,
    mock_dist_matrix_data,
    mock_contigs_ranked,
    mock_sample_groups,
    mock_classifications,
    mock_seq_records
):
    """Test determine_duplicates with explicit path parameters."""
    # Set up mocks
    mock_root_tree.return_value = mock_tree
    mock_to_distance_matrix.return_value = mock_dist_matrix_data
    mock_sort_table.return_value = mock_contigs_ranked
    mock_group_sequences.return_value = mock_sample_groups
    mock_find_duplicates.return_value = mock_classifications

    # Explicit paths
    tree_path = Path("explicit_tree.treefile")
    seq_path = Path("explicit_sequences.fasta")
    table_path = Path("explicit_table.tsv")
    prefix_path = Path("explicit_prefix")

    # Run the function with explicit paths
    with patch('Bio.SeqIO.to_dict', return_value=mock_seq_records):
        with patch('Bio.SeqIO.parse'):
            result = determine_duplicates(
                config=mock_config,
                tree=tree_path,
                sequences=seq_path,
                table=table_path,
                prefix=prefix_path,
                segment="L"
            )

    # Assertions
    assert result == mock_classifications
    mock_makedirs.assert_called_once_with(prefix_path, exist_ok=True)
    mock_root_tree.assert_called_once_with(tree_path)
    mock_write_results.assert_called_once()
