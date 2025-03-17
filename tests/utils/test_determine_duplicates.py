import os
import tempfile
import pytest
import yaml
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
from Bio import SeqIO, Phylo
from io import StringIO

from lasvdedup.utils.determine_duplicates import determine_duplicates
from lasvdedup.utils.classification import Classification, ClassificationType, DecisionCategory


# Create test fixtures
@pytest.fixture
def mock_tree_file():
    """Create a temporary mock Newick tree file."""
    tree_content = "(((A:0.1,B:0.2):0.3,C:0.4):0.5,(D:0.6,E:0.7):0.8);"
    with tempfile.NamedTemporaryFile(suffix='.treefile', delete=False) as tmp:
        tmp.write(tree_content.encode('utf-8'))
        tmp_path = tmp.name
    yield tmp_path
    os.unlink(tmp_path)


@pytest.fixture
def mock_sequence_file():
    """Create a temporary mock FASTA file."""
    fasta_content = """>A
ACGTACGTACGT
>B
ACGTACGTACGA
>C
ACGTACGTAGGT
>D
ACGTACGTACCC
>E
ACGCACGTACGT
"""
    with tempfile.NamedTemporaryFile(suffix='.fasta', delete=False) as tmp:
        tmp.write(fasta_content.encode('utf-8'))
        tmp_path = tmp.name
    yield tmp_path
    os.unlink(tmp_path)


@pytest.fixture
def mock_table_file():
    """Create a temporary mock table file."""
    # Fix table format to ensure proper TSV structure
    table_content = "index\tlength\tcoverage\nA\t1000\t20.5\nB\t950\t15.2\nC\t1200\t18.7\nD\t800\t22.3\nE\t1100\t19.8\n"
    with tempfile.NamedTemporaryFile(suffix='.tsv', delete=False) as tmp:
        tmp.write(table_content.encode('utf-8'))
        tmp_path = tmp.name
    yield tmp_path
    os.unlink(tmp_path)


@pytest.fixture
def mock_config():
    """Create a mock config dictionary."""
    return {
        "tree": "tree_path.treefile",
        "sequences": "sequences.fasta",
        "table": "contigs.tsv",
        "prefix": "output_directory",
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


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.sort_table')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
@patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
@patch('os.makedirs')
def test_determine_duplicates_with_mocks(
    mock_makedirs,
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_sort_table,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_config
):
    """Test determine_duplicates with mocked dependencies."""
    # Set up mocks
    mock_tree = {'biophylo': MagicMock(), 'phylodm': MagicMock()}
    mock_root_tree.return_value = mock_tree

    mock_tips = ['A', 'B', 'C', 'D', 'E']
    mock_dist_matrix = np.zeros((5, 5))
    mock_to_distance_matrix.return_value = (mock_tips, mock_dist_matrix)

    mock_contigs_ranked = {
        'A': {'length': 1000, 'coverage': 20.5},
        'B': {'length': 950, 'coverage': 15.2},
        'C': {'length': 1200, 'coverage': 18.7},
        'D': {'length': 800, 'coverage': 22.3},
        'E': {'length': 1100, 'coverage': 19.8}
    }
    mock_sort_table.return_value = mock_contigs_ranked

    mock_samples = {'sample1': ['A', 'B'], 'sample2': ['C', 'D', 'E']}
    mock_group_sequences.return_value = mock_samples

    # Mock classifications with contig_stats instead of read_count
    mock_classifications = {
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
    mock_find_duplicates.return_value = mock_classifications

    # Mock sequence records
    mock_seq_records = {
        'A': MagicMock(),
        'B': MagicMock(),
        'C': MagicMock(),
        'D': MagicMock(),
        'E': MagicMock()
    }

    # Run the function
    with patch('Bio.SeqIO.to_dict', return_value=mock_seq_records):
        result = determine_duplicates(config=mock_config)

    # Assertions
    assert result == mock_classifications
    mock_makedirs.assert_called_once_with(mock_config["prefix"], exist_ok=True)
    mock_root_tree.assert_called_once_with(mock_config["tree"])
    mock_to_distance_matrix.assert_called_once_with(mock_tree['phylodm'])
    mock_sort_table.assert_called_once_with(
        mock_config["table"],
        mock_config["length_column"],
        mock_config["selection_column"],
        expected_length=mock_config["DEDUPLICATE"]["THRESHOLDS"]["L"]["target_length"]
    )
    mock_group_sequences.assert_called_once_with(mock_tips, mock_config["sample_regex"])
    mock_find_duplicates.assert_called_once_with(
        mock_samples,
        mock_tips,
        mock_dist_matrix,
        mock_contigs_ranked,
        mock_tree["biophylo"],
        mock_config["segment"],
        mock_config["DEDUPLICATE"]["THRESHOLDS"]["L"]
    )
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
def test_determine_duplicates_with_config_file_loading(
    mock_makedirs,
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_sort_table,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_config
):
    """Test determine_duplicates with config file loading."""
    # Set up mocks similar to previous test
    mock_tree = {'biophylo': MagicMock(), 'phylodm': MagicMock()}
    mock_root_tree.return_value = mock_tree

    mock_tips = ['A', 'B', 'C', 'D', 'E']
    mock_dist_matrix = np.zeros((5, 5))
    mock_to_distance_matrix.return_value = (mock_tips, mock_dist_matrix)

    mock_contigs_ranked = {
        'A': {'length': 1000, 'coverage': 20.5},
        'B': {'length': 950, 'coverage': 15.2}
    }
    mock_sort_table.return_value = mock_contigs_ranked

    mock_samples = {'sample1': ['A', 'B'], 'sample2': ['C', 'D', 'E']}
    mock_group_sequences.return_value = mock_samples

    mock_classifications = {
        'A': Classification(
            sequence_name='A',
            classification_type=ClassificationType.GOOD,
            reason='Test reason',
            sample_id='sample1',
            group_members=['A', 'B'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            contig_stats={'length': 1000, 'coverage': 20.5}
        )
    }
    mock_find_duplicates.return_value = mock_classifications

    # Mock sequence records
    mock_seq_records = {
        'A': MagicMock(),
        'B': MagicMock(),
        'C': MagicMock(),
        'D': MagicMock(),
        'E': MagicMock()
    }

    # Create a temporary config file
    with tempfile.NamedTemporaryFile(suffix='.yaml', delete=False) as tmp:
        yaml.dump(mock_config, tmp)
        tmp_path = tmp.name

    try:
        # Run the function with config file
        with patch('Bio.SeqIO.to_dict', return_value=mock_seq_records):
            result = determine_duplicates(config=tmp_path)

        # Assertions
        assert result == mock_classifications
        mock_find_duplicates.assert_called_once()
    finally:
        os.unlink(tmp_path)


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.sort_table')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
@patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
@patch('Bio.SeqIO.to_dict')
@patch('Bio.SeqIO.parse')
@patch('os.makedirs')
def test_determine_duplicates_missing_segment(
    mock_makedirs,
    mock_seqio_parse,
    mock_seqio_to_dict,
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_sort_table,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_config
):
    """Test error handling when segment is missing."""
    # Remove segment from config
    config_without_segment = mock_config.copy()
    del config_without_segment["segment"]

    # Expect ValueError
    with pytest.raises(ValueError) as excinfo:
        determine_duplicates(config=config_without_segment)

    assert "Segment not provided" in str(excinfo.value)

    # Verify that none of the processing functions were called
    mock_root_tree.assert_not_called()
    mock_to_distance_matrix.assert_not_called()
    mock_sort_table.assert_not_called()


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.sort_table')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
@patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
@patch('Bio.SeqIO.to_dict')
@patch('Bio.SeqIO.parse')
@patch('os.makedirs')
def test_determine_duplicates_missing_species(
    mock_makedirs,
    mock_seqio_parse,
    mock_seqio_to_dict,
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_sort_table,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_config
):
    """Test error handling when species is missing."""
    # Remove species from config
    config_without_species = mock_config.copy()
    del config_without_species["species"]

    # Expect ValueError
    with pytest.raises(ValueError) as excinfo:
        determine_duplicates(config=config_without_species)

    assert "Species not provided" in str(excinfo.value)
