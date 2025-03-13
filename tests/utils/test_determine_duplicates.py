import os
import tempfile
import pytest
import yaml
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
from Bio import SeqIO, Phylo
from io import StringIO

from lasvdedup.utils.determine_duplicates import determine_duplicates, get_segment_thresholds
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
    table_content = "index\treads\nA\t100\nB\t200\nC\t150\nD\t300\nE\t250\n"
    with tempfile.NamedTemporaryFile(suffix='.tsv', delete=False) as tmp:
        tmp.write(table_content.encode('utf-8'))
        tmp_path = tmp.name
    yield tmp_path
    os.unlink(tmp_path)


@pytest.fixture
def mock_config_file():
    """Create a temporary mock config YAML file."""
    # Fix the escape sequence for regex in YAML - use raw string to avoid double escaping
    config_content = r"""
SPECIES: LASV
DEDUPLICATE:
  SAMPLE_REGEX: "(\\w+)_.*"
  READS_COLUMN: "reads"
  DEFAULT_THRESHOLD:
    LOWER: 0.02
    UPPER: 0.05
    CLADE_SIZE: 10
    Z_THRESHOLD: 2.0
  THRESHOLDS:
    L:
      LOWER: 0.01
      UPPER: 0.03
      CLADE_SIZE: 15
      Z_THRESHOLD: 2.5
    S:
      LOWER: 0.03
      UPPER: 0.07
      CLADE_SIZE: 8
      Z_THRESHOLD: 1.8
"""
    with tempfile.NamedTemporaryFile(suffix='.yaml', delete=False) as tmp:
        tmp.write(config_content.encode('utf-8'))
        tmp_path = tmp.name
    yield tmp_path
    os.unlink(tmp_path)


# Test functions
def test_get_segment_thresholds():
    """Test the get_segment_thresholds function."""
    # Test with dictionary config
    config = {
        'DEDUPLICATE': {
            'THRESHOLDS': {
                'L': {
                    'LOWER': 0.01,
                    'UPPER': 0.03,
                    'CLADE_SIZE': 15,
                    'Z_THRESHOLD': 2.5
                }
            },
            'DEFAULT_THRESHOLD': {
                'LOWER': 0.02,
                'UPPER': 0.05,
                'CLADE_SIZE': 10,
                'Z_THRESHOLD': 2.0
            }
        }
    }

    # Test segment-specific thresholds
    lower, upper, clade_size, z_threshold = get_segment_thresholds(config, 'L')
    assert lower == 0.01
    assert upper == 0.03
    assert clade_size == 15
    assert z_threshold == 2.5

    # Test default thresholds for missing segment
    lower, upper, clade_size, z_threshold = get_segment_thresholds(config, 'M')
    assert lower == 0.02
    assert upper == 0.05
    assert clade_size == 10
    assert z_threshold == 2.0

    # Test kwargs fallback
    kwargs = {'lowerthreshold': 0.03, 'upperthreshold': 0.06, 'clade_size': 12, 'z_threshold': 1.5}
    lower, upper, clade_size, z_threshold = get_segment_thresholds(kwargs, 'X')
    assert lower == 0.03
    assert upper == 0.06
    assert clade_size == 12
    assert z_threshold == 1.5


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.load_read_counts')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
@patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
def test_determine_duplicates_with_mocks(
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_load_reads,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_tree_file,
    mock_sequence_file,
    mock_table_file,
    tmp_path
):
    """Test determine_duplicates with mocked dependencies."""
    # Set up mocks
    mock_tree = {'biophylo': MagicMock(), 'phylodm': MagicMock()}
    mock_root_tree.return_value = mock_tree

    mock_tips = ['A', 'B', 'C', 'D', 'E']
    mock_dist_matrix = np.zeros((5, 5))
    mock_to_distance_matrix.return_value = (mock_tips, mock_dist_matrix)

    mock_reads = {'A': {'reads': 100}, 'B': {'reads': 200}}
    mock_load_reads.return_value = mock_reads

    mock_samples = {'sample1': ['A', 'B'], 'sample2': ['C', 'D', 'E']}
    mock_group_sequences.return_value = mock_samples

    # Mock classifications
    mock_classifications = {
        'A': Classification(
            sequence_name='A',
            classification_type=ClassificationType.GOOD,
            reason='Test reason',
            sample_id='sample1',
            group_members=['A', 'B'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            read_count=100
        ),
        'B': Classification(
            sequence_name='B',
            classification_type=ClassificationType.BAD,
            reason='Test reason',
            sample_id='sample1',
            group_members=['A', 'B'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            read_count=200
        )
    }
    mock_find_duplicates.return_value = mock_classifications

    # Create output directory
    output_dir = Path(tmp_path) / "output"
    output_dir.mkdir(exist_ok=True)

    # Run the function
    result = determine_duplicates(
        tree=mock_tree_file,
        sequences=mock_sequence_file,
        table=mock_table_file,
        prefix=str(output_dir),
        species="LASV",
        segment="L",
        sample_regex=r"(\w+)_.*",
        reads_column="reads",
        lowerthreshold=0.01,
        upperthreshold=0.03
    )

    # Assertions
    assert result == mock_classifications
    mock_root_tree.assert_called_once_with(mock_tree_file)
    mock_to_distance_matrix.assert_called_once_with(mock_tree['phylodm'])
    mock_load_reads.assert_called_once()
    mock_group_sequences.assert_called_once_with(mock_tips, r"(\w+)_.*")
    mock_find_duplicates.assert_called_once()
    mock_write_distance_matrix.assert_called_once()
    mock_write_results.assert_called_once()


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.load_read_counts')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
@patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
def test_determine_duplicates_with_config_file(
    mock_write_distance_matrix,
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_load_reads,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_tree_file,
    mock_sequence_file,
    mock_table_file,
    mock_config_file,
    tmp_path
):
    """Test determine_duplicates with config file."""
    # Set up mocks similar to previous test
    mock_tree = {'biophylo': MagicMock(), 'phylodm': MagicMock()}
    mock_root_tree.return_value = mock_tree

    mock_tips = ['A', 'B', 'C', 'D', 'E']
    mock_dist_matrix = np.zeros((5, 5))
    mock_to_distance_matrix.return_value = (mock_tips, mock_dist_matrix)

    mock_reads = {'A': {'reads': 100}, 'B': {'reads': 200}}
    mock_load_reads.return_value = mock_reads

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
            read_count=100
        )
    }
    mock_find_duplicates.return_value = mock_classifications

    # Create output directory
    output_dir = Path(tmp_path) / "output"
    output_dir.mkdir(exist_ok=True)

    # Mock config loading with actual YAML parsing to ensure it works
    with open(mock_config_file, 'r') as f:
        loaded_config = yaml.safe_load(f)

    # Patch yaml.safe_load to return our verified config
    with patch('yaml.safe_load', return_value=loaded_config):
        # Run the function with config file
        determine_duplicates(
            config=mock_config_file,
            tree=mock_tree_file,
            sequences=mock_sequence_file,
            table=mock_table_file,
            prefix=str(output_dir),
            segment="L"  # Should use L-specific thresholds from config
        )

    # Verify L-specific thresholds were used
    mock_find_duplicates.assert_called_once()
    call_args = mock_find_duplicates.call_args[0]
    # The 4th, 5th, and 6th arguments should be contig_to_reads, lowerthreshold, upperthreshold
    assert call_args[3] == mock_reads  # contig_to_reads
    assert call_args[4] == 0.01  # lowerthreshold from config for L segment
    assert call_args[5] == 0.03  # upperthreshold from config for L segment


@patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint')
@patch('lasvdedup.utils.determine_duplicates.to_distance_matrix')
@patch('lasvdedup.utils.determine_duplicates.load_read_counts')
@patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample')
@patch('lasvdedup.utils.determine_duplicates.find_duplicates')
@patch('lasvdedup.utils.determine_duplicates.write_results')
def test_determine_duplicates_with_config_dict(
    mock_write_results,
    mock_find_duplicates,
    mock_group_sequences,
    mock_load_reads,
    mock_to_distance_matrix,
    mock_root_tree,
    mock_tree_file,
    mock_sequence_file,
    mock_table_file,
    tmp_path
):
    """Test determine_duplicates with config dictionary."""
    # Set up mocks
    mock_tree = {'biophylo': MagicMock(), 'phylodm': MagicMock()}
    mock_root_tree.return_value = mock_tree

    mock_tips = ['A', 'B', 'C', 'D', 'E']
    mock_dist_matrix = np.zeros((5, 5))
    mock_to_distance_matrix.return_value = (mock_tips, mock_dist_matrix)

    mock_reads = {'A': {'reads': 100}, 'B': {'reads': 200}}
    mock_load_reads.return_value = mock_reads

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
            read_count=100
        )
    }
    mock_find_duplicates.return_value = mock_classifications

    # Create config dict
    config_dict = {
        'SPECIES': 'LASV',
        'DEDUPLICATE': {
            'SAMPLE_REGEX': r"(\w+)_.*",
            'READS_COLUMN': 'reads',
            'THRESHOLDS': {
                'S': {
                    'LOWER': 0.03,
                    'UPPER': 0.07,
                    'CLADE_SIZE': 8,
                    'Z_THRESHOLD': 1.8
                }
            }
        }
    }

    # Create output directory
    output_dir = Path(tmp_path) / "output"
    output_dir.mkdir(exist_ok=True)

    # Run the function with config dict
    determine_duplicates(
        config=config_dict,
        tree=mock_tree_file,
        sequences=mock_sequence_file,
        table=mock_table_file,
        prefix=str(output_dir),
        segment="S"  # Should use S-specific thresholds from config
    )

    # Verify S-specific thresholds were used
    mock_find_duplicates.assert_called_once()
    call_args = mock_find_duplicates.call_args[0]
    assert call_args[3] == mock_reads  # contig_to_reads
    assert call_args[4] == 0.03  # lowerthreshold from config for S segment
    assert call_args[5] == 0.07  # upperthreshold from config for S segment


def test_integration(mock_tree_file, mock_sequence_file, mock_table_file, tmp_path):
    """Integration test with minimal mocking."""
    with patch('Bio.Phylo.read') as mock_phylo_read, \
         patch('Bio.Phylo.write') as mock_phylo_write, \
         patch('phylodm.PhyloDM.load_from_newick_path') as mock_phylodm_load, \
         patch('lasvdedup.utils.determine_duplicates.load_read_counts') as mock_load_reads, \
         patch('os.makedirs') as mock_makedirs:

        # Mock read counts directly instead of parsing file
        mock_load_reads.return_value = {
            'A': {'reads': 100},
            'B': {'reads': 200},
            'C': {'reads': 150},
            'D': {'reads': 300},
            'E': {'reads': 250}
        }

        # Mock the phylo tree
        mock_tree = MagicMock()
        mock_tree.root_at_midpoint.return_value = None
        mock_tree.get_terminals.return_value = [MagicMock(name=name) for name in ['A', 'B', 'C', 'D', 'E']]
        mock_tree.find_any.side_effect = lambda name: MagicMock(name=name)
        mock_tree.common_ancestor.return_value = mock_tree
        mock_phylo_read.return_value = mock_tree

        # Mock the phylodm
        mock_dm = MagicMock()
        mock_dm.dm.return_value = np.array([
            [0.0, 0.01, 0.02, 0.03, 0.04],
            [0.01, 0.0, 0.05, 0.06, 0.07],
            [0.02, 0.05, 0.0, 0.08, 0.09],
            [0.03, 0.06, 0.08, 0.0, 0.10],
            [0.04, 0.07, 0.09, 0.10, 0.0]
        ])
        mock_dm.taxa.return_value = ['A', 'B', 'C', 'D', 'E']
        mock_phylodm_load.return_value = mock_dm

        # Create output directory
        output_dir = Path(tmp_path) / "output"
        output_dir.mkdir(exist_ok=True)

        # Run the function with minimal arguments
        with patch('builtins.open', mock_open()) as m:
            result = determine_duplicates(
                tree=mock_tree_file,
                sequences=mock_sequence_file,
                table=mock_table_file,
                prefix=str(output_dir),
                species="LASV",
                segment="L",
                sample_regex=r"^([A-E])$",  # Regex that extracts single letter as sample ID
                reads_column="reads",
                lowerthreshold=0.01,
                upperthreshold=0.03
            )

        # Basic assertions
        assert isinstance(result, dict)
        assert len(result) > 0  # Should have some classifications

        # Verify expected function calls
        mock_phylo_read.assert_called_once()
        mock_phylo_write.assert_called_once()
        mock_phylodm_load.assert_called_once()
