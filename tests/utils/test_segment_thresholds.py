import pytest
from lasvdedup.utils.determine_duplicates import get_segment_thresholds, determine_duplicates

def test_get_segment_thresholds_from_config():
    """Test extracting segment-specific thresholds from config."""
    # Test with segment-specific thresholds
    config = {
        'DEDUPLICATE': {
            'THRESHOLDS': {
                'L': {'LOWER': 0.03, 'UPPER': 0.05},
                'S': {'LOWER': 0.02, 'UPPER': 0.04}
            }
        }
    }

    # Check L segment thresholds
    lower, upper, clade_size, z_threshold = get_segment_thresholds(config, 'L')
    assert lower == 0.03
    assert upper == 0.05
    assert clade_size == 8  # Should get default value
    assert z_threshold == 2.0  # Should get default value

    # Check S segment thresholds
    lower, upper, clade_size, z_threshold = get_segment_thresholds(config, 'S')
    assert lower == 0.02
    assert upper == 0.04
    assert clade_size == 8  # Should get default value
    assert z_threshold == 2.0  # Should get default value

def test_get_default_thresholds():
    """Test fallback to default thresholds when segment not specified."""
    # Test with default thresholds
    config = {
        'DEDUPLICATE': {
            'THRESHOLDS': {
                'L': {'LOWER': 0.03, 'UPPER': 0.05}
            },
            'DEFAULT_THRESHOLD': {
                'LOWER': 0.02,
                'UPPER': 0.04
            }
        }
    }

    # Check fallback for segment not in config
    lower, upper, clade_size, z_threshold = get_segment_thresholds(config, 'Z')
    assert lower == 0.02
    assert upper == 0.04
    assert clade_size == 8  # Should get default value
    assert z_threshold == 2.0  # Should get default value

def test_backward_compatibility():
    """Test backward compatibility with direct threshold parameters."""
    # Test with direct parameters
    params = {
        'lowerthreshold': 0.01,
        'upperthreshold': 0.03
    }

    # Should use provided thresholds directly
    lower, upper, clade_size, z_threshold = get_segment_thresholds(params, 'L')
    assert lower == 0.01
    assert upper == 0.03
    assert clade_size == 8  # Should get default value
    assert z_threshold == 2.0  # Should get default value

def test_config_integration(mocker):
    """Test integration with determine_duplicates function."""
    # Mock the full determine_duplicates function
    mock_find_duplicates = mocker.patch(
        'lasvdedup.utils.determine_duplicates.find_duplicates',
        return_value={}
    )

    # Mock setup_logging to prevent file creation
    mocker.patch('lasvdedup.utils.determine_duplicates.setup_logging')

    # Create test config with segment-specific thresholds
    config = {
        'DEDUPLICATE': {
            'THRESHOLDS': {
                'L': {'LOWER': 0.03, 'UPPER': 0.05},
                'S': {'LOWER': 0.02, 'UPPER': 0.04}
            }
        }
    }

    # Mock other required functions and parameters
    mocker.patch('lasvdedup.utils.determine_duplicates.root_tree_at_midpoint', return_value={'biophylo': None, 'phylodm': None})
    mocker.patch('lasvdedup.utils.determine_duplicates.to_distance_matrix', return_value=([], []))
    mocker.patch('lasvdedup.utils.determine_duplicates.write_distance_matrix')
    mocker.patch('lasvdedup.utils.determine_duplicates.load_read_counts', return_value={})
    mocker.patch('lasvdedup.utils.determine_duplicates.group_sequences_by_sample', return_value={})

    # Fix the SeqIO.parse issue by mocking the file opening
    mock_seq_parser = mocker.MagicMock()
    mock_to_dict = mocker.patch('Bio.SeqIO.to_dict', return_value={})
    mocker.patch('Bio.SeqIO.parse', return_value=mock_seq_parser)

    mocker.patch('lasvdedup.utils.determine_duplicates.write_results')
    mocker.patch('os.makedirs')

    # Call determine_duplicates with S segment
    determine_duplicates(
        config=config,
        tree="tree.nwk",
        sequences="seqs.fasta",
        prefix="output",
        table="table.tsv",
        segment='S'
    )

    # Verify find_duplicates was called with correct thresholds for S segment
    mock_find_duplicates.assert_called_once()
    args, kwargs = mock_find_duplicates.call_args

    # Check positional arguments - lowerthreshold is at index 4, upperthreshold at index 5
    assert len(args) >= 6, "Not enough positional arguments"
    assert args[4] == 0.02, "Lower threshold not passed correctly as positional argument"
    assert args[5] == 0.04, "Upper threshold not passed correctly as positional argument"

    # Check keyword arguments
    assert kwargs.get('clade_size') == 8
    assert kwargs.get('z_threshold') == 2.0

def test_get_all_segment_thresholds():
    """Test extracting all segment-specific thresholds including clade size and z-threshold."""
    # Test with all possible segment-specific thresholds
    config = {
        'DEDUPLICATE': {
            'THRESHOLDS': {
                'L': {'LOWER': 0.03, 'UPPER': 0.05, 'CLADE_SIZE': 12, 'Z_THRESHOLD': 2.5},
                'S': {'LOWER': 0.02, 'UPPER': 0.04, 'CLADE_SIZE': 8, 'Z_THRESHOLD': 2.0}
            }
        }
    }

    # Check L segment thresholds including clade size and z-threshold
    lower, upper, clade_size, z_threshold = get_segment_thresholds(config, 'L')
    assert lower == 0.03
    assert upper == 0.05
    assert clade_size == 12
    assert z_threshold == 2.5

    # Check S segment thresholds including clade size and z-threshold
    lower, upper, clade_size, z_threshold = get_segment_thresholds(config, 'S')
    assert lower == 0.02
    assert upper == 0.04
    assert clade_size == 8
    assert z_threshold == 2.0
