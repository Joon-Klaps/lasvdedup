import os
import pytest
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock

from lasvdedup.utils.io_utils import (
    sort_table,
    write_results,
    write_distance_matrix
)
from lasvdedup.utils.classification import Classification, ClassificationType, DecisionCategory

def test_sort_table():
    """Test sorting contigs table by length and selection criteria."""
    # Create a test DataFrame
    data = {
        'index': ['seq1', 'seq2', 'seq3'],
        'length': [7000, 6800, 7200],
        'coverage': [20.5, 15.2, 18.7]
    }
    df = pd.DataFrame(data)

    # Save to temporary file
    with tempfile.NamedTemporaryFile(suffix='.tsv', mode='w+', delete=False) as tmp:
        df.to_csv(tmp.name, sep='\t', index=False)
        tmp_path = Path(tmp.name)

    try:
        # Test with length column and coverage as selection column
        result = sort_table(tmp_path, 'length', ['coverage'], expected_length=7000)
        assert len(result) == 3

        # seq1 should be ranked first (exact match to expected length and highest coverage)
        assert result['seq1']['rank'] == 1
        assert result['seq1']['distance_to_expectation'] == 0
        assert result['seq1']['coverage'] == 20.5

        # seq3 should be ranked lower due to distance from expected length
        assert 'seq3' in result

        # Test with invalid length column
        with pytest.raises(ValueError):
            sort_table(tmp_path, 'nonexistent_column', ['coverage'], expected_length=7000)

        # Test with invalid selection column
        with pytest.raises(ValueError):
            sort_table(tmp_path, 'length', ['nonexistent_column'], expected_length=7000)

    finally:
        # Clean up
        os.unlink(tmp_path)

def test_write_results():
    """Test writing classification results to file."""
    # Create test classifications with contig_stats instead of read_count
    classifications = {
        'seq1': Classification(
            sequence_name='seq1',
            classification_type=ClassificationType.GOOD,
            reason='test good',
            sample_id='sample1',
            group_members=['seq1', 'seq2'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            contig_stats={'length': 7000, 'coverage': 20.5}
        ),
        'seq2': Classification(
            sequence_name='seq2',
            classification_type=ClassificationType.BAD,
            reason='test bad',
            sample_id='sample1',
            group_members=['seq1', 'seq2'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            contig_stats={'length': 6800, 'coverage': 15.2}
        ),
        'seq3': Classification(
            sequence_name='seq3',
            classification_type=ClassificationType.COINFECTION,
            reason='test coinfection',
            sample_id='sample2',
            group_members=['seq3'],
            decision_category=DecisionCategory.TRUE_COINFECTION,
            contig_stats={'length': 7200, 'coverage': 18.7}
        )
    }

    # Mock sequences
    sequences = {
        'seq1': MagicMock(),
        'seq2': MagicMock(),
        'seq3': MagicMock()
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        # Mock file writing to avoid actual file operations
        with patch('builtins.open', MagicMock()) as mock_open, \
             patch('os.makedirs') as mock_makedirs, \
             patch('Bio.SeqIO.write') as mock_write:

            # Call write_results
            write_results(
                classifications=classifications,
                sequences=sequences,
                species='LASV',
                segment='L',
                prefix=tmpdir,
                sample_regex=r'sample\d+'
            )

            # Check that directories were created
            mock_makedirs.assert_called()

            # Check that sequences were written
            assert mock_write.call_count == 3

            # Verify the paths used for writing sequences
            # Extract paths from mock_open calls
            output_paths = [call[0][0] for call in mock_open.call_args_list if len(call[0]) > 0]

            # Filter only the paths for FASTA files (exclude the classification summary files)
            fasta_paths = [path for path in output_paths if path.endswith('.fasta')]

            # Check that files were written to the correct directories
            assert any('sample1/good' in path for path in fasta_paths)
            assert any('sample1/bad' in path for path in fasta_paths)
            assert any('sample2/good' in path for path in fasta_paths)  # coinfection goes to 'good'

def test_write_distance_matrix():
    """Test writing a distance matrix to file."""
    labels = ['A', 'B', 'C']
    matrix = np.array([
        [0.0, 0.1, 0.2],
        [0.1, 0.0, 0.3],
        [0.2, 0.3, 0.0]
    ])

    with tempfile.NamedTemporaryFile(suffix='.mldist', delete=False) as tmp:
        output_path = tmp.name

    try:
        # Call the function
        write_distance_matrix(matrix, labels, output_path)

        # Verify file contents
        with open(output_path, 'r') as f:
            lines = f.readlines()

        assert lines[0].strip() == '3'  # Should have the number of labels
        assert lines[1].startswith('A\t0.000000\t0.100000\t0.200000')
        assert lines[2].startswith('B\t0.100000\t0.000000\t0.300000')
        assert lines[3].startswith('C\t0.200000\t0.300000\t0.000000')
    finally:
        # Cleanup
        if os.path.exists(output_path):
            os.unlink(output_path)