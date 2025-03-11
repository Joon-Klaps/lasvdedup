import os
import pytest
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock

from lasvdedup.utils.io_utils import (
    load_read_counts,
    write_results,
    write_distance_matrix
)
from lasvdedup.utils.classification import Classification, ClassificationType, DecisionCategory

def test_load_read_counts():
    """Test loading read counts from a contigs table."""
    # Create a test DataFrame
    data = {
        'index': ['seq1', 'seq2', 'seq3'],
        'reads': [100, 200, 300],
        'another_reads_column': [50, 150, 250]
    }
    df = pd.DataFrame(data)

    # Save to temporary file
    with tempfile.NamedTemporaryFile(suffix='.tsv', mode='w+', delete=False) as tmp:
        df.to_csv(tmp.name, sep='\t', index=False)
        tmp_path = Path(tmp.name)

    try:
        # Test with direct column name
        result = load_read_counts(tmp_path, 'reads')
        assert len(result) == 3
        assert result['seq1']['reads'] == 100
        assert result['seq2']['reads'] == 200
        assert result['seq3']['reads'] == 300

        # Test with different column name
        result = load_read_counts(tmp_path, 'another_reads_column')
        assert len(result) == 3
        assert result['seq1']['reads'] == 50
        assert result['seq2']['reads'] == 150
        assert result['seq3']['reads'] == 250

        # Test with column not found
        with pytest.raises(ValueError):
            load_read_counts(tmp_path, 'nonexistent_column')

    finally:
        # Clean up
        os.unlink(tmp_path)

def test_write_results():
    """Test writing classification results to file."""
    # Create test classifications
    classifications = {
        'seq1': Classification(
            sequence_name='seq1',
            classification_type=ClassificationType.GOOD,
            reason='test good',
            sample_id='sample1',
            group_members=['seq1', 'seq2'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            read_count=100
        ),
        'seq2': Classification(
            sequence_name='seq2',
            classification_type=ClassificationType.BAD,
            reason='test bad',
            sample_id='sample1',
            group_members=['seq1', 'seq2'],
            decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
            read_count=50
        ),
        'seq3': Classification(
            sequence_name='seq3',
            classification_type=ClassificationType.COINFECTION,
            reason='test coinfection',
            sample_id='sample2',
            group_members=['seq3'],
            decision_category=DecisionCategory.TRUE_COINFECTION,
            read_count=200
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
        with patch('builtins.open', MagicMock()), \
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