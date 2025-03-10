#!/usr/bin/env python3
## I have little clue to what is going on here, 99% was written by claude 3.7 sonnet preview thinking.

import os
import pytest
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from phylodm import PhyloDM
from Bio import SeqIO
from io import StringIO

import lasvdedup.utils.determine_duplicates as dedup_module

from lasvdedup.utils.determine_duplicates import (
    to_distance_matrix,
    write_distance_matrix,
    load_read_counts,
    group_sequences_by_sample,
    find_duplicates,
    get_distances,
    get_reads,
    write_results,
    determine_duplicates,
    root_tree_at_midpoint  # Add import for new function
)

# Fixtures for test data
@pytest.fixture
def example_tree():
    """Create a simple newick tree for testing."""
    tree_str = "(((A:0.1,B:0.2)Node1:0.3,C:0.4)Node2:0.5,D:0.6);"
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".tree", delete=False) as f:
        f.write(tree_str)
        filepath = f.name
    try:
        # Load and return the tree
        yield PhyloDM.load_from_newick_path(filepath)
    finally:
        # Clean up the temporary file
        if os.path.exists(filepath):
            os.unlink(filepath)


@pytest.fixture
def example_terminal_tips(example_tree):
    """Return terminal tips from the example tree."""
    return [l.replace(" ", "_") for l in example_tree.taxa()]


@pytest.fixture
def example_distance_matrix(example_tree):
    """Create a distance matrix from the example tree."""
    tips, matrix = to_distance_matrix(example_tree)
    return tips, matrix


@pytest.fixture
def example_read_counts_df():
    """Create a sample contigs table with read counts."""
    data = {
        "index": ["sampleA_1", "sampleA_2", "sampleB_1", "sampleB_2"],
        "reads": [100, 200, 150, 50],
    }
    return pd.DataFrame(data)


@pytest.fixture
def example_read_counts_file(example_read_counts_df):
    """Create a temporary file with read counts data."""
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".tsv", delete=False) as f:
        example_read_counts_df.to_csv(f.name, sep="\t", index=False)
        filename = f.name

    yield Path(filename)

    # Cleanup
    os.unlink(filename)


@pytest.fixture
def sample_to_seqs_mapping():
    """Create a sample-to-sequence mapping."""
    return {
        "sampleA": ["sampleA_1", "sampleA_2"],
        "sampleB": ["sampleB_1"],
    }


@pytest.fixture
def contig_to_reads_mapping():
    """Create a contig-to-reads mapping."""
    return {
        "sampleA_1": {"reads": 100},
        "sampleA_2": {"reads": 200},
        "sampleB_1": {"reads": 150},
        "sampleB_2": {"reads": 50},
    }


@pytest.fixture
def example_sequences():
    """Create sample SeqIO records for testing."""
    seq_records = {
        'sampleA_1': SeqIO.SeqRecord(seq='ACGT', id='sampleA_1', description=''),
        'sampleA_2': SeqIO.SeqRecord(seq='ACTT', id='sampleA_2', description=''),
        'sampleB_1': SeqIO.SeqRecord(seq='GCTA', id='sampleB_1', description=''),
        'sampleB_2': SeqIO.SeqRecord(seq='GCAT', id='sampleB_2', description=''),
    }
    return seq_records


# Tests for each function
def test_to_distance_matrix(example_tree):
    """Test that to_distance_matrix creates correct matrix structure."""
    tips, matrix = to_distance_matrix(example_tree)

    # Check structure
    assert len(tips) > 0
    assert matrix.shape == (len(tips), len(tips))

    # Check diagonal elements (distance to self)
    for i in range(len(tips)):
        assert matrix[i, i] == 0  # Distance to self is 0

    # Find indices of specific terminals - directly use the taxon names
    a_idx = tips.index('A')
    b_idx = tips.index('B')
    c_idx = tips.index('C')
    d_idx = tips.index('D')

    # Check specific pairwise distances based on the tree structure:
    # (((A:0.1,B:0.2)Node1:0.3,C:0.4)Node2:0.5,D:0.6);

    # Exact or approximate distance checks, depending on how the matrix is calculated
    assert 0.3 <= float(matrix[a_idx, b_idx]) <= 0.3001  # A to B: 0.1 + 0.2 = 0.3
    assert 0.8 <= float(matrix[a_idx, c_idx]) <= 0.8001  # A to C: 0.1 + 0.3 + 0.4 = 0.8
    assert 1.5 <= float(matrix[a_idx, d_idx]) <= 1.5001  # A to D: 0.1 + 0.3 + 0.5 + 0.6 = 1.5
    assert 0.9 <= float(matrix[b_idx, c_idx]) <= 0.9001  # B to C: 0.2 + 0.3 + 0.4 = 0.9
    assert 1.6 <= float(matrix[b_idx, d_idx]) <= 1.6001  # B to D: 0.2 + 0.3 + 0.5 + 0.6 = 1.6
    assert 1.5 <= float(matrix[c_idx, d_idx]) <= 1.5001  # C to D: 0.4 + 0.5 + 0.6 = 1.5

    # Verify symmetry of distance matrix
    assert float(matrix[a_idx, b_idx]) == float(matrix[b_idx, a_idx])
    assert float(matrix[a_idx, c_idx]) == float(matrix[c_idx, a_idx])
    assert float(matrix[a_idx, d_idx]) == float(matrix[d_idx, a_idx])


def test_write_distance_matrix(example_terminal_tips, example_distance_matrix):
    """Test that distance matrix is written correctly to file."""
    _, matrix = example_distance_matrix

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
        filename = f.name

    try:
        write_distance_matrix(matrix, example_terminal_tips, filename)

        # Verify file content
        with open(filename, "r") as f:
            lines = f.readlines()

            # Check header
            assert len(lines) > 0
            assert lines[0].strip() == str(len(example_terminal_tips))

            # Check that there's a line for each terminal
            assert len(lines) == len(example_terminal_tips) + 1

    finally:
        os.unlink(filename)


def test_load_read_counts(example_read_counts_file):
    """Test loading read counts from a contigs table."""
    # Test with standard column name
    contig_to_reads = load_read_counts(example_read_counts_file, "reads")

    assert len(contig_to_reads) == 4
    assert contig_to_reads["sampleA_1"]["reads"] == 100
    assert contig_to_reads["sampleA_2"]["reads"] == 200

    # Test error when column not found
    with pytest.raises(ValueError):
        load_read_counts(example_read_counts_file, "nonexistent_column")


def test_group_sequences_by_sample():
    """Test grouping sequences by sample ID using regex."""
    mock_tips = [
        'sampleA_1_contig',
        'sampleA_2_contig',
        'sampleB_1_contig',
        'sampleC_xyz',
    ]

    # Test with sample regex
    sample_regex = r'sample[A-C]_\d+'
    sample_to_seqs = group_sequences_by_sample(mock_tips, sample_regex)

    assert len(sample_to_seqs) == 3  # sampleA, sampleB, sampleC
    assert len(sample_to_seqs['sampleA_1']) == 1
    assert len(sample_to_seqs['sampleA_2']) == 1

    # Test with different regex
    sample_regex = r'sample[A-C]'
    sample_to_seqs = group_sequences_by_sample(mock_tips, sample_regex)

    assert len(sample_to_seqs) == 3
    assert len(sample_to_seqs['sampleA']) == 2  # Both sampleA sequences


def test_get_distances():
    """Test getting distances between sequences."""
    names = ["seq1", "seq2", "seq3"]
    tips_lookup = {"seq1": 0, "seq2": 1, "seq3": 2}
    matrix = np.array([
        [0, 0.1, 0.2],
        [0.1, 0, 0.3],
        [0.2, 0.3, 0]
    ])

    distances = get_distances(names, tips_lookup, matrix)

    # Should return 3 distances for pairs: (seq1,seq2), (seq1,seq3), (seq2,seq3)
    assert len(distances) == 3
    assert 0.1 in distances  # seq1-seq2 distance
    assert 0.2 in distances  # seq1-seq3 distance
    assert 0.3 in distances  # seq2-seq3 distance


def test_get_reads():
    """Test retrieving read counts for sequences."""
    contig_to_reads = {
        "sampleA_1": {"reads": 100},
        "sampleB_1": {"reads": 150},
    }

    # Test direct matching
    assert get_reads("sampleA_1", contig_to_reads) == 100

    # Test with _R_ prefix matching
    assert get_reads("_R_sampleB_1", contig_to_reads) == 150

    # Test error on missing sequence
    with pytest.raises(ValueError):
        get_reads("nonexistent", contig_to_reads)


def test_find_duplicates():
    """Test duplicate identification based on thresholds."""
    sample_to_seqs = {
        "sampleA": ["sampleA_1", "sampleA_2"],  # Close sequences
        "sampleB": ["sampleB_1", "sampleB_2"],  # Distant sequences
    }

    tips = [
        'sampleA_1',
        'sampleA_2',
        'sampleB_1',
        'sampleB_2',
    ]

    tips_lookup = {t: i for i, t in enumerate(tips)}

    # Distance matrix: close for sampleA (0.01), distant for sampleB (0.5)
    dist_matrix = np.array([
        [0, 0.01, 0.3, 0.3],
        [0.01, 0, 0.3, 0.3],
        [0.3, 0.3, 0, 0.5],
        [0.3, 0.3, 0.5, 0]
    ])

    contig_to_reads = {
        "sampleA_1": {"reads": 100},
        "sampleA_2": {"reads": 200},  # Higher read count
        "sampleB_1": {"reads": 150},
        "sampleB_2": {"reads": 50},
    }

    threshold = 0.02  # Threshold to determine duplicates

    classifications = find_duplicates(
        sample_to_seqs, tips, dist_matrix, contig_to_reads, threshold
    )

    # Check classifications
    assert classifications['sampleA_1'] == 'bad'
    assert classifications['sampleA_2'] == 'good'  # Highest read count
    assert classifications['sampleB_1'] == 'coinfection'  # Distant sequences
    assert classifications['sampleB_2'] == 'coinfection'  # Distant sequences


def test_write_results(example_sequences):
    """Test writing results to organized directories and files."""
    classifications = {
        'sampleA_1': 'good',
        'sampleA_2': 'bad',
        'sampleB_1': 'coinfection',
    }

    sample_regex = r'sample[A-C]'
    species = 'virus'
    segment = 'L'

    with tempfile.TemporaryDirectory() as tmpdir:
        # Call write_results
        write_results(classifications, example_sequences, species, segment, tmpdir, sample_regex)

        # Check classifications summary file
        parent_dir = os.path.dirname(tmpdir)
        (f"{parent_dir}/{species}-{segment}-classifications.txt", "w")
        with open(f"{parent_dir}/{species}-{segment}-classifications.tsv", "r") as f:
            lines = f.readlines()
            assert len(lines) == 4  # Header + 3 sequences
            assert "tip name\tclassification\n" == lines[0]
            assert "sampleA_1\tgood\n" in lines
            assert "sampleA_2\tbad\n" in lines
            assert "sampleB_1\tcoinfection\n" in lines

        # Check directory structure and files
        assert os.path.exists(f"{tmpdir}/sampleA/good/sampleA_1_{species}_{segment}.fasta")
        assert os.path.exists(f"{tmpdir}/sampleA/bad/sampleA_2_{species}_{segment}.fasta")
        assert os.path.exists(f"{tmpdir}/sampleB/good/sampleB_1_{species}_{segment}.fasta")  # coinfection goes in good dir
        assert os.path.exists(f"{tmpdir}/sampleB/good/sampleB_1_{species}_{segment}.fasta")  # coinfection goes in good dir

        # Verify content of a FASTA file
        with open(f"{tmpdir}/sampleA/good/sampleA_1_{species}_{segment}.fasta", "r") as f:
            content = f.read()
            assert "sampleA_1" in content
            assert "ACGT" in content


def test_root_tree_at_midpoint(mocker):
    """Test that a tree is correctly rooted at midpoint."""
    # Create a test tree file
    unrooted_tree_content = "(((A:0.1,B:0.2)Node1:0.3,C:0.4)Node2:0.5,D:0.6);"

    with tempfile.NamedTemporaryFile(mode="w+", suffix=".nwk", delete=False) as f:
        f.write(unrooted_tree_content)
        tree_path = Path(f.name)

    try:
        # Set up mocks
        mock_tree = mocker.Mock()
        mock_phylo_read = mocker.patch("Bio.Phylo.read", return_value=mock_tree)
        mock_phylo_write = mocker.patch("Bio.Phylo.write")

        mock_phylodm_tree = mocker.Mock()
        mock_phylodm_load = mocker.patch(
            "phylodm.PhyloDM.load_from_newick_path",
            return_value=mock_phylodm_tree
        )

        # Execute the function
        result = root_tree_at_midpoint(tree_path)

        # Verify the function calls
        mock_phylo_read.assert_called_once_with(str(tree_path), "newick")
        mock_tree.root_at_midpoint.assert_called_once()
        mock_phylo_write.assert_called_once()
        mock_phylodm_load.assert_called_once()

        # Check result
        assert result == mock_phylodm_tree

    finally:
        # Clean up
        if os.path.exists(tree_path):
            os.unlink(tree_path)


def test_determine_duplicates(mocker):
    """Test the main determine_duplicates function."""
    # Mock tree path
    mock_tree_path = Path("tree.nwk")
    mock_sequences = Path("sequences.fasta")
    mock_table = Path("contigs.tsv")
    species = 'virus'
    segment = 'L'

    # Mock the various functions that determine_duplicates calls
    mocker.patch("pathlib.Path.exists", return_value=True)
    mock_makedirs = mocker.patch("os.makedirs")

    # Mock the root_tree_at_midpoint function instead of Phylo.read
    mock_tree_obj = mocker.Mock()
    mock_root_tree = mocker.patch.object(
        dedup_module, "root_tree_at_midpoint",
        return_value=mock_tree_obj
    )

    # Create mock sequence records
    mock_seq_records = {
        'seq1': mocker.Mock(),
        'seq2': mocker.Mock()
    }

    # Fix the sequence parsing mocks to properly handle the iterator pattern
    mock_seq_iter = list(mock_seq_records.keys())  # Use a list instead of iterator
    mock_seqio = mocker.patch("Bio.SeqIO.parse", return_value=mock_seq_iter)
    mocker.patch("Bio.SeqIO.to_dict", return_value=mock_seq_records)

    # Patch the module functions directly
    mock_to_distance_matrix = mocker.patch.object(
        dedup_module, "to_distance_matrix",
        return_value=(["tip1", "tip2"], np.array([[0, 0.1], [0.1, 0]]))
    )
    mock_write_distance_matrix = mocker.patch.object(
        dedup_module, "write_distance_matrix"
    )
    mock_load_read_counts = mocker.patch.object(
        dedup_module, "load_read_counts",
        return_value={"tip1": {"reads": 100}}
    )
    mock_group_sequences = mocker.patch.object(
        dedup_module, "group_sequences_by_sample",
        return_value={"sample1": ["tip1", "tip2"]}
    )
    mock_find_duplicates = mocker.patch.object(
        dedup_module, "find_duplicates",
        return_value={"tip1": "good", "tip2": "bad"}
    )
    mock_write_results = mocker.patch.object(
        dedup_module, "write_results"
    )

    # Call the function
    determine_duplicates(
        tree=mock_tree_path,
        sequences=mock_sequences,
        prefix="output",
        table=mock_table,
        sample_regex=r"sample\d+",
        reads_column="reads",
        species=species,
        segment=segment,
        threshold=0.02,
        log_level="DEBUG"
    )

    # Verify that the appropriate functions were called
    mock_makedirs.assert_called()
    mock_root_tree.assert_called_once_with(mock_tree_path)
    mock_to_distance_matrix.assert_called_once_with(mock_tree_obj)
    mock_write_distance_matrix.assert_called_once()
    mock_load_read_counts.assert_called_once()
    mock_group_sequences.assert_called_once()
    mock_find_duplicates.assert_called_once()

    # Check that write_results was called with the correct parameters
    mock_write_results.assert_called_once_with(
        {"tip1": "good", "tip2": "bad"},  # classifications
        mock_seq_records,  # sequences
        species,  # species
        segment,  # segment
        "output",  # prefix
        r"sample\d+"  # sample_regex
    )
