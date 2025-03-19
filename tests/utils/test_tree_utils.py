import os
import pytest
import tempfile
from pathlib import Path
from io import StringIO
import Bio.Phylo

from lasvdedup.utils.tree_utils import (
    root_tree_at_midpoint,
    get_mrca_clade
)

def test_root_tree_at_midpoint(mocker):
    """Test that root_tree_at_midpoint correctly processes a tree."""
    # Create a simple tree file
    tree_content = "(((A:0.1,B:0.2)Node1:0.3,C:0.4)Node2:0.5,D:0.6);"

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.nwk', delete=False) as tmp:
        tmp.write(tree_content)
        tree_path = Path(tmp.name)

    try:
        # Mock the intermediate functions to isolate the test
        mock_phylo_read = mocker.patch('Bio.Phylo.read', autospec=True)
        mock_tree = mock_phylo_read.return_value
        mock_phylo_write = mocker.patch('Bio.Phylo.write')
        mock_phylodm = mocker.patch('phylodm.PhyloDM.load_from_newick_path')
        mock_os_unlink = mocker.patch('os.unlink')

        # Call the function
        result = root_tree_at_midpoint(tree_path)

        # Verify the function called the expected dependencies
        mock_phylo_read.assert_called_once_with(str(tree_path), "newick")
        mock_tree.root_at_midpoint.assert_called_once()
        mock_phylo_write.assert_called_once()
        mock_phylodm.assert_called_once()
        mock_os_unlink.assert_called_once()
    finally:
        # Cleanup
        if os.path.exists(tree_path):
            os.unlink(tree_path)

def test_get_mrca_clade():
    """Test getting MRCA clade members."""
    # Create a simple tree for testing
    tree_str = '(((A:0.1,B:0.2)Node1:0.3,(C:0.2,E:0.3)Node3:0.2)Node2:0.5,D:0.6);'
    tree = Bio.Phylo.read(StringIO(tree_str), 'newick')

    # Test with two sequences
    _, clade_members = get_mrca_clade(['A', 'B'], tree)
    assert sorted(clade_members) == ['A', 'B']

    # Test with sequences in different clades
    _, clade_members = get_mrca_clade(['A', 'C'], tree)
    assert sorted(clade_members) == ['A', 'B', 'C', 'E']

    # Test with all sequences
    _,clade_members = get_mrca_clade(['A', 'B', 'C', 'D', 'E'], tree)
    assert sorted(clade_members) == ['A', 'B', 'C', 'D', 'E']
