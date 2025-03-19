import os
import tempfile
import logging
from pathlib import Path
from phylodm import PhyloDM
from Bio import Phylo

logger = logging.getLogger("lasvdedup.duplicates")

def root_tree_at_midpoint(tree_path: Path) -> PhyloDM:
    """Root a phylogenetic tree at midpoint and return a PhyloDM object."""
    logger.info("Loading phylogenetic tree from %s", tree_path)

    tree = {}
    # Read the tree with Bio.Phylo
    bio_tree = Phylo.read(str(tree_path), "newick")

    # Root the tree at midpoint
    logger.info("Rooting tree at midpoint")
    bio_tree.root_at_midpoint()
    tree["BIOPHYLO"] = bio_tree

    # Write the rooted tree to a temporary file
    with tempfile.NamedTemporaryFile(suffix='.treefile', delete=False) as tmp_tree_file:
        Phylo.write(bio_tree, tmp_tree_file.name, "newick")
        rooted_tree_path = tmp_tree_file.name

    # Load the rooted tree with PhyloDM
    logger.info("Loading rooted tree into PhyloDM")
    tree["PHYLODM"] = PhyloDM.load_from_newick_path(rooted_tree_path)

    # Clean up the temporary file
    try:
        os.unlink(rooted_tree_path)
    except:
        logger.warning("Could not remove temporary tree file: %s", rooted_tree_path)

    logger.info("Phylogenetic tree loaded and rooted successfully")
    return tree

def get_mrca_clade(seq_names, tree):
    """Get all members of the clade containing the most recent common ancestor."""
    logger.debug("Finding MRCA clade for %d sequences", len(seq_names))

    # Find the MRCA of the sequences
    mrca = tree.common_ancestor([tree.find_any(name) for name in seq_names])

    # Get all terminals in this clade
    return mrca, [terminal.name for terminal in mrca.get_terminals()]