import os
import pytest
from pathlib import Path

@pytest.fixture(scope="session", autouse=True)
def create_test_data_dirs():
    """Create the test data directories if they don't exist."""
    # Define directory paths
    data_dir = Path(__file__).parent / "data"
    seq_data_dir = data_dir / "seq_data"
    raw_dir = data_dir / "raw"

    # Create directories
    for dir_path in [data_dir, seq_data_dir, raw_dir]:
        dir_path.mkdir(exist_ok=True)

    # Create minimal test files if they don't exist
    create_test_files(seq_data_dir, raw_dir)

def create_test_files(seq_data_dir, raw_dir):
    """Create minimal test files for testing the workflow."""
    # Create sample sequence files
    for sample_id in ["LVE00001", "LVE00002"]:
        for segment in ["L", "S"]:
            seq_file = seq_data_dir / f"{sample_id}_{segment}.fasta"
            if not seq_file.exists():
                with open(seq_file, "w") as f:
                    f.write(f">{sample_id}_{segment}_contig1\n")
                    f.write("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n")

    # Create base alignment and tree files
    for segment in ["L", "S"]:
        # Base alignment
        aln_file = raw_dir / f"LASV-base-{segment}.aln"
        if not aln_file.exists():
            with open(aln_file, "w") as f:
                f.write(">Reference_LASV1\n")
                f.write("ATGCATGCATGCATGCATGCATGC\n")
                f.write(">Reference_LASV2\n")
                f.write("ATGCATGCATGCATGCATGCATGC\n")

        # Base tree file
        tree_file = raw_dir / f"LASV-base-{segment}.treefile"
        if not tree_file.exists():
            with open(tree_file, "w") as f:
                f.write("(Reference_LASV1:0.01,Reference_LASV2:0.01);\n")
