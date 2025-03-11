import os
import pytest
import shutil
import subprocess
from pathlib import Path

@pytest.fixture(scope="session")
def test_workdir(tmp_path_factory):
    """Create a temporary working directory for the tests."""
    workdir = tmp_path_factory.mktemp("snakemake_test")
    return workdir

@pytest.fixture(scope="session")
def test_data_dir():
    """Return the path to the test data directory."""
    return Path(__file__).parent / "data"

@pytest.fixture(scope="session")
def setup_test_data(test_workdir, test_data_dir):
    """Set up the test data for the workflow."""
    # Create necessary directories
    out_dir = test_workdir / "out"
    os.makedirs(out_dir, exist_ok=True)

    # Copy test data
    seq_data_dir = test_workdir / "seq_data"
    os.makedirs(seq_data_dir, exist_ok=True)

    # Copy contigs table
    shutil.copy(test_data_dir / "contigs-test.tsv", test_workdir / "contigs-test.tsv")

    # Copy sequence files
    for seq_file in (test_data_dir / "seq_data").glob("*.fasta"):
        shutil.copy(seq_file, seq_data_dir / seq_file.name)

    # Set up base data directory
    base_dir = test_workdir / "base_data"
    os.makedirs(base_dir, exist_ok=True)

    for data_file in (test_data_dir / "raw").glob("*"):
        shutil.copy(data_file, base_dir / data_file.name)

    # Copy test config
    shutil.copy(test_data_dir / "test-config.yaml", test_workdir / "test-config.yaml")

    return test_workdir

def test_full_workflow(setup_test_data):
    """Test that the full workflow runs successfully."""
    test_workdir = setup_test_data

    # Path to the Snakefile
    snakefile_path = Path(__file__).parents[2] / "lasvdedup" / "Snakefile"

    # Run snakemake
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile", str(snakefile_path),
            "--configfile", str(test_workdir / "test-config.yaml"),
            "--cores", "2",
            "--directory", str(test_workdir),
            "--printshellcmds",
        ],
        capture_output=True,
        text=True,
    )

    # Print output for debugging if the command failed
    if result.returncode != 0:
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")

    # Check that snakemake completed successfully
    assert result.returncode == 0, "Snakemake workflow failed"

    # Check that expected output files exist
    output_dir = test_workdir / "test-output"

    assert (output_dir / "dedup/LASV-L-classifications.tsv").exists(), \
        "L segment classifications file not created"
    assert (output_dir / "dedup/LASV-S-classifications.tsv").exists(), \
        "S segment classifications file not created"

def test_individual_rules(setup_test_data):
    """Test individual rules from the workflow."""
    test_workdir = setup_test_data
    snakefile_path = Path(__file__).parents[2] / "lasvdedup" / "Snakefile"

    # Test rules to run individually
    rules = [
        "prepare_base_files",
        "extract_sequences",
        "align_sequences",
        "trim_alignment"
    ]

    for rule in rules:
        result = subprocess.run(
            [
                "snakemake",
                "--snakefile", str(snakefile_path),
                "--configfile", str(test_workdir / "test-config.yaml"),
                "--directory", str(test_workdir),
                "--cores", "1",
                rule,
                "-F",  # Force run the rule
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(f"RULE {rule} FAILED")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")

        assert result.returncode == 0, f"Rule {rule} failed to execute"
