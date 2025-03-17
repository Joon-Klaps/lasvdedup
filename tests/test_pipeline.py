import os
import pytest
import shutil
import yaml
import logging
from pathlib import Path

# Import our pipeline runner
from lasvdedup.pipeline import run_pipeline
logger = logging.getLogger("lasvdedup")
logger.setLevel(logging.DEBUG)

@pytest.fixture(scope="session")
def test_workdir(tmp_path_factory):
    """Create a temporary working directory for the tests."""
    workdir = tmp_path_factory.mktemp("snakemake_test")
    logger.info(f"Test working directory: {workdir}")
    return workdir

@pytest.fixture(scope="session")
def test_data_dir():
    """Return the path to the test data directory."""
    return Path(__file__).parent / "data"

@pytest.fixture(scope="session")
def check_dependencies():
    """Check if required external tools are available."""
    deps = ["mafft", "trimal", "iqtree"]
    missing = []

    for tool in deps:
        if shutil.which(tool) is None:
            missing.append(tool)

    if missing:
        pytest.skip(f"Missing required tools: {', '.join(missing)}")

@pytest.fixture(scope="session")
def setup_test_data(test_workdir, test_data_dir):
    """Set up the test data for the workflow."""
    logger.info("Setting up test data...")

    # Create necessary directories
    out_dir = test_workdir / "out"
    os.makedirs(out_dir, exist_ok=True)

    # Copy test data
    seq_data_dir = test_workdir / "seq_data"
    os.makedirs(seq_data_dir, exist_ok=True)

    # Verify source paths exist
    contigs_file = test_data_dir / "contigs-test.tsv"
    if not contigs_file.exists():
        logger.error(f"Test data file not found: {contigs_file}")
        pytest.skip(f"Missing test data file: {contigs_file}")

    # Copy contigs table
    shutil.copy(contigs_file, test_workdir / "contigs-test.tsv")
    logger.info(f"Copied {contigs_file} to {test_workdir / 'contigs-test.tsv'}")

    # Copy sequence files
    seq_data_path = test_data_dir / "seq_data"

    if not seq_data_path.exists():
        logger.error(f"Sequence data directory not found: {seq_data_path}")

    if not any(seq_data_path.rglob("*.fasta")):
        logger.error(f"No sequence data found in {seq_data_path}")

    for seq_file in seq_data_path.rglob("*.fasta"):
        dest_path = seq_data_dir / seq_file.name
        shutil.copy(seq_file, dest_path)
        logger.info(f"Copied {seq_file} to {dest_path}")

    # Set up base data directory
    base_dir = test_workdir / "base_data"
    os.makedirs(base_dir, exist_ok=True)

    for data_file in (test_data_dir / "raw").rglob("*"):
        shutil.copy(data_file, base_dir / data_file.name)

    # Copy test config
    shutil.copy(test_data_dir / "test-config.yaml", test_workdir / "test-config.yaml")

    return test_workdir

@pytest.fixture(scope="session")
def test_config(test_workdir):
    """Load test configuration from file."""
    config_path = test_workdir / "test-config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Set working directory to test directory
    config["WORKDIR"] = str(test_workdir)

    # Set absolute paths to the test data files in temp directory
    config["CONTIGS_TABLE"] = str(test_workdir / "contigs-test.tsv")
    config["SEQ_DATA_DIR"] = str(test_workdir / "seq_data")

    # Make sure BASE_DATA_DIR is correctly set
    # If it's a URL, leave it as is, otherwise set to the local base_data directory
    if not config["BASE_DATA_DIR"].startswith(("http://", "https://")):
        config["BASE_DATA_DIR"] = str(test_workdir / "base_data")

    logger.debug(f"Test configuration: {config}")
    return config

def test_full_workflow(setup_test_data, test_config):
    """Test that the full workflow runs successfully."""
    test_workdir = setup_test_data

    # Run the pipeline using our API
    success = run_pipeline(config=test_config)

    # Check that pipeline completed successfully
    assert success, "Pipeline execution failed"

    # Check that expected output files exist
    output_dir = test_workdir / "test-output"

    assert (output_dir / "dedup/LASV-L-classifications.tsv").exists(), \
        "L segment classifications file not created"
    assert (output_dir / "dedup/LASV-S-classifications.tsv").exists(), \
        "S segment classifications file not created"