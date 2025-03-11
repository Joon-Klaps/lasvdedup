#!/usr/bin/env python3

import pytest
import logging
from pathlib import Path

@pytest.fixture(autouse=True)
def set_debug_logging():
    """Set logger to debug level for all tests."""
    # Configure the main logger for the package
    logger = logging.getLogger("lasvdedup")
    logger.setLevel(logging.DEBUG)

    # Make sure we have at least one handler
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    else:
        # Set existing handlers to DEBUG level
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)

    # Also set the duplicates logger specifically
    duplicates_logger = logging.getLogger("lasvdedup.duplicates")
    duplicates_logger.setLevel(logging.DEBUG)

@pytest.fixture
def test_data_dir():
    """Return the path to the test data directory."""
    return Path(__file__).parent / 'data'

@pytest.fixture
def test_tree_path(test_data_dir):
    """Return the path to the test tree file."""
    return test_data_dir / 'test.tree'

@pytest.fixture
def test_fasta_path(test_data_dir):
    """Return the path to the test fasta file."""
    return test_data_dir / 'test.fasta'

@pytest.fixture
def test_tsv_path(test_data_dir):
    """Return the path to the test tsv file."""
    return test_data_dir / 'test.tsv'
