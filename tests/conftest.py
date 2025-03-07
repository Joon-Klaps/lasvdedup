#!/usr/bin/env python3

import pytest
import os
import tempfile
import shutil
from pathlib import Path

@pytest.fixture
def test_data_dir():
    """Create a temporary directory for test data."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)
