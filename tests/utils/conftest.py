import pytest
import numpy as np

@pytest.fixture
def mock_phylodm():
    """Create a mock PhyloDM object."""
    class MockPhyloDM:
        def __init__(self, taxa, distances):
            self._taxa = taxa
            self._distances = distances

        def dm(self, norm=False):
            return self._distances

        def taxa(self):
            return self._taxa

    return MockPhyloDM
