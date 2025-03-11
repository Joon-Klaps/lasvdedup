# Snakemake Workflow Tests

This directory contains tests for the LASV deduplication Snakemake workflow.

## Running Tests

To run the tests:

```bash
# From the repository root directory
pytest tests/snakemake -v
```

## Test Structure

- `test_workflow.py`: Main test script that runs the Snakemake workflow
- `conftest.py`: Pytest configuration and fixtures
- `data/`: Contains test data files
  - `test-config.yaml`: Test configuration
  - `contigs-test.tsv`: Minimal contigs table for testing
  - `seq_data/`: Test sequence files
  - `raw/`: Base alignments and tree files

## Troubleshooting

If tests fail, check the test output to identify issues. The most common causes of failure:

1. Missing dependencies (snakemake, mafft, iqtree, etc.)
2. Path issues (the tests expect the Snakefile to be in a specific location)
3. Permissions issues when creating temporary files
