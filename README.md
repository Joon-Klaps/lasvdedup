# LASV Deduplication Pipeline

![Python Tests](https://github.com/Joon-Klaps/lasv-deduplicate/actions/workflows/python-tests.yml/badge.svg)
![Conda Test](https://github.com/Joon-Klaps/lasv-deduplicate/actions/workflows/conda-test.yml/badge.svg)

A small snakemake pet project for deduplicating LASV (Lassa virus) sequences, from the output of the [viralgenie pipeline](https://github.com/Joon-Klaps/viralgenie)

## Installation

```bash
# Installation with pip
pip install -e .

# Or using conda
conda env create -f environment.yml
conda activate lasvdedup-env
pip install -e .
```

## Usage

```bash
# Example command
lasvdedup --help
```

## Development

To run tests:

```bash
pytest tests/
```

