#!/usr/bin/env python3
"""Command-line interface for LASV deduplication pipeline."""

import argparse
import os
import sys
from pathlib import Path
import yaml

from .pipeline import run_pipeline


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='LASV deduplication pipeline')

    # Required arguments
    parser.add_argument('--input', '--contigs-table', dest='contigs_table', required=True,
                      help='Path to contigs table (CSV/TSV format)')

    # Optional arguments
    parser.add_argument('--seq-dir', '--seq-data-dir', dest='seq_data_dir',
                      help='Directory containing sequence data')
    parser.add_argument('--ref-dir', '--base-data-dir', dest='base_data_dir',
                      help='Base URL or path for reference data')
    parser.add_argument('--workdir', '-w', help='Working directory', default='.')
    parser.add_argument('--config', '-c', help='Path to config file')
    parser.add_argument('--threads', '-t', type=int, help='Number of threads to use')
    parser.add_argument('--force', '-f', action='store_true', help='Force rerun all jobs')
    parser.add_argument('--dry-run', '-n', action='store_true', help='Perform a dry run')

    return parser.parse_args()


def build_config(args):
    """
    Build a complete configuration dictionary from CLI args and config file.

    Args:
        args: Parsed command-line arguments

    Returns:
        dict: Complete configuration dictionary
    """
    # Determine config file path
    package_dir = Path(__file__).parent
    default_config = package_dir / "config.yaml"
    config_path = args.config or default_config

    # Ensure config file exists
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    # Load config from file
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Update paths to absolute paths
    contigs_table = os.path.abspath(args.contigs_table)

    # Infer seq_data_dir if not provided
    if not args.seq_data_dir:
        seq_data_dir = os.path.dirname(contigs_table)
        print(f"Sequence data directory not specified, using: {seq_data_dir}")
    else:
        seq_data_dir = os.path.abspath(args.seq_data_dir)

    # Override config with CLI arguments
    cli_overrides = {
        'CONTIGS_TABLE': contigs_table,
        'SEQ_DATA_DIR': seq_data_dir,
        'WORKDIR': os.path.abspath(args.workdir) if args.workdir else ".",
    }

    if args.base_data_dir:
        cli_overrides['BASE_DATA_DIR'] = args.base_data_dir

    if args.threads:
        cli_overrides['THREADS'] = args.threads

    if args.force:
        cli_overrides['FORCE'] = True

    # Update config with CLI overrides
    config.update(cli_overrides)

    return config


def main():
    """Main entry point for the CLI."""
    args = parse_args()

    try:
        # Build complete configuration
        config = build_config(args)

        # Run pipeline with complete config
        success = run_pipeline(
            config=config,
            dry_run=args.dry_run
        )

        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
