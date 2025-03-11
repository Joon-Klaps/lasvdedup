#!/usr/bin/env python3
"""Command-line interface for LASV deduplication pipeline."""

import argparse
import os
import sys
from pathlib import Path
import yaml

from .pipeline import run_pipeline
from .utils.determine_duplicates import determine_duplicates


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='LASV deduplication pipeline')

    # Check if we should use backward compatibility mode
    if len(sys.argv) > 1 and not sys.argv[1] in ['run', 'deduplicate']:
        # No subcommand provided, use 'run' for backward compatibility
        cmd_args = ['run'] + sys.argv[1:]
    else:
        cmd_args = sys.argv[1:]

    subparsers = parser.add_subparsers(dest='command', help='Command to run')

    # Main pipeline command
    run_parser = subparsers.add_parser('run', help='Run the full pipeline')

    # Required arguments for the run command
    run_parser.add_argument('--input', '--contigs-table', dest='contigs_table', required=True,
                      help='Path to contigs table (CSV/TSV format)')

    # Optional arguments for the run command
    run_parser.add_argument('--seq-dir', '--seq-data-dir', dest='seq_data_dir',
                      help='Directory containing sequence data')
    run_parser.add_argument('--ref-dir', '--base-data-dir', dest='base_data_dir',
                      help='Base URL or path for reference data')
    run_parser.add_argument('--outdir', '-o', help='Output directory')
    run_parser.add_argument('--workdir', '-w', help='Working directory')
    run_parser.add_argument('--config', '-c', help='Path to config file')
    run_parser.add_argument('--threads', '-t', type=int, help='Number of threads to use')
    run_parser.add_argument('--force', '-f', action='store_true', help='Force rerun all jobs')
    run_parser.add_argument('--dry-run', '-n', action='store_true', help='Perform a dry run')

    # Deduplicate command
    dedup_parser = subparsers.add_parser('deduplicate', help='Run just the deduplication step')

    # Required arguments for the deduplicate command
    dedup_parser.add_argument('--tree', '-t', required=True, type=Path,
                       help='Path to the phylogenetic tree file')
    dedup_parser.add_argument('--sequences', '-s', required=True, type=Path,
                       help='Path to sequences FASTA file')
    dedup_parser.add_argument('--table', '-i', required=True, type=Path,
                       help='Path to the contigs table')
    dedup_parser.add_argument('--prefix', '-o', required=True, type=str,
                       help='Output directory prefix')

    # Optional arguments for the deduplicate command
    dedup_parser.add_argument('--sample-regex', '-r', type=str, default=r'sample\d+',
                       help='Regular expression to extract sample identifiers')
    dedup_parser.add_argument('--reads-column', type=str, default='reads',
                       help='Name of the column containing read counts')
    dedup_parser.add_argument('--species', type=str, default='LASV',
                       help='Species name for output files')
    dedup_parser.add_argument('--segment', type=str, default='L',
                       help='Segment name for output files')
    dedup_parser.add_argument('--lowerthreshold', type=float,
                       help='Override distance threshold to identify duplicates')
    dedup_parser.add_argument('--upperthreshold', type=float,
                       help='Override distance threshold for intrahost variation')
    dedup_parser.add_argument('--config', '-c', type=str,
                       help='Path to configuration file with segment-specific thresholds')
    dedup_parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )

    # Parse the modified args
    args = parser.parse_args(cmd_args)

    # If no command was provided and we didn't insert 'run', show help
    if not args.command:
        parser.print_help()
        sys.exit(1)

    return args


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
    if not os.path.isfile(contigs_table):
        raise FileNotFoundError(f"Contigs table not found: {contigs_table}")

    # Infer seq_data_dir if not provided
    if not args.seq_data_dir:
        seq_data_dir = os.path.dirname(contigs_table)
        print(f"Sequence data directory not specified, using: {seq_data_dir}")
    else:
        seq_data_dir = os.path.abspath(args.seq_data_dir)

    if not os.path.isdir(seq_data_dir):
        raise FileNotFoundError(f"Sequence data directory not found: {seq_data_dir}")

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

    if args.outdir:
        cli_overrides['OUTDIR'] =  os.path.abspath(args.outdir)

    # Update config with CLI overrides
    config.update(cli_overrides)

    return config


def main():
    """Main entry point for the CLI."""
    args = parse_args()

    try:
        if args.command == 'run':
            # Build complete configuration
            config = build_config(args)

            # Run pipeline with complete config
            success = run_pipeline(
                config=config,
                dry_run=args.dry_run
            )
            sys.exit(0 if success else 1)

        elif args.command == 'deduplicate':
            # Create config dict if no config file provided but thresholds specified
            config_dict = None
            if args.config:
                with open(args.config, 'r') as f:
                    config_dict = yaml.safe_load(f)
            elif args.lowerthreshold is not None or args.upperthreshold is not None:
                # Create minimal config if thresholds were specified directly
                config_dict = {
                    'DEDUPLICATE': {
                        'THRESHOLDS': {
                            args.segment: {
                                'LOWER': args.lowerthreshold,
                                'UPPER': args.upperthreshold
                            }
                        }
                    }
                }

            # Run determine_duplicates with config
            determine_duplicates(
                config=config_dict,
                tree=args.tree,
                sequences=args.sequences,
                prefix=args.prefix,
                table=args.table,
                sample_regex=args.sample_regex,
                reads_column=args.reads_column,
                species=args.species,
                segment=args.segment,
                log_level=args.log_level
            )
            print(f"Deduplication completed successfully. Results saved to: {args.prefix}")
            sys.exit(0)

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
