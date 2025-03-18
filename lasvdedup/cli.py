#!/usr/bin/env python3
"""Command-line interface for LASV deduplication pipeline."""

import argparse
import sys
import logging
from pathlib import Path

from .utils.determine_duplicates import determine_duplicates
from .utils.config_setup import build_config
from .pipeline import run_pipeline

# Set up logger
logger = logging.getLogger("lasvdedup.cli")

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
    run_parser.add_argument("-l", "--log-level", help="The desired log level (default WARNING).", choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"), default="WARNING")

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
    dedup_parser.add_argument('--sample-regex', '-r', type=str,
                       help='Regular expression to extract sample identifiers')
    dedup_parser.add_argument('--length-column', type=str,
                       help='Name of the column containing consensus length')
    dedup_parser.add_argument('--species', type=str,
                       help='Species name for output files')
    dedup_parser.add_argument('--segment', type=str,
                       help='Segment name for output files')
    dedup_parser.add_argument('--pairwise-distance', type=float,
                       help='Override pairwise distance threshold for intrahost variation')
    dedup_parser.add_argument('--z-threshold', type=float,
                       help='Override z threshold for intrahost variation')
    dedup_parser.add_argument('--clade-size', type=float,
                       help='Override clade size threshold for intrahost variation')
    dedup_parser.add_argument('--clade-size', type=float,
                       help='Override clade size threshold for intrahost variation')
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

def main():
    """Main entry point for the CLI."""
    args = parse_args()

    # Setup basic logging
    logging.basicConfig(
        level=getattr(logging, args.log_level if hasattr(args, 'log_level') else 'WARNING'),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    try:
        if args.command == 'run':
            # Build complete configuration
            config = build_config(args)

            # Run pipeline with complete config
            success = run_pipeline( config=config, dry_run=args.dry_run)
            sys.exit(0 if success else 1)

        elif args.command == 'deduplicate':
            # Create config dict with proper priority order
            config = build_config(args)

            # Run determine_duplicates with config
            determine_duplicates(config=config)
            print(f"Deduplication completed successfully. Results saved to: {args.prefix}")
            sys.exit(0)

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
