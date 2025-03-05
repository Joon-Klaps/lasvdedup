"""Core pipeline functionality for LASV deduplication."""

import os
from pathlib import Path
import snakemake


def run_pipeline(config, dry_run=False):
    """
    Run the LASV deduplication pipeline using Snakemake.

    Args:
        config: Complete configuration dictionary with all required parameters
        dry_run: Whether to perform a dry run

    Returns:
        bool: True if the pipeline ran successfully, False otherwise
    """
    # Get Snakefile path from package
    package_dir = Path(__file__).parent
    snakefile = package_dir / "Snakefile"

    # Extract some basic parameters from config
    workdir = config.get("WORKDIR", ".")
    threads = config.get("THREADS", 1)
    force = config.get("FORCE", False)

    # Prepare snakemake arguments
    snakemake_args = {
        'snakefile': str(snakefile),
        'config': config,
        'cores': threads,
        'forceall': force,
        'dryrun': dry_run,
        'printshellcmds': True,
        'workdir': workdir,

    }

    try:
        # Run snakemake using the correct API
        success = snakemake.snakemake(**snakemake_args)
        return success
    except Exception as e:
        print(f"Error running snakemake: {e}")
        return False
