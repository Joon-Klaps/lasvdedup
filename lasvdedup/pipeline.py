"""Core pipeline functionality for LASV deduplication."""

import os
from pathlib import Path
import snakemake
from .utils.resources import get_snakefile_path

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
    snakefile = get_snakefile_path()

    # Extract some basic parameters from config
    workdir = config.get("WORKDIR", ".")
    threads = config.get("THREADS", 1)
    force = config.get("FORCE", False)

    # Convert workdir to absolute path
    workdir = os.path.abspath(workdir)

    # Resolve relative paths in config for input files
    for key in ["CONTIGS_TABLE", "SEQ_DATA_DIR"]:
        if key in config and isinstance(config[key], (str, Path)) and not os.path.isabs(str(config[key])):
            config[key] = os.path.abspath(str(config[key]))

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
