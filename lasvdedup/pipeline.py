"""Core pipeline functionality for LASV deduplication."""

import os
import logging
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

    # Make a copy of the config to avoid modifying the original
    config_copy = config.copy()

    # Get the current working directory (where the CLI is invoked)
    cwd = os.getcwd()

    # Helper function to resolve paths
    def resolve_path(path_value, is_url_allowed=False):
        path_str = str(path_value)

        # Skip URLs if allowed
        if is_url_allowed and path_str.startswith(("http://", "https://")):
            return path_str

        # Return if already absolute
        if os.path.isabs(path_str):
            return path_str

        # Try relative to CWD first
        cwd_path = os.path.join(cwd, path_str)
        if os.path.exists(cwd_path):
            return cwd_path

        # Default to workdir
        return os.path.join(workdir, path_str)

    # Resolve file and directory paths
    for key in ["CONTIGS_TABLE", "SEQ_DATA_DIR"]:
        if key in config_copy and config_copy[key]:
            config_copy[key] = resolve_path(config_copy[key])

    # Handle BASE_DATA_DIR with URL support
    if "BASE_DATA_DIR" in config_copy and config_copy["BASE_DATA_DIR"]:
        config_copy["BASE_DATA_DIR"] = resolve_path(config_copy["BASE_DATA_DIR"], is_url_allowed=True)

    # Create workdir if it doesn't exist
    os.makedirs(workdir, exist_ok=True)

    # Prepare snakemake arguments
    snakemake_args = {
        'snakefile': str(snakefile),
        'config': config_copy,
        'cores': threads,
        'forceall': force,
        'dryrun': dry_run,
        'printshellcmds': True,
        'workdir': workdir,
    }


    # Run snakemake using the correct API
    success = snakemake.snakemake(**snakemake_args)
    return success

