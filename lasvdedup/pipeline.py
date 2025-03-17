"""Core pipeline functionality for LASV deduplication."""

import os
import logging
from pathlib import Path
import snakemake
from .utils.resources import get_snakefile_path

# Set up logger
logger = logging.getLogger("lasvdedup.pipeline")

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

    # Log the working directory being used
    logger.debug(f"Using working directory: {workdir}")

    # Make a copy of the config to avoid modifying the original
    config_copy = config.copy()

    # Get the current working directory (where the CLI is invoked)
    cwd = os.getcwd()
    logger.debug(f"Current working directory: {cwd}")

    # Resolve paths but preserve URLs - Yeah idk what's going on here anymore
    for key in ["CONTIGS_TABLE", "SEQ_DATA_DIR"]:
        if key in config_copy and config_copy[key]:
            path_value = str(config_copy[key])

            # Skip if already an absolute path
            if os.path.isabs(path_value):
                logger.debug(f"Using absolute {key}: {path_value}")
            else:
                # Check if file exists relative to current working directory
                cwd_path = os.path.join(cwd, path_value)
                if os.path.exists(cwd_path):
                    config_copy[key] = cwd_path
                    logger.debug(f"Found {key} in current directory: {cwd_path}")
                else:
                    # If not in current directory, try workdir
                    work_path = os.path.join(workdir, path_value)
                    config_copy[key] = work_path
                    logger.debug(f"Using {key} relative to workdir: {work_path}")

                    # Check if path exists and warn if not
                    if not os.path.exists(work_path):
                        logger.warning(f"Path for {key} does not exist: {work_path}")

    # Special handling for BASE_DATA_DIR to preserve URLs
    if "BASE_DATA_DIR" in config_copy and config_copy["BASE_DATA_DIR"]:
        base_dir = str(config_copy["BASE_DATA_DIR"])
        # Only convert to absolute path if it's not a URL
        if not base_dir.startswith(("http://", "https://")):
            if not os.path.isabs(base_dir):
                # Check if directory exists relative to current working directory
                cwd_path = os.path.join(cwd, base_dir)
                if os.path.exists(cwd_path):
                    config_copy["BASE_DATA_DIR"] = cwd_path
                    logger.debug(f"Found BASE_DATA_DIR in current directory: {cwd_path}")
                else:
                    # If not in current directory, use workdir
                    config_copy["BASE_DATA_DIR"] = os.path.join(workdir, base_dir)
                    logger.debug(f"Using BASE_DATA_DIR relative to workdir: {config_copy['BASE_DATA_DIR']}")

        logger.debug(f"Using BASE_DATA_DIR: {config_copy['BASE_DATA_DIR']}")

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

    try:
        # Run snakemake using the correct API
        logger.debug(f"Running snakemake with config: {config_copy}")
        success = snakemake.snakemake(**snakemake_args)
        return success
    except Exception as e:
        logger.error(f"Error running snakemake: {e}")
        return False
