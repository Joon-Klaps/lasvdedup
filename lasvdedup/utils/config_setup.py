#!/usr/bin/env python3
"""Set of functions for configuration handling."""

import os
import logging
import yaml
from pathlib import Path

from .resources import get_config_path

# Set up logger
logger = logging.getLogger("lasvdedup.duplicates")

# Define parameter mapping once at module level
# Format: 'cli_arg_name': ('config_path', transform_function or None)
# Where config_path is dot notation for nested dict access (empty string for root level)
PARAMETER_MAPPING = {

    # Run command parameters
    'contigs_table': ('CONTIGS_TABLE', Path),
    'seq_data_dir': ('SEQ_DATA_DIR', Path),
    'base_data_dir': ('BASE_DATA_DIR', Path),
    'outdir': ('OUTDIR', Path),
    'workdir': ('WORKDIR', Path),
    'threads': ('THREADS', None),
    'force': ('FORCE', None),

    # Deduplicate command parameters
    'tree': ('TREE', Path),
    'sequences': ('SEQUENCES', Path),
    'table': ('CONTIGS_TABLE', Path),
    'prefix': ('PREFIX', None),  # Don't convert prefix to abspath
    'sample_regex': ('DEDUPLICATE.SAMPLE_REGEX', None),
    'selection_column': ('DEDUPLICATE.SELECTION_COLUMNS', None),
    'length_column': ('DEDUPLICATE.LENGTH_COLUMN', None),
    'species': ('SPECIES', None),
    'segment': ('segment', None),
    'log_level': ('LOGLEVEL', None),
}

THRESHOLD_MAPPING = {
    'pairwise-distance': 'PWD',
    'quantile': 'QUANTILE',
    'clade-size': 'CLADE_SIZE',
    'target-length': 'TARGET_LENGTH',
}

def build_config(args):
    """
    Build a complete configuration dictionary respecting the priority order:
    1. CLI arguments (highest priority)
    2. Custom config file
    3. Default configuration (lowest priority)

    Uses a parameter mapping to automatically structure the configuration.

    Args:
        args: Parsed command-line arguments

    Returns:
        dict: Complete configuration dictionary
    """
    # Start with default configuration (lowest priority)
    default_config_path = get_config_path()
    config = {}

    # Load default configuration if available
    if os.path.isfile(default_config_path):
        try:
            with open(default_config_path) as f:
                config = yaml.safe_load(f) or {}
            logger.debug("Loaded default configuration from %s", default_config_path)
        except Exception as e:
            logger.warning("Failed to load default config: %s", e)

    # Load custom config if provided (overrides defaults)
    if args.config and isinstance(args.config, dict):
        config = deep_update(config, args.config)

    elif args.config and os.path.isfile(args.config):
        try:
            with open(args.config) as f:
                custom_config = yaml.safe_load(f) or {}
            logger.debug("Loaded custom configuration from %s", args.config)

            # Merge custom config on top of defaults
            config = deep_update(config, custom_config)
        except Exception as e:
            logger.warning("Failed to load custom config from %s: %s", args.config, e)

    # Extract CLI arguments as a dictionary
    cli_config = vars(args).copy()

    # Remove None values and command from CLI config
    cli_config = {k: v for k, v in cli_config.items() if v is not None and k != 'command' and k != 'config'}

    # Process CLI args using parameter mapping
    normalized_cli_config = {}

    # Special handling for threshold parameters when segment is specified
    if 'segment' in cli_config and any([k in THRESHOLD_MAPPING.items() for k in cli_config]):
        segment = cli_config['segment']
        thresholds = {}

        if 'pairwise-distance' in cli_config:
            thresholds['PWD'] = cli_config.pop('pairwise-distance')
        if 'quantile' in cli_config:
            thresholds['QUANTILE'] = cli_config.pop('quantile')
        if 'clade-size' in cli_config:
            thresholds['CLADE_SIZE'] = cli_config.pop('clade-size')
        if 'target-length' in cli_config:
            thresholds['TARGET_LENGTH'] = cli_config.pop('target-length')

        if thresholds:
            normalized_cli_config['DEDUPLICATE']['THRESHOLDS'][segment] = thresholds

    # Process all other arguments using the parameter mapping
    for cli_arg, (config_path, transform_func) in PARAMETER_MAPPING.items():
        if cli_arg in cli_config and cli_arg != 'threshold':
            value = cli_config[cli_arg]

            # Apply transformation function if provided
            if transform_func and callable(transform_func):
                value = transform_func(value)

            # Set the value in the normalized config structure
            if '.' in config_path:
                set_nested_value(normalized_cli_config, config_path, value)
            else:
                normalized_cli_config[config_path] = value

    # Final merge with CLI args taking highest priority
    final_config = deep_update(config, normalized_cli_config)

    logger.debug("Final merged configuration: %s", final_config)
    return final_config

def set_nested_value(config_dict, path, value):
    """
    Set a value in a nested dictionary using dot notation path.

    Args:
        config_dict: Dictionary to update
        path: Dot notation path (e.g., 'DEDUPLICATE.SAMPLE_REGEX')
        value: Value to set

    Returns:
        Updated dictionary
    """
    if not path or '.' not in path:
        # Base case - set directly on the dictionary
        config_dict[path] = value
        return config_dict

    # Split the path into current key and remaining path
    current, remaining = path.split('.', 1)

    # Ensure the current level exists
    if current not in config_dict:
        config_dict[current] = {}
    elif not isinstance(config_dict[current], dict):
        # Convert to dict if it wasn't already
        config_dict[current] = {}

    # Recurse to set the value at the nested level
    set_nested_value(config_dict[current], remaining, value)
    return config_dict

def deep_update(source, updates):
    """
    Recursively update a nested dictionary.
    Args:
        source: Original dictionary to update
        updates: Dictionary with updates to apply

    Returns:
        dict: Updated dictionary
    """
    result = source.copy()

    for key, value in updates.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            # Recursively update nested dictionaries
            result[key] = deep_update(result[key], value)
        else:
            # Replace or add values
            result[key] = value

    return result