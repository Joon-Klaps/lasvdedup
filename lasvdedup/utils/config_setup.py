#!/usr/bin/env python3
"""Set of functions for configuration handling."""

import os
import logging
import yaml
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Type, Union

from .resources import get_config_path

# Set up logger
logger = logging.getLogger("lasvdedup.duplicates")

# Define parameter mapping once at module level
# Format: 'cli_arg_name': ('config_path', transform_function)
# Where config_path is dot notation for nested dict access (empty string for root level)
PARAMETER_MAPPING = {
    # Run command parameters
    'contigs_table': ('CONTIGS_TABLE', Path),
    'seq_data_dir': ('SEQ_DATA_DIR', Path),
    'base_data_dir': ('BASE_DATA_DIR', Path),
    'outdir': ('OUTDIR', Path),
    'workdir': ('WORKDIR', Path),
    'threads': ('THREADS', int),
    'force': ('FORCE', bool),
    'dry_run': ('DRY_RUN', bool),

    # Deduplicate command parameters
    'tree': ('TREE', Path),
    'sequences': ('SEQUENCES', Path),
    'table': ('CONTIGS_TABLE', Path),
    'prefix': ('PREFIX', str),
    'sample_regex': ('DEDUPLICATE.SAMPLE_REGEX', str),
    'selection_column': ('DEDUPLICATE.SELECTION_COLUMNS', lambda x: x.split(',') if isinstance(x, str) else x),
    'length_column': ('DEDUPLICATE.LENGTH_COLUMN', str),
    'species': ('SPECIES', str),
    'segment': ('segment', str),
    'log_level': ('LOGLEVEL', str),

    # Threshold parameters
    'pairwise_distance': ('pairwise_distance', float),
    'z_threshold': ('z_threshold', float),
    'clade_size': ('clade_size', int),
    'target_length': ('target_length', int),
}

THRESHOLD_MAPPING = {
    'pairwise_distance': 'PWD',
    'z_threshold': 'Z_THRESHOLD',
    'clade_size': 'CLADE_SIZE',
    'target_length': 'TARGET_LENGTH',
}

def validate_and_cast(value: Any, transform_func: Callable) -> Any:
    """
    Validate and cast a value using the provided transform function.

    Args:
        value: The value to validate and cast
        transform_func: Function to use for transformation

    Returns:
        The transformed value, or the original value if transformation failed
    """
    if value is None:
        return None

    # Check if value is already the right type
    if transform_func is Path and isinstance(value, Path):
        return value
    elif transform_func is bool and isinstance(value, bool):
        return value
    elif transform_func is int and isinstance(value, int):
        return value
    elif transform_func is float and isinstance(value, float):
        return value
    elif transform_func is str and isinstance(value, str):
        return value
    elif callable(transform_func) and not isinstance(transform_func, type):
        # For custom functions like lambdas
        try:
            return transform_func(value)
        except Exception as e:
            logger.warning(f"Failed to transform {value} with custom function: {e}")
            return value

    try:
        # Handle special case for bool since bool("False") is True
        if transform_func is bool and isinstance(value, str):
            if value.lower() in ("false", "0", "no", "n", "off"):
                return False
            elif value.lower() in ("true", "1", "yes", "y", "on"):
                return True

        # General transformation
        return transform_func(value)
    except (ValueError, TypeError) as e:
        logger.warning(f"Failed to transform {value} to {transform_func.__name__}: {e}")
        return value

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
    if hasattr(args, 'config') and args.config and isinstance(args.config, dict):
        config = deep_update(config, args.config)

    elif hasattr(args, 'config') and args.config and os.path.isfile(args.config):
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

    # Handle threshold parameters when segment is specified
    if 'segment' in cli_config:
        segment = cli_config['segment']
        thresholds = {}

        # Check for threshold parameters
        for cli_arg, config_key in THRESHOLD_MAPPING.items():
            if cli_arg in cli_config:
                # Get the transform function from PARAMETER_MAPPING
                transform_func = next((val[1] for key, val in PARAMETER_MAPPING.items()
                                      if key == cli_arg), None)

                # Transform the value
                if transform_func:
                    transformed_value = validate_and_cast(cli_config[cli_arg], transform_func)
                    thresholds[config_key] = transformed_value
                else:
                    thresholds[config_key] = cli_config[cli_arg]

                cli_config.pop(cli_arg)

        # Add thresholds to config if any were specified
        if thresholds:
            # Initialize the nested structure if not present
            if 'DEDUPLICATE' not in normalized_cli_config:
                normalized_cli_config['DEDUPLICATE'] = {}
            if 'THRESHOLDS' not in normalized_cli_config['DEDUPLICATE']:
                normalized_cli_config['DEDUPLICATE']['THRESHOLDS'] = {}

            # Set segment-specific thresholds
            normalized_cli_config['DEDUPLICATE']['THRESHOLDS'][segment] = thresholds

    # Process all other arguments using the parameter mapping
    for cli_arg, (config_path, transform_func) in PARAMETER_MAPPING.items():
        if cli_arg in cli_config:
            value = cli_config[cli_arg]

            # Transform the value
            value = validate_and_cast(value, transform_func)

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
    """
    parts = path.split('.')

    # Navigate to the deepest level
    current = config_dict
    for i, part in enumerate(parts[:-1]):
        if part not in current:
            current[part] = {}
        current = current[part]

    # Set the value at the deepest level
    current[parts[-1]] = value

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