import os
import importlib.resources

def get_config_path():
    """Return the path to the default config.yaml file."""
    return str(importlib.resources.files('lasvdedup').joinpath('config.yaml'))

def get_snakefile_path():
    """Return the path to the default Snakefile."""
    return str(importlib.resources.files('lasvdedup').joinpath('Snakefile'))

def get_data_dir():
    """Return the path to the data directory."""
    return str(importlib.resources.files('lasvdedup').joinpath('data'))
