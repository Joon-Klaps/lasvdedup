import os
import pkg_resources

def get_config_path():
    """Return the path to the default config.yaml file."""
    return pkg_resources.resource_filename('lasvdedup', 'config.yaml')

def get_snakefile_path():
    """Return the path to the default Snakefile."""
    return pkg_resources.resource_filename('lasvdedup', 'Snakefile')

def get_data_dir():
    """Return the path to the data directory."""
    return pkg_resources.resource_filename('lasvdedup', 'data')
