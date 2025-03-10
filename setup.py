from setuptools import setup, find_packages

setup(
    name="lasvdedup",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "biopython",
        "numpy",
        "pandas",
        "pyyaml",
        "phylodm",
        # Other dependencies
    ],
    entry_points={
        "console_scripts": [
            "lasvdedup=lasvdedup.cli:main",
        ],
    },
    # Other metadata
)
