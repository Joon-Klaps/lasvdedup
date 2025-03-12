from setuptools import setup, find_packages

setup(
    name="lasvdedup",
    version="0.1.1",
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.9, <3.13',
    install_requires=[
        "biopython",
        "numpy",
        "pandas",
        "pyyaml",
        "phylodm",
    ],
    entry_points={
        "console_scripts": [
            "lasvdedup=lasvdedup.cli:main",
        ],
    },
)
