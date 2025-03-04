from setuptools import setup, find_packages

setup(
    name="lasv-dedup",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "snakemake>=7.0.0",
        "pyyaml>=6.0",
    ],
    entry_points={
        "console_scripts": [
            "lasvdedup=lasvdedup.cli:main",
        ],
    },
    author="Joon-Klaps",
    author_email="your.email@example.com",
    description="LASV deduplication pipeline",
    keywords="bioinformatics, LASV, deduplication",
    python_requires=">=3.7",
)
