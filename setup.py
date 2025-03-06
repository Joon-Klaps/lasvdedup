from setuptools import setup, find_packages

setup(
    name="lasv-dedup",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "pyyaml>=6.0",
        "biopython>=1.79",
        "pandas>=1.3.3",
        "numpy>=1.21.2",
        "baltic>=0.3.0"
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
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
