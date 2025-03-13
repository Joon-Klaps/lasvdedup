from setuptools import setup, find_packages

setup(
    name="lasvdedup",
    version="0.1.2",
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

    author="Your Name",
    description="A pipeline for deduplicating LASV sequences",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
