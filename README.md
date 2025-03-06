# lasvdedup

A small utility to deduplicate LASV consensus genomes from the pipeline [viralgenie](https://github.com/Joon-Klaps/viralgenie).

## Installation

```bash
git clone https://github.com/Joon-Klaps/lasvdedup.git
cd lasvdedup
```

```bash
conda env create -f environment.yml
conda activate lasvdedup
```


## Usage

```bash
lasvdedup --contigs viralgenie/overview-tables/contigs_overview.tsv --seq-dir viralgenie/consensus/seq
```

