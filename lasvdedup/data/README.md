# Reference Data for LASV Deduplication

This directory contains reference data files used by the LASV deduplication pipeline:

- Reference alignments (*.aln files)
- Reference phylogenetic trees (*.treefile files)

These files are automatically used when no explicit `BASE_DATA_DIR` is specified.

Alignments are based on all available sequences on nucleotide Genbank that matched the search term ‘Lassa’ on 19/10/2021. Sequences were excluded if they did not correspond to the species Lassa mammarenavirus (taxonomy ID: 3052310), contained the term "hypothetical protein - modified" in their descriptor, or were shorter than 5000 bp for the L-segment or 2000 bp for the S-segment.

## File Format

Files should be named according to the pattern `{species}-{segment}.aln` and `{species}-{segment}.treefile`.

For example:
- `LASV-L.aln`: LASV L segment alignment
- `LASV-L.treefile`: LASV L segment phylogenetic tree
- `LASV-S.aln`: LASV S segment alignment
- `LASV-S.treefile`: LASV S segment phylogenetic tree

## Adding Custom Reference Data

You can specify a custom reference data directory using the `--ref-dir` argument:

```bash
lasvdedup --input contigs.tsv --seq-dir sequences/ --ref-dir /path/to/references/
```

Or by setting `BASE_DATA_DIR` in your config file.
