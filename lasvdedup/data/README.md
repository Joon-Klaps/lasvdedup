# Reference Data for LASV Deduplication

This directory contains reference data files used by the LASV deduplication pipeline:

- Reference alignments (*.aln files)
- Reference phylogenetic trees (*.treefile files)

These files are automatically used when no explicit `BASE_DATA_DIR` is specified.

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
