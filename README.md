# Tbio/UBC: Bringing functional dark matter annotations into the fold

This repository provides a protein analysis pipeline focused on identifying structural homologs for specified proteins, hence enhancing the coverage and recall of protein function annotations. It mainly leverages ESMFold, Foldseek, and MMseqs2, with support for structural alignment, parameter grid search, and precision-recall evaluation.

- `scripts/` contains all core Python scripts for structure prediction, Foldseek submission, MMseqs2 runs, small helper utilities and visualization tools.
- `fasta_files/` includes example FASTA inputs for SQR proteins.
- `tsvs/` contains the UniProtKB SQR search results as a TSV file.

## Create and activate environment
```
conda env create -f env.yml
conda activate functional-fold-env
```