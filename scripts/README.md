## `run_foldseek_and_visualize.py`
This script performs end-to-end protein structure prediction, similarity search, and visualization using **ESMFold** and **Foldseek**.
### Features
1. Predicts all 3D protein structures from FASTA sequences using the ESMFold API.
2. Submits predicted PDBs to Foldseek for structural similarity search.
3. Saves Foldseek results as JSON.
4. (Optional) Downloads top hit structures from ESMAtlas.
5. Aligns and visualizes structures in PyMOL, saving a `.pse` session.
### Usage
```
python ./scripts/run_foldseek_and_visualize.py \
  --fasta_file ./fasta_files/test_sequence.fasta \
  --search_mode tmalign \
  --num_structures 5 \
  --num_structures_json 50 \
  --output_dir output_dir \
  --retry 3
```
### Arguments

- `--fasta_file`: Path to input FASTA file.
- `--search_mode`: Foldseek mode (`tmalign` or `3diaa`).
- `--num_structures`: Number of top PDB hits to download and align.
- `--num_structures_json`: Number of top Foldseek hits to save in JSON.
- `--output_dir`: Output folder to save structures and results.
- `--retry`: Max retries for failed API requests.


## `build_combined_database.py`

This script automates the process of building the **base** of a database expanded with ESM Atlas by:

1. **Downloading AlphaFold PDBs** from a list of UniProt IDs (these would be the query proteins).
2. **Submitting PDBs to Foldseek and Search ESM Atlas** using either the `tmalign` or `3diaa` search mode.
3. **Retrieving and parsing Foldseek results** into a structured TSV format.
4. **(Optional)** Filtering results based on a user-defined score threshold.
5. **(Optional)** Exporting sequences to FASTA files for downstream analysis.

The functionality of Number 4. and 5. are included and automated in the `param_search` scripts series.

### Usage
```
python ./scripts/build_combined_database.py \
  --input_file ./tsvs/uniprotkb_sqr_reviewed.tsv \
  --output_folder_name sqr \
  --search_mode tmalign \
  --num_structures_json 100 \
  --score_threshold 0.7 \
  --retry 5
```
### Arguments
- `--input_file`: Path to TSV file with UniProt `Entry` and `Sequence` columns.
- `--output_folder_name`: Output folder to save results (this is the base folder name, for example 'sqr').
- `--search_mode`: Foldseek mode (`tmalign` or `3diaa`).
- `--num_structures_json`: Number of top Foldseek results to save per query.
- `--score_threshold`: Score threshold for filtering (TM-score or E-value).
- `--retry`: Max number of retries for API requests.


# Running Grid Parameter Search

Running the three scripts in order performs a full parametric evaluation of MMseqs2 search quality:

1. `param_search_run_mmseq.py` builds databases and performs MMseqs2 similarity searches using various TM-score and E-value thresholds.
2. `param_search_process_mmseq_res.py` filters results using various sequence coverages.
3. `param_search_calculate_PR.py` evaluates performance (precision and recall) by comparing predictions to a known positive set.
4. `scripts/plot.ipynb` to obtain the heatmap of the parameter search.

**NOTE**: A base foldseek results file for the combined database will be needed, so run `build_combined_database.py` first and set the thresholds as loose as possible. Then this series of scripts will perform further filtering based on the `build_combined_database.py` results using tigher TM-score thresholds.

## `param_search_run_mmseq.py`
This script prepares and runs **MMseqs2** similarity searches across multiple thresholds:
### Features
1. Filters Foldseek results by TM-score or E-value thresholds and exports filtered FASTA files (building combined database using different thresholds).
2. Builds MMseqs2 databases for each thresholded FASTA set.
3. Performs MMseqs2 similarity searches for each database using a given query protein set.
4. Converts results to a human-readable `.m8` format.
### Usage
```
python ./scripts/param_search_run_mmseq.py \
  --foldseek_input_file ./sqr/foldseek_results_top_100_all.tsv \
  --MAG_input_file ./fasta_files/sqr_combined_100.fasta \
  --MAG_input_name sqr_combined_100 \
  --search_mode tmalign \
  --tm_thresholds 0.7 0.75 0.8 0.85 0.9 \
  --eval_thresholds 1e-3 1e-5
```

## `param_search_process_mmseq_res.py`
This script parses MMseqs2 `.m8` output and filters by **sequence coverage**:
### Features
1. Extracts aligned sequence info from `.m8` MMseqs2 results.
2. Computes coverage of alignment over hit sequence length.
3. Saves filtered results to `.csv` files based on coverage thresholds.
### Usage
```
python ./scripts/param_search_process_mmseq_res.py \
  --base_folder ./sqr/foldseek_results_top_100_all.tsv \
  --coverages 0.7 0.75 0.8 0.85 0.9 \
  --query_sequences_path ./fasta_files/sqr_combined_100.fasta
```

## `param_search_calculate_PR.py`
This script calculates **precision** and **recall** from MMseqs2 results based on a true positive set:
### Features
1. Extracts predicted query IDs from coverage-filtered `.csv` files.
2. Compares them against a known true positive FASTA set.
3. Computes and summarizes precision/recall per result file.
### Usage
```
python ./scripts/param_search_calculate_PR.py \
  --csv_folder ./sqr/sqr_combined_100/coverage_filtered_csv \
  --true_positive_fasta ./fasta_files/sqr_positive_50.fasta \
  --out precision_recall_summary.csv
```

## `utils.py`

A general-purpose utility module containing helper functions for sequence processing, CSV/FASTA conversion, structure visualization, and API calls.

- `combine_fasta_folders(...)`  
  Combines FASTA files from different folders into a single directory with labeled filenames.

- `count_up_to_stop(seq)`  
  Trims a protein sequence at the first stop codon.

- `filter_and_trim_fasta(...)`  
  Trims sequences and filters them based on a max length threshold.

- `batch_process_folder(...)`  
  Applies `filter_and_trim_fasta` to all FASTA files in a directory.

- `convert_csv_to_fasta(...)`  
  Converts a CSV with sequence info into a FASTA file.

- `batch_convert_csvs(...)`  
  Batch version of `convert_csv_to_fasta`.

- `csv_to_fasta_folder(...)`  
  Extracts query sequences from a result CSV into a FASTA file.

- `compute_sequence_identity(seq1, seq2)`  
  Computes global sequence identity via pairwise alignment.

- `extract_unique_rows(...)`  
  Returns unique rows based on a specified column in a CSV.

- `process_all_csvs(...)`  
  Processes CSVs to generate query-unique and hit-unique files.

- `extract_unique_first_column(file_path)`  
  Extracts all unique values from the first column of a file.

- `align_and_visualize_structures(...)`  
  Uses PyMOL to align and visualize structures in a directory.

- `read_fasta(fasta_path)`  
  Loads sequences from a FASTA file into a Python dictionary.

- `submit_esmfold(...)`  
  Submits a sequence to the ESMFold API to obtain a predicted PDB.

- `submit_foldseek(...)`  
  Submits a structure to Foldseek and retrieves a job ID.

- `check_foldseek_status(...)`  
  Monitors a Foldseek job until completion and fetches the results.