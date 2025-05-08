import requests
import time
import json
import os
import glob
import argparse
from pymol import cmd
from tqdm import tqdm

from scripts.utils import read_fasta, submit_esmfold, submit_foldseek, check_foldseek_status, align_and_visualize_structures

# ESM and Foldseek APIs (NOTE: Hardcoded - do not change, as well as the url requests status codes)
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb"
ESMATLAS_FETCH_URL = "https://api.esmatlas.com/fetchPredictedStructure/"
FOLDSEEK_SUBMIT_URL = "https://search.foldseek.com/api/ticket"
FOLDSEEK_CHECK_STATUS_URL = "https://search.foldseek.com/api/ticket/"
FOLDSEEK_RESULT_URL = "https://search.foldseek.com/api/result/"

def fetch_pdb_from_esmatlas(target_name, seq_id, save_path):
    """
    Fetches a PDB structure from ESMAtlas and saves it
    """
    pdb_fetch_url = f"{ESMATLAS_FETCH_URL}{target_name}"
    response = requests.get(pdb_fetch_url)

    if response.status_code == 200:
        with open(save_path, "w") as pdb_file:
            pdb_file.write(response.text)
        return True
    else:
        print(f"Failed to fetch {target_name} from ESMAtlas for {seq_id}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", type=str, help="Path to the input FASTA file.", default="./fasta_files/test_sequence.fasta")
    parser.add_argument("--search_mode", type=str, choices=["3diaa", "tmalign"], help="Foldseek search mode: '3diaa' or 'tmalign'.", default="tmalign")
    parser.add_argument("--num_structures", type=int, help="Number of top pdb structures to save from Foldseek results.", default=5)
    parser.add_argument("--num_structures_json", type=int, help="Number of top resutls to save in json.", default=50)
    parser.add_argument("--output_dir", type=str, help="Output directory to save the results.", default="output_dir")
    parser.add_argument("--retry", type=int, help="Maximum number of retries for failed requests.", default=5)

    args = parser.parse_args()

    FASTA_FILE = args.fasta_file
    SEARCH_MODE = args.search_mode
    NUM_STRUCTURES = args.num_structures
    NUM_STRUCTURES_JSON = args.num_structures_json
    OUTPUT_BASE_DIR = args.output_dir
    MAX_RETRIES = args.retry

    sequences = read_fasta(FASTA_FILE)

    for seq_name, sequence in tqdm(sequences.items()):
        seq_id = seq_name.split("|")[2].split()[0] # FIXME: this way of extracting the seq ID might not work for some cases
        print(f"\nProcessing sequence: {seq_name}")

        output_dir = os.path.join(OUTPUT_BASE_DIR, seq_id)
        os.makedirs(output_dir, exist_ok=True)

        predicted_pdb_esmfold = submit_esmfold(sequence, seq_id, MAX_RETRIES)
        if not predicted_pdb_esmfold:
            print(f"ESMFold failed for {seq_id}.")
            continue

        # save predicted structure
        pdb_filename = os.path.join(output_dir, f"ESMFold_{seq_id}.pdb")
        with open(pdb_filename, "w") as pdb_file:
            pdb_file.write(predicted_pdb_esmfold)
        print(f"ESMFold predicted structure saved to: {pdb_filename}")

        # submit to Foldseek
        job_id = submit_foldseek(predicted_pdb_esmfold, seq_id, SEARCH_MODE, MAX_RETRIES)
        if not job_id:
            print(f"Foldseek submission failed for {seq_id}.")
            continue

        results = check_foldseek_status(job_id, seq_id)
        if not results:
            print(f"Foldseek running failed for {seq_id}.")
            continue

        fetched_results = []
        for result in results.get("results", []): # NOTE: return results structure are hardcoded - do not change
            for alignment_group in result.get("alignments", []):
                for i, alignment in enumerate(alignment_group):
                    target_name = alignment.get("target", "").replace(".pdb.gz", "")

                    if len(fetched_results) < NUM_STRUCTURES_JSON:
                        fetched_results.append(alignment)
                        if len(fetched_results) < NUM_STRUCTURES:
                            save_path = os.path.join(output_dir, f"{target_name}_rank{i+1}.pdb")

        output_json = os.path.join(output_dir, f"filtered_results_top_{NUM_STRUCTURES_JSON}.json")
        with open(output_json, "w") as outfile:
            json.dump(fetched_results, outfile, indent=4)
        print(f"Saved top results JSON to: {output_json}")

        align_and_visualize_structures(output_dir, seq_id)