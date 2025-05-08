import argparse
import os
import pandas as pd
import requests
import json
from tqdm import tqdm
from scripts.utils import submit_foldseek, check_foldseek_status

TESTING = False #FIXME: set to False for actual runs, when set to true, only max of 5 sequences will be processed

FOLDSEEK_SUBMIT_URL = "https://search.foldseek.com/api/ticket"
FOLDSEEK_CHECK_STATUS_URL = "https://search.foldseek.com/api/ticket/"
FOLDSEEK_RESULT_URL = "https://search.foldseek.com/api/result/"


def download_alphafold_structures(folder_name, input_file, max_retries):

    os.makedirs(folder_name, exist_ok=True)
    os.makedirs(os.path.join(folder_name, "query_pdbs"), exist_ok=True)

    df = pd.read_csv(input_file, sep='\t')

    if 'Entry' not in df.columns:
        raise ValueError("Input file does not contain 'Entry' column.")

    base_url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb"

    for i, uniprot_id in enumerate(tqdm(df['Entry'].unique())):

        if TESTING:
            if i >= 5:
                break

        url = base_url.format(uniprot_id)
        output_path = os.path.join(folder_name, "query_pdbs", f"{uniprot_id}.pdb")

        for i in range(max_retries):

            response = requests.get(url)
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded: {uniprot_id}")
                break
            else:
                print(f"Failed to download {uniprot_id} (status code {response.status_code}), retrying # {i}...")


def run_foldseek_on_pdbs(folder_name, search_mode, num_structures_json, max_retries):

    pdb_folder = os.path.join(folder_name, "query_pdbs")
    json_folder = os.path.join(folder_name, "foldseek_jsons")
    os.makedirs(json_folder, exist_ok=True)

    pdb_files = [f for f in os.listdir(pdb_folder) if f.endswith(".pdb")]

    for pdb_file in tqdm(pdb_files, desc="Running Foldseek on query pdbs"):
        seq_id = os.path.basename(pdb_file).split(".")[0]
        pdb_path = os.path.join(pdb_folder, pdb_file)

        with open(pdb_path, "r") as f:
            pdb_content = f.read()

        job_id = submit_foldseek(pdb_content, seq_id, search_mode, max_retries)

        if not job_id:
            print(f"Submission failed for {seq_id}. Skipping.")
            continue

        results = check_foldseek_status(job_id, seq_id)
        if not results:
            print(f"Foldseek running failed for {seq_id}.")
            continue

        fetched_results = []
        for result in results.get("results", []):
            for alignment_group in result.get("alignments", []):
                for i, alignment in enumerate(alignment_group):
                    if len(fetched_results) < num_structures_json:
                        fetched_results.append(alignment)

        output_json_path = os.path.join(json_folder, f"{seq_id}_{num_structures_json}.json")
        with open(output_json_path, "w") as f:
            json.dump(fetched_results, f, indent=4)
        print(f"Saved top {len(fetched_results)} Foldseek hits for {seq_id}")


def parse_foldseek_jsons_to_dataframe(folder_name, num_structures_json):

    json_folder = os.path.join(folder_name, "foldseek_jsons")
    all_records = []

    for file_name in os.listdir(json_folder):

        if not file_name.endswith(".json"):
            continue

        query_id = file_name.replace(".json", "").split("_")[0]
        file_path = os.path.join(json_folder, file_name)

        with open(file_path, "r") as f:
            alignments = json.load(f)

        for alignment in alignments:
            record = {
                "QueryID": query_id,
                "MGnifyID": alignment.get("target").replace(".pdb.gz", ""),
                "eval": alignment.get("eval"),
                "SeqIdentity": alignment.get("seqId"),
                "Sequence": alignment.get("tSeq")
            }
            all_records.append(record)

    df = pd.DataFrame(all_records)
    output_path = os.path.join(folder_name, f"foldseek_results_top_{num_structures_json}_all.tsv")
    df.to_csv(output_path, sep='\t', index=False)

    print(f"Saved parsed results csv to: {output_path}")


def filter_foldseek_dataframe(folder_name, search_mode, num_top, threshold):

    input_path = os.path.join(folder_name, f"foldseek_results_top_{num_top}_all.tsv")
    output_path = os.path.join(folder_name, f"foldseek_results_top_{num_top}_filtered_{threshold}.tsv")

    if not os.path.exists(input_path):
        print(f"Foldseek result tsv file not found at: {input_path}")
        return

    df = pd.read_csv(input_path, sep='\t')

    if search_mode == "3diaa":
        filtered_df = df[df["eval"] <= threshold]
    elif search_mode == "tmalign":
        filtered_df = df[df["eval"] >= threshold]
    else:
        raise ValueError("Invalid search mode for filtering.")

    filtered_df.to_csv(output_path, sep='\t', index=False)

    print(f"Saved filtered results to: {output_path}.")


def export_sequences_to_fastas(folder_name, input_file, threshold, num_top):
    """
    Save the UniProt and Foldseek sequences to FASTA files.
    """

    df_uniprot = pd.read_csv(input_file, sep='\t')
    if 'Entry' not in df_uniprot.columns or 'Sequence' not in df_uniprot.columns:
        raise ValueError("Input TSV must have 'Entry' and 'Sequence' columns.")

    # Load Foldseek filtered results
    filtered_foldseek_path = os.path.join(folder_name, f"foldseek_results_top_{num_top}_filtered_{threshold}.tsv")
    df_foldseek = pd.read_csv(filtered_foldseek_path, sep='\t')

    # Output paths
    uniprot_fasta_path = os.path.join(folder_name, "uniprot.fasta")
    foldseek_fasta_path = os.path.join(folder_name, f"foldseek_esm_top_{num_top}_filtered_{threshold}.fasta")
    combined_fasta_path = os.path.join(folder_name, f"combined_top_{num_top}_filtered_{threshold}.fasta")

    # Write uniprot.fasta
    with open(uniprot_fasta_path, "w") as f:
        for _, row in df_uniprot.iterrows():
            f.write(f">UNIPROT_{row['Entry']}\n{row['Sequence']}\n")

    # Write foldseek_esm.fasta
    with open(foldseek_fasta_path, "w") as f:
        for _, row in df_foldseek.iterrows():
            f.write(f">{row['MGnifyID']}\n{row['Sequence']}\n")

    # Write combined.fasta
    with open(combined_fasta_path, "w") as fout:
        with open(uniprot_fasta_path, "r") as f1:
            fout.write(f1.read())
        with open(foldseek_fasta_path, "r") as f2:
            fout.write(f2.read())

    print(f"Saved UniProt sequences to: {uniprot_fasta_path}")
    print(f"Saved Foldseek sequences to: {foldseek_fasta_path}")
    print(f"Saved combined sequences to: {combined_fasta_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download AlphaFold PDBs from UniProt ID list.")
    parser.add_argument("--input_file", type=str, help="Path to the input TSV file.", default="./tsvs/uniprotkb_sqr_reviewed.tsv")
    parser.add_argument("--output_folder_name", type=str, help="Name of the output folder to save PDB files.", default="sqr")

    parser.add_argument("--search_mode", type=str, choices=["3diaa", "tmalign"], help="Foldseek search mode: '3diaa' or 'tmalign'.", default="tmalign")
    parser.add_argument("--num_structures_json", type=int, help="Number of top resutls to save in json.", default=100)
    parser.add_argument("--score_threshold", type=float, help="Threshold to filter Foldseek results.", default=0.7) #NOTE: use E-value for 3diaa and TM-score for tmalign

    parser.add_argument("--retry", type=int, help="Maximum number of retries for failed requests.", default=3)

    args = parser.parse_args()

    if TESTING: print(f"WARNING: TESTING IS TURNED TO TRUE. ONLY MAX OF 5 SEQUENCES WILL BE PROCESSED.")

    # get all alphafold structure pdbs for thr given uniprot ids
    download_alphafold_structures(args.output_folder_name, args.input_file, args.retry)

    # run foldseek on the downloaded pdbs
    run_foldseek_on_pdbs(args.output_folder_name, args.search_mode, args.num_structures_json, args.retry)

    # parse foldseek data and save to a df/tsv
    parse_foldseek_jsons_to_dataframe(args.output_folder_name, args.num_structures_json)

    # filter df using specified threshold
    # filter_foldseek_dataframe(args.output_folder_name, args.search_mode, args.num_structures_json, args.score_threshold)

    # export sequences to fastas
    # export_sequences_to_fastas(args.output_folder_name, args.sinput_file, args.score_threshold, args.num_structures_json)