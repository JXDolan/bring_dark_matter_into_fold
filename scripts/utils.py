import os
import csv
import shutil
from Bio import SeqIO
from Bio import pairwise2
from pymol import cmd
import requests
import time
import glob
from pymol import cmd

def combine_fasta_folders(input_folders, output_folder):
    """
    Combines FASTA files from multiple folders into one, renaming them with labels like 'unique-query' or 'unique-hit'.
    """
    os.makedirs(output_folder, exist_ok=True)
    total = 0

    for folder in input_folders:
        label = "unique-query" if "query" in folder.lower() else "unique-hit" if "hit" in folder.lower() else "unknown"
        for root, _, files in os.walk(folder):
            for file in files:
                if file.endswith(".fasta") or file.endswith(".fa"):
                    src_path = os.path.join(root, file)
                    base, ext = os.path.splitext(file)
                    dst_path = os.path.join(output_folder, f"{base}-{label}{ext}")
                    shutil.copy(src_path, dst_path)
                    total += 1

    print(f"\nCopied and renamed {total} FASTA files into: {output_folder}")

def count_up_to_stop(seq):
    """
    Trims sequence at the first stop codon '*'.
    """
    return str(seq).split("*")[0]

def filter_and_trim_fasta(input_path, output_path, max_length=400):
    """
    Filters sequences that exceed max_length and trims at stop codons.
    """
    filtered_records = []
    original_count = 0

    for record in SeqIO.parse(input_path, "fasta"):
        original_count += 1
        trimmed_seq = count_up_to_stop(record.seq)
        record.seq = trimmed_seq
        if len(trimmed_seq) <= max_length:
            filtered_records.append(record)

    if filtered_records:
        SeqIO.write(filtered_records, output_path, "fasta")

    return original_count, len(filtered_records)

def batch_process_folder(input_folder, output_folder, max_length=400, label=""):
    """
    Runs `filter_and_trim_fasta` on all FASTA files in a folder.
    """
    os.makedirs(output_folder, exist_ok=True)
    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            in_path = os.path.join(input_folder, filename)
            out_path = os.path.join(output_folder, f"{os.path.splitext(filename)[0]}_filtered.fasta")
            total, kept = filter_and_trim_fasta(in_path, out_path, max_length)
            print(f"[{label}] {filename}: kept {kept}/{total} (≤ {max_length})")


def convert_csv_to_fasta(input_csv_path, output_fasta_path, id_column, seq_column):
    """
    Converts a CSV with ID and sequence columns to a FASTA file.
    """
    with open(input_csv_path, "r") as infile, open(output_fasta_path, "w") as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            outfile.write(f">{row[id_column]}\n{row[seq_column]}\n")

def batch_convert_csvs(input_folder, output_folder, id_column, seq_column):
    """
    Batch converts all CSVs in a folder to FASTA format.
    """
    os.makedirs(output_folder, exist_ok=True)
    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            in_path = os.path.join(input_folder, filename)
            out_path = os.path.join(output_folder, filename.replace(".csv", ".fasta"))
            convert_csv_to_fasta(in_path, out_path, id_column, seq_column)
            print(f"Wrote: {out_path}")

def csv_to_fasta_folder(input_file, output_file):
    """
    Converts only the query sequences from a result CSV to FASTA.
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            outfile.write(f">{row['query_protein_ID']}\n{row['query_protein_sequence']}\n")


def compute_sequence_identity(seq1, seq2):
    """
    Computes global sequence identity using a naive pairwise alignment.
    """
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    aligned_seq1 = best_alignment.seqA
    aligned_seq2 = best_alignment.seqB

    # NOTE: this is naive, but results are very close to mmseq outputs even though this is global
    matches = sum(res1 == res2 for res1, res2 in zip(aligned_seq1, aligned_seq2))
    identity = matches / len(aligned_seq1)
    return identity

def extract_unique_rows(input_file, id_column):
    """
    Extracts rows with unique entries in the specified column.
    """
    seen, unique = set(), []
    with open(input_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            uid = row[id_column]
            if uid not in seen:
                seen.add(uid)
                unique.append(row)
    return reader.fieldnames, unique

def process_all_csvs(input_folder, query_output_folder, hit_output_folder):
    """
    Splits CSVs into two new CSVs with unique query and hit protein IDs respectively.
    """
    os.makedirs(query_output_folder, exist_ok=True)
    os.makedirs(hit_output_folder, exist_ok=True)

    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            path = os.path.join(input_folder, filename)

            # Query processing
            q_fields, q_rows = extract_unique_rows(path, "query_protein_ID")
            with open(os.path.join(query_output_folder, filename), 'w', newline='') as qf:
                writer = csv.DictWriter(qf, fieldnames=q_fields)
                writer.writeheader()
                writer.writerows(q_rows)

            # Hit processing
            h_fields, h_rows = extract_unique_rows(path, "hit_protein_ID")
            with open(os.path.join(hit_output_folder, filename), 'w', newline='') as hf:
                writer = csv.DictWriter(hf, fieldnames=h_fields)
                writer.writeheader()
                writer.writerows(h_rows)

def extract_unique_first_column(file_path):
    """
    Returns a set of unique values from the first column of a TSV/CSV file.
    """
    unique = set()
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip():
                unique.add(line.split('\t')[0])
    return unique


def align_and_visualize_structures(output_dir, seq_id):
    """
    Loads and aligns PDBs using PyMOL, saving a `.pse` session.
    """
    cmd.reinitialize()

    pdb_files = sorted(glob.glob(os.path.join(output_dir, "*.pdb")))

    # TODO: here it uses the query structure (from ESMFold) as the reference structure
    ref_file = next((pdb for pdb in pdb_files if not os.path.basename(pdb).startswith("MGY") and not os.path.basename(pdb).startswith("ESMFold")), None)
    if not ref_file:
        print(f"No reference structure found for {seq_id}, skipping PyMOL visualization.")
        return
    else:
        print(f"Using {ref_file} as reference structure.")

    colors = ["cyan", "blue", "green", "yellow", "gray"]

    ref_name = os.path.splitext(os.path.basename(ref_file))[0]
    cmd.load(ref_file, ref_name)
    cmd.color("red", ref_name) # TODO: uses red as the reference structure color

    for i, pdb in enumerate(pdb_files):
        if pdb == ref_file:
            continue

        protein_name = os.path.splitext(os.path.basename(pdb))[0]
        cmd.load(pdb, protein_name)
        cmd.align(protein_name, ref_name)
        cmd.color(colors[i % len(colors)], protein_name) # colors looping supported

    cmd.show("cartoon")
    cmd.zoom()

    pymol_save_path = os.path.join(output_dir, "aligned_proteins_with_gt.pse")
    cmd.save(pymol_save_path)
    print(f"Aligned structures for {seq_id} saved to {pymol_save_path}")


def read_fasta(fasta_path):
    """
    Reads sequences from a FASTA file into a dictionary.
    """
    sequences = {}
    with open(fasta_path, "r") as file:
        name = None
        seq = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name:
            sequences[name] = "".join(seq)
    return sequences

def submit_esmfold(sequence, seq_id, max_retries):
    """
    Submits a sequence to the ESMFold API and retrieves the predicted PDB.
    """
    ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb"
    for attempt in range(max_retries):
        response = requests.post(ESMFOLD_API, data=sequence)
        if response.status_code != 200:
            print(f"ESMFold failed for {seq_id}: {response.text} (Attempt {attempt+1}/{max_retries})")
            time.sleep(30) # FIXME: wait time before trying again
            if attempt == max_retries - 1:
                return None
            else:
                continue

        return response.text

def submit_foldseek(pdb_content, seq_id, search_mode, max_retries):
    """
    Submits a PDB structure to Foldseek and returns a job ID.
    """
    FOLDSEEK_SUBMIT_URL = "https://search.foldseek.com/api/ticket"
    for attempt in range(max_retries):
        response = requests.post(FOLDSEEK_SUBMIT_URL, data={
            "q": pdb_content,
            "mode": search_mode,
            "database[]": "mgnify_esm30"
        })
        if response.status_code != 200:
            print(f"Foldseek submission failed for {seq_id}: {response.text} (Attempt {attempt+1}/{max_retries})")
            time.sleep(30) # FIXME: wait time before trying again
            if attempt == max_retries - 1:
                return None
            else:
                continue
    
    return response.json().get("id")

def check_foldseek_status(job_id, seq_id):
    """
    Monitors the Foldseek job until completion and returns the results.
    """
    FOLDSEEK_CHECK_STATUS_URL = "https://search.foldseek.com/api/ticket/"
    while True:
        status_response = requests.get(FOLDSEEK_CHECK_STATUS_URL + job_id)
        if status_response.status_code == 200:
            status_data = status_response.json()
            status = status_data.get("status")
            # print(f"Foldseek Job Status: {status}")
            if status == "COMPLETE":
                print(f"Foldseek job finished for {seq_id}. Fetching results...")
                break
            elif status == "ERROR":
                print(f"Foldseek job failed for {seq_id}.")
                return None
        else:
            print(f"Error checking Foldseek status for {seq_id}: {status_response.text}")
            return None
        
        time.sleep(1) # check frequency

    FOLDSEEK_RESULT_URL = "https://search.foldseek.com/api/result/"
    result_response = requests.get(FOLDSEEK_RESULT_URL + job_id + "/0") # NOTE: 0 is hardcoded, do not change

    if result_response.status_code != 200:
        print(f"Error fetching Foldseek results for {seq_id}: {result_response.text}")
        return None
    
    return result_response.json()