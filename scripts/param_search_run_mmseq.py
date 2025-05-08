import pandas as pd
import argparse
import os
import subprocess

from tqdm import tqdm

def export_sequences_to_fastas_with_filter(base_folder_name, input_uniprot_fasta, input_foldseek_tsv, search_mode, thresholds, num_top = 100):
    """
    Save the Filtered Foldseek sequences to different FASTA files.
    """
    output_folder = os.path.join(base_folder_name, "filtered_fasta")
    os.makedirs(output_folder, exist_ok=True)

    if not os.path.exists(input_foldseek_tsv) or not os.path.exists(input_uniprot_fasta):
        print("Foldseek or iniprot result tsv file not found.")
        return

    df_foldseek = pd.read_csv(input_foldseek_tsv, sep='\t')

    for threshold in thresholds:
        print(threshold)

        if search_mode == "3diaa":
            filtered_df = df_foldseek[df_foldseek["eval"] <= threshold]
        elif search_mode == "tmalign":
            filtered_df = df_foldseek[df_foldseek["eval"] >= threshold]
        else:
            raise ValueError("Invalid search mode for filtering.")

        # Output paths
        foldseek_fasta_path = os.path.join(base_folder_name, "filtered_fasta", f"foldseek_esm_top_{num_top}_filtered_{threshold}.fasta")
        combined_fasta_path = os.path.join(base_folder_name, "filtered_fasta", f"combined_top_{num_top}_filtered_{threshold}.fasta")

        # Write foldseek_esm.fasta
        with open(foldseek_fasta_path, "w") as f:
            for _, row in filtered_df.iterrows():
                f.write(f">{row['MGnifyID']}\n{row['Sequence']}\n")

        # Write combined.fasta
        with open(combined_fasta_path, "w") as fout:
            with open(input_uniprot_fasta, "r") as f1:
                fout.write(f1.read())
            with open(foldseek_fasta_path, "r") as f2:
                fout.write(f2.read())

        print(f"Saved Foldseek sequences to: {foldseek_fasta_path}")
        print(f"Saved combined sequences with tm score threshold {threshold} to: {combined_fasta_path}")

def run_mmseq(base_folder, tm_thresholds, eval_thresholds, protein_file_path, num_top = 100):

    mmseq_db_base_path = os.path.join(base_folder, "mmseq_dbs")
    os.makedirs(mmseq_db_base_path, exist_ok=True)
    mmseq_tmps_path = os.path.join(mmseq_db_base_path, "tmps")
    os.makedirs(mmseq_tmps_path, exist_ok=True)
    
    try:
        os.makedirs(os.path.join(mmseq_db_base_path, "db_query"), exist_ok=True)
        mmseq_query_db = os.path.join(mmseq_db_base_path, "db_query", "query")
        subprocess.run(["mmseqs", "createdb", f"{protein_file_path}", f"{mmseq_query_db}"])
    except Exception as e:
        print(f"Creating mmseq db for query proteins failed for {protein_file_path} at {mmseq_query_db}: {e}")

    for tm_threshold in tqdm(tm_thresholds):
        os.makedirs(os.path.join(mmseq_db_base_path, f"db_{tm_threshold}"), exist_ok=True)
        mmseq_combined_db_path = os.path.join(mmseq_db_base_path, f"db_{tm_threshold}", f"db_{tm_threshold}")
        try:
            subprocess.run(["mmseqs", "createdb", f"{base_folder}/filtered_fasta/combined_top_{num_top}_filtered_{tm_threshold}.fasta", f"{mmseq_combined_db_path}"])
        except Exception as e:
            print(f"Creating mmseq db for database failed at {mmseq_combined_db_path}: {e}")

        # loop through e value thresholds and run, this is fast. 
        # if efficiency is a concern, run it once and then manually filter the results using the the thresholds would be faster
        for eval_threshold in eval_thresholds:
            try:
                mmseq_tmp_path = os.path.join(mmseq_tmps_path, f"tmp_{tm_threshold}_{eval_threshold}")
                
                mmseq_result_tmp_base_path = os.path.join(base_folder, "mmseq_results_non_readable")
                os.makedirs(mmseq_result_tmp_base_path, exist_ok=True)
                mmseq_result_tmp_path = os.path.join(mmseq_result_tmp_base_path, f"mmseq_result_{tm_threshold}_{eval_threshold}")
                os.makedirs(mmseq_result_tmp_path, exist_ok=True)

                subprocess.run([f"mmseqs", "search", "-e", f"{eval_threshold}", f"{mmseq_query_db}", f"{mmseq_combined_db_path}", 
                                f"{mmseq_result_tmp_path}/res_{tm_threshold}_{eval_threshold}", f"{mmseq_tmp_path}"])
            except Exception as e:
                print(f"Running mmseq search failed for {mmseq_query_db}/{mmseq_combined_db_path}: {e}")
            
            try:
                readble_results_base = os.path.join(base_folder, "mmseq_results")
                os.makedirs(readble_results_base, exist_ok=True)
                save_result_path = os.path.join(readble_results_base, f"result_tm_{tm_threshold}_eval_{eval_threshold}.m8")
                subprocess.run(["mmseqs", "convertalis", f"{mmseq_query_db}", f"{mmseq_combined_db_path}",
                                 f"{mmseq_result_tmp_path}/res_{tm_threshold}_{eval_threshold}", f"{save_result_path}"])
            except Exception as e:
                print(f"Converting mmseq results to readble format failed: {mmseq_query_db}/{mmseq_combined_db_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #NOTE: Assuming you already have a file like sqr/foldseek_results_top_100_all.tsv file, if not, run `build_combined_database.py` first.
    parser.add_argument("--foldseek_input_file", type=str, help="Path to the input foldseek TSV file.", default="sqr/foldseek_results_top_100_all.tsv")
    parser.add_argument("--MAG_input_file", type=str, help="Path to the input MAG fasta file.", default = "./fasta_files/sqr_combined_100.fasta")
    parser.add_argument("--MAG_input_name", type=str, help="input MAG file name.", default = "") #NOTE: If not passed in, it uses the folder name of the protein file as the name

    parser.add_argument("--search_mode", type=str, choices=["3diaa", "tmalign"], help="Foldseek search mode: '3diaa' or 'tmalign'.", default="tmalign")

    parser.add_argument('--tm_thresholds', nargs='+', type=float, help="List of tm align score thresholds for building combined db.", 
                        default=[0.7, 0.725, 0.75, 0.775, 0.8, 0.805, 0.81, 0.815, 0.82, 0.825, 0.83, 0.835, 0.84, 0.845, 0.85, 0.855, 0.86, 0.865, 
                                 0.87, 0.875, 0.88, 0.885, 0.89, 0.895, 0.9, 0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 
                                 0.965, 0.97, 0.975, 0.98, 0.985, 0.99])
    parser.add_argument('--eval_thresholds', nargs='+', type=float, help="List of E value thresholds for filtering mmseq results.",
                        default=[1e-3, 1e-4, 1e-5])
    args = parser.parse_args()

    base_folder = os.path.join(os.path.dirname(args.foldseek_input_file), args.MAG_input_name)
    if len(args.MAG_input_name) == 0:
        args.MAG_input_name = args.MAG_input_file.split('/')[-2]

    #FIXME: This is the ground truth proteins from uniprot and is hard-coded here as the 3 swissprot sqr proteins
    input_uniprot_fasta = os.path.join(os.path.dirname(args.foldseek_input_file), "uniprot.fasta")

    export_sequences_to_fastas_with_filter(base_folder_name = base_folder, input_uniprot_fasta = input_uniprot_fasta, 
                                           input_foldseek_tsv = args.foldseek_input_file, search_mode = args.search_mode, thresholds = args.tm_thresholds)

    run_mmseq(base_folder = base_folder, tm_thresholds = args.tm_thresholds, eval_thresholds = args.eval_thresholds, protein_file_path = args.MAG_input_file)
    
    