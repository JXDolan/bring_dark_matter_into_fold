import os
import csv
import argparse

def extract_ids(fasta_text):
    ids = []
    for line in fasta_text.splitlines():
        if line.startswith('>'):
            # Take the first "word" after ">", split by spaces
            id_part = line.split()[0][1:]  # remove '>' and keep ID
            ids.append(id_part)
    return ids

def compute_precision_recall(predicted_ids, true_ids):
    if not true_ids:  #NOTE: no true positives
        if not predicted_ids:
            return 1.0, 1.0  # both scores perfect
        else:
            return 0.0, 1.0  # predicted something, but there were no true positives

    if not predicted_ids:
        return 0.0, 0.0  # no predictions when there were positives

    true_positives = predicted_ids & true_ids
    precision = len(true_positives) / len(predicted_ids)
    recall = len(true_positives) / len(true_ids)
    return precision, recall

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv_folder", type=str, required=True, help="Folder containing filtered CSV files.")
    parser.add_argument("--true_positive_fasta", type=str, required=True, help="FASTA file of true positive query IDs.")
    parser.add_argument("--out", type=str, default="precision_recall_summary.csv")
    args = parser.parse_args()

    csv_files = [f for f in os.listdir(args.csv_folder) if f.endswith(".csv") and not f.startswith("precision")]
    summary_path = os.path.join(args.csv_folder, args.out)

    with open(args.true_positive_fasta, 'r') as f:
        fasta_text = f.read()

    true_positive_ids = set(extract_ids(fasta_text))
    print(f"Loaded {len(true_positive_ids)} true positive IDs.")

    with open(summary_path, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(["filename", "num_predicted_queries", "num_true_positives", "precision", "recall"])

        for file in csv_files:
            file_path = os.path.join(args.csv_folder, file)
            try:
                with open(file_path, newline='') as f:
                    reader = csv.DictReader(f)
                    predicted_query_ids = {row['query_protein_ID'] for row in reader}
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
                continue

            tp_matches = predicted_query_ids & true_positive_ids
            precision, recall = compute_precision_recall(predicted_query_ids, true_positive_ids)

            writer.writerow([
                file,
                len(predicted_query_ids),
                len(tp_matches),
                f"{precision:.3f}",
                f"{recall:.3f}"
            ])

    print(f"Evaluation complete. Results saved to {summary_path}")

if __name__ == "__main__":

    main()
