import argparse
import os
import csv
import sys

def find_sequence_in_fasta(fasta_path, target_id):
    seq = []
    found = False
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                header_id = line[1:].split()[0]
                if header_id == target_id:
                    found = True
                    seq = []
                elif found:
                    break
            elif found:
                seq.append(line)
    return ''.join(seq) if found else "MISSING"

def extract_m8_to_rows(m8_path, query_fa, hit_fa):
    rows = []
    with open(m8_path) as f:
        for raw in f:
            raw = raw.strip()
            if not raw or raw.startswith('#'):
                continue
            cols = raw.split('\t')
            if len(cols) < 11:
                sys.stderr.write(f"Skipping malformed line: {raw}\n")
                continue

            qid, hid = cols[0], cols[1]
            seqid = cols[2]
            alilen = cols[3]
            evalue = cols[10]

            qseq = find_sequence_in_fasta(query_fa, qid)
            hseq = find_sequence_in_fasta(hit_fa, hid)

            rows.append({
                "query_protein_ID": qid,
                "query_protein_sequence": qseq,
                "hit_protein_ID": hid,
                "hit_protein_sequence": hseq,
                "sequence_id": seqid,
                "alignment_length": alilen,
                "e_value": evalue
            })
    return rows

def filter_rows_by_coverage(rows, threshold):
    filtered = []
    for row in rows:
        try:
            aln_len = float(row['alignment_length'])
            hit_seq = row['hit_protein_sequence']
            hit_len = len(hit_seq)
            if hit_len == 0:
                raise ValueError("Hit sequence length is zero")
            coverage = aln_len / hit_len
            if coverage >= threshold:
                row['coverage'] = f"{coverage:.3f}"
                filtered.append(row)
        except Exception as e:
            sys.stderr.write(f"Skipping row due to error ({e}): {row}\n")
    return filtered

def write_csv(rows, output_path):
    fieldnames = [
        "query_protein_ID",
        "query_protein_sequence",
        "hit_protein_ID",
        "hit_protein_sequence",
        "sequence_id",
        "alignment_length",
        "e_value",
        "coverage"
    ]
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_folder', type=str, help="Folder with mmseq_results and filtered_fasta subfolders.", default="sqr/GCA_020348965.1")
    parser.add_argument('--coverages', nargs='+', type=float, 
                        help="Coverage thresholds to filter by.", 
                        default=[0.7, 0.725, 0.75, 0.775,
                                 0.8, 0.805, 0.81, 0.815, 0.82, 0.825, 0.83, 0.835, 0.84, 0.845, 0.85, 0.855, 0.86, 0.865, 
                                 0.87, 0.875, 0.88, 0.885, 0.89, 0.895,
                                 0.9, 0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 
                                 0.965, 0.97, 0.975, 0.98, 0.985, 0.99])
    parser.add_argument('--query_sequences_path', help="FASTA file of query sequences", default="fasta_files/sqr_combined_100.fasta")
    args = parser.parse_args()

    mmseq_result_folder = os.path.join(args.base_folder, "mmseq_results")
    fasta_folder = os.path.join(args.base_folder, "filtered_fasta")
    output_folder = os.path.join(args.base_folder, "coverage_filtered_csv")
    os.makedirs(output_folder, exist_ok=True)

    all_m8_files = [os.path.join(mmseq_result_folder, f) for f in os.listdir(mmseq_result_folder) if f.endswith(".m8")]

    for m8_file in all_m8_files:
        m8_basename = os.path.splitext(os.path.basename(m8_file))[0]

        hit_fasta_path = os.path.join(fasta_folder, f"combined_top_100_filtered_0.7.fasta")
        rows = extract_m8_to_rows(m8_file, args.query_sequences_path, hit_fasta_path)

        for cov in args.coverages:
            filtered_rows = filter_rows_by_coverage(rows, cov)
            out_name = f"{m8_basename}_cov_{cov}.csv"
            out_path = os.path.join(output_folder, out_name)
            write_csv(filtered_rows, out_path)
            print(f"Saved {len(filtered_rows)} rows to {out_name}")

if __name__ == "__main__":
    main()