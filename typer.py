import argparse
import os
import csv
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

def run_blast(reference, genome_file, output_xml):
    blast_cline = NcbiblastnCommandline(query=reference, subject=genome_file, outfmt=5, out=output_xml)
    stdout, stderr = blast_cline()
    return output_xml

def find_best_references(blast_output):
    best_refs = set()
    best_score = 0
    with open(blast_output) as handle:
        records = NCBIXML.parse(handle)
        for record in records:
            query_id = record.query.split()[0]
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    score = hsp.score
                    if score > best_score:
                        best_score = score
                        best_refs = {query_id}
                    elif score == best_score:
                        best_refs.add(query_id)
    return list(best_refs), best_score

def extract_sequence(blast_output, genome_file, genome_name, ref_name):
    extracted = []
    with open(blast_output) as result_handle:
        records = NCBIXML.parse(result_handle)
        for record in records:
            if not record.query.startswith(ref_name):  # avoid mismatches
                continue
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    aligned_seq = hsp.sbjct
                    extracted.append((f"{genome_name}_{ref_name}", aligned_seq))
    return extracted

def run_mafft(input_fasta, output_fasta):
    subprocess.run(["mafft", "--auto", input_fasta], stdout=open(output_fasta, 'w'), check=True)

def run_fasttree(alignment_fasta, tree_output):
    subprocess.run(["fasttree", "-nt", alignment_fasta], stdout=open(tree_output, 'w'), check=True)

def filter_short_sequences(aligned_fasta, filtered_fasta, min_fraction=0.33):
    sequences = list(SeqIO.parse(aligned_fasta, "fasta"))
    if not sequences:
        raise ValueError("No sequences found in alignment!")

    reference_len = max(len(rec.seq) for rec in sequences)
    min_len = int(reference_len * min_fraction)

    filtered = [rec for rec in sequences if len(rec.seq.ungap('-')) >= min_len]

    if len(filtered) < 2:
        raise ValueError("Not enough sequences remaining after filtering for tree building.")

    SeqIO.write(filtered, filtered_fasta, "fasta")

def main(input_dir, reference_db, output_dir, build_tree=False):
    os.makedirs(output_dir, exist_ok=True)
    aligned_dir = os.path.join(output_dir, "aligned_gp60_sequences")
    blast_dir = os.path.join(output_dir, "blast_results")
    os.makedirs(aligned_dir, exist_ok=True)
    os.makedirs(blast_dir, exist_ok=True)
    log_csv = os.path.join(output_dir, "typing_log.csv")

    genome_files = [f for f in os.listdir(input_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    typing_results = {}
    all_gp60_seqs = []

    for genome_file in genome_files:
        genome_path = os.path.join(input_dir, genome_file)
        genome_name = os.path.splitext(genome_file)[0]

        # Step 1: Run multi-ref BLAST
        tmp_xml = os.path.join(blast_dir, f"{genome_name}_vs_refs.xml")
        run_blast(reference_db, genome_path, tmp_xml)

        best_refs, best_score = find_best_references(tmp_xml)
        typing_results[genome_file] = best_refs

        for ref in best_refs:
            blast_output = os.path.join(blast_dir, f"{genome_name}_vs_{ref}.xml")
            tmp_ref_file = os.path.join(blast_dir, f"{ref}.fasta")

            with open(tmp_ref_file, 'w') as rf_out:
                for record in SeqIO.parse(reference_db, "fasta"):
                    if record.id == ref:
                        SeqIO.write(record, rf_out, "fasta")
                        break

            run_blast(tmp_ref_file, genome_path, blast_output)
            aligned = extract_sequence(blast_output, genome_path, genome_name, ref)

            output_fasta = os.path.join(aligned_dir, f"{genome_name}.fasta")
            with open(output_fasta, 'w') as fasta_out:
                for header, seq in aligned:
                    fasta_out.write(f">{header}\n{seq}\n")

            all_gp60_seqs.extend(aligned)

    # Step 2: Write typing log
    with open(log_csv, 'w', newline='') as log_file:
        writer = csv.writer(log_file)
        writer.writerow(["Genome", "Best Reference(s)"])
        for genome, refs in typing_results.items():
            writer.writerow([genome, "|".join(refs)])

    # Step 3: Build tree if requested
    if build_tree:
        all_gp60_fasta = os.path.join(output_dir, "all_gp60s.fasta")
        aligned_fasta = os.path.join(output_dir, "aligned_gp60s.fasta")
        filtered_fasta = os.path.join(output_dir, "filtered_aligned_gp60s.fasta")
        tree_file = os.path.join(output_dir, "gp60_cladogram.nwk")

        with open(all_gp60_fasta, 'w') as out_all:
            for header, seq in all_gp60_seqs:
                out_all.write(f">{header}\n{seq}\n")
            for record in SeqIO.parse(reference_db, "fasta"):
                out_all.write(f">{record.id}\n{str(record.seq)}\n")

        run_mafft(all_gp60_fasta, aligned_fasta)
        filter_short_sequences(aligned_fasta, filtered_fasta)
        run_fasttree(filtered_fasta, tree_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract gp60 gene from genomes using best-matching reference.")
    parser.add_argument("-d", "--input_dir", required=True, help="Directory with genome FASTA files")
    parser.add_argument("-r", "--reference_db", required=True, help="Multi-fasta file of gp60 reference sequences")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory for output files")
    parser.add_argument("--tree", action="store_true", help="Build phylogenetic tree using MAFFT and FastTree")

    args = parser.parse_args()
    main(args.input_dir, args.reference_db, args.output_dir, args.tree)
