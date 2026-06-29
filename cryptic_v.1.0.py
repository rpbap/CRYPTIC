#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


DNA_CHARS = set("ACGTRYSWKMBDHVN-")


@dataclass
class RefMeta:
    ref_id: str
    species: str
    species_raw: str
    gene: str
    gp60_family: str
    gp60_variant: str
    gp60_subtype: str
    accession: str
    length_bp: int
    original_header: str


@dataclass
class FastaRecord:
    id: str
    desc: str
    seq: str


@dataclass
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int
    qseq: str
    sseq: str

    @property
    def qcov(self) -> float:
        return 100.0 * self.length / self.qlen if self.qlen else 0.0

    @property
    def scov(self) -> float:
        return 100.0 * self.length / self.slen if self.slen else 0.0

    @property
    def strand(self) -> str:
        return "minus" if self.sstart > self.send else "plus"


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def log(msg: str, verbose: bool = False) -> None:
    if verbose:
        print(f"[gp60typer] {msg}", file=sys.stderr, flush=True)


def require_executable(name: str) -> None:
    if shutil.which(name) is None:
        die(f"Required executable not found in PATH: {name}")


def clean_seq(seq: str) -> str:
    seq = str(seq).upper().replace("U", "T")
    return "".join(ch for ch in seq if ch in DNA_CHARS and ch != "-")


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    return clean_seq(seq).translate(comp)[::-1]


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def parse_fasta(path: Path) -> Iterable[FastaRecord]:
    header = None
    seq_lines: List[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    first = header.split()[0]
                    yield FastaRecord(first, header, clean_seq("".join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        first = header.split()[0]
        yield FastaRecord(first, header, clean_seq("".join(seq_lines)))


def write_fasta(records: Iterable[Tuple[str, str, Optional[str]]], path: Path) -> None:
    with path.open("w") as out:
        for rid, seq, desc in records:
            seq = clean_seq(seq)
            if not seq:
                continue
            if desc:
                out.write(f">{rid} {desc}\n{wrap(seq)}\n")
            else:
                out.write(f">{rid}\n{wrap(seq)}\n")


def load_single_json_db(db_path: Path, work_dir: Path, verbose: bool = False) -> Tuple[Path, Dict[str, RefMeta]]:
    """
    Load one JSON database file and materialize a temporary reference FASTA for BLAST.
    The JSON must contain a top-level records list.
    """
    log(f"Loading DB JSON: {db_path}", verbose)
    payload = json.loads(db_path.read_text())
    records = payload.get("records")
    if not isinstance(records, list) or not records:
        die("Invalid DB JSON: expected a non-empty top-level 'records' list.")

    ref_fasta = work_dir / "gp60_refs_from_json.fasta"
    meta: Dict[str, RefMeta] = {}

    with ref_fasta.open("w") as out:
        for rec in records:
            ref_id = str(rec.get("ref_id", "")).strip()
            seq = clean_seq(rec.get("sequence", ""))
            if not ref_id or not seq:
                continue

            species = str(rec.get("species", "NA"))
            gene = str(rec.get("gene", "gp60"))
            family = str(rec.get("gp60_family", "NA"))
            subtype = str(rec.get("gp60_subtype", "NA"))
            accession = str(rec.get("accession", "NA"))

            out.write(
                f">{ref_id} species={species} gene={gene} "
                f"family={family} subtype={subtype} accession={accession}\n{wrap(seq)}\n"
            )

            try:
                length_bp = int(rec.get("length_bp", len(seq)))
            except Exception:
                length_bp = len(seq)

            meta[ref_id] = RefMeta(
                ref_id=ref_id,
                species=species,
                species_raw=str(rec.get("species_raw", species)),
                gene=gene,
                gp60_family=family,
                gp60_variant=str(rec.get("gp60_variant", "NA")),
                gp60_subtype=subtype,
                accession=accession,
                length_bp=length_bp,
                original_header=str(rec.get("original_header", "NA")),
            )

    if not meta:
        die("No valid records were loaded from DB JSON.")
    log(f"Loaded {len(meta)} gp60 reference records from DB.", verbose)
    log(f"Materialized temporary BLAST reference FASTA: {ref_fasta}", verbose)
    return ref_fasta, meta


def find_fasta_files(input_dir: Path) -> List[Path]:
    exts = {".fa", ".fna", ".fasta", ".fas"}
    return [p for p in sorted(input_dir.iterdir()) if p.is_file() and p.suffix.lower() in exts]


def sanitize_id(name: str) -> str:
    name = Path(name).name
    for suffix in [".fasta", ".fna", ".fa", ".fas"]:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", name)
    return name.strip("_") or "sample"


def make_blast_db(reference_fasta: Path, db_prefix: Path, verbose: bool = False) -> None:
    require_executable("makeblastdb")
    log("Building temporary BLAST database.", verbose)
    cmd = ["makeblastdb", "-in", str(reference_fasta), "-dbtype", "nucl", "-parse_seqids", "-out", str(db_prefix)]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def run_blast_query_against_db(query_fasta: Path, db_prefix: Path, out_tsv: Path, threads: int, verbose: bool = False) -> None:
    require_executable("blastn")
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq"
    log(f"Running BLAST query mode: {query_fasta}", verbose)
    cmd = [
        "blastn", "-query", str(query_fasta), "-db", str(db_prefix),
        "-outfmt", outfmt, "-max_target_seqs", "50",
        "-num_threads", str(threads), "-dust", "no", "-soft_masking", "false",
        "-out", str(out_tsv),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def run_blast_refs_against_genome(ref_fasta: Path, genome_fasta: Path, out_tsv: Path, threads: int, verbose: bool = False) -> None:
    require_executable("blastn")
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq"
    log(f"Running BLAST against genome: {genome_fasta.name}", verbose)
    cmd = [
        "blastn", "-query", str(ref_fasta), "-subject", str(genome_fasta),
        "-outfmt", outfmt, "-max_target_seqs", "50",
        "-num_threads", str(threads), "-dust", "no", "-soft_masking", "false",
        "-out", str(out_tsv),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def parse_blast_tsv(path: Path) -> List[BlastHit]:
    hits: List[BlastHit] = []
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 16:
                continue
            hits.append(BlastHit(
                qseqid=parts[0], sseqid=parts[1], pident=float(parts[2]),
                length=int(parts[3]), mismatch=int(parts[4]), gapopen=int(parts[5]),
                qstart=int(parts[6]), qend=int(parts[7]), sstart=int(parts[8]), send=int(parts[9]),
                evalue=float(parts[10]), bitscore=float(parts[11]),
                qlen=int(parts[12]), slen=int(parts[13]), qseq=parts[14], sseq=parts[15],
            ))
    return hits


def hit_rank_score(hit: BlastHit) -> float:
    return (
        hit.pident * 1.0
        + hit.qcov * 1.2
        + min(hit.length, hit.qlen) * 0.01
        + hit.bitscore * 0.001
        - hit.gapopen * 0.2
    )


def pick_best_hits(hits: List[BlastHit], max_delta: float = 1.0) -> Tuple[Optional[BlastHit], List[BlastHit]]:
    if not hits:
        return None, []
    ranked = sorted(hits, key=hit_rank_score, reverse=True)
    best = ranked[0]
    best_score = hit_rank_score(best)
    close = [h for h in ranked if best_score - hit_rank_score(h) <= max_delta]
    return best, close


def load_contigs(genome_fasta: Path) -> Dict[str, str]:
    return {rec.id: rec.seq for rec in parse_fasta(genome_fasta)}


def extract_region_from_hit(hit: BlastHit, contigs: Dict[str, str], pad: int = 0) -> str:
    if hit.sseqid not in contigs:
        return ""
    contig = contigs[hit.sseqid]
    start = min(hit.sstart, hit.send)
    end = max(hit.sstart, hit.send)
    start = max(1, start - pad)
    end = min(len(contig), end + pad)
    seq = contig[start - 1:end]
    if hit.strand == "minus":
        seq = revcomp(seq)
    return clean_seq(seq)


def classify_status(best: Optional[BlastHit], close_hits: List[BlastHit], min_qcov: float, min_pident: float) -> Tuple[str, str, str]:
    if best is None:
        return "not_found", "none", "no gp60 hit"

    notes = []
    if best.qcov < min_qcov:
        status = "partial"
        confidence = "low" if best.qcov < 50 else "medium"
        notes.append(f"reference coverage below threshold: {best.qcov:.1f}%")
    elif best.pident < min_pident:
        status = "low_identity"
        confidence = "low"
        notes.append(f"identity below threshold: {best.pident:.2f}%")
    else:
        status = "complete"
        confidence = "high"

    if len(close_hits) > 1:
        top_ids = [h.qseqid for h in close_hits[:5]]
        notes.append("close competing references: " + ",".join(top_ids))

    return status, confidence, "; ".join(notes) if notes else "-"


def collapse_close_hit_labels(close_hits: List[BlastHit], meta: Dict[str, RefMeta]) -> Tuple[str, str, str, str]:
    species = sorted({meta[h.qseqid].species for h in close_hits if h.qseqid in meta})
    families = sorted({meta[h.qseqid].gp60_family for h in close_hits if h.qseqid in meta})
    subtypes = sorted({meta[h.qseqid].gp60_subtype for h in close_hits if h.qseqid in meta})
    accessions = sorted({meta[h.qseqid].accession for h in close_hits if h.qseqid in meta})
    return "|".join(species), "|".join(families), "|".join(subtypes), "|".join(accessions)


def type_assemblies(
    input_dir: Path,
    db_path: Path,
    output_dir: Path,
    threads: int,
    min_qcov: float,
    min_pident: float,
    ambiguity_delta: float,
    verbose: bool = False,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = output_dir / "blast_results"
    blast_dir.mkdir(exist_ok=True)

    genomes = find_fasta_files(input_dir)
    if not genomes:
        die(f"No FASTA files found in {input_dir}")
    log(f"Found {len(genomes)} genome FASTA file(s) in {input_dir}.", verbose)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        ref_fasta, meta = load_single_json_db(db_path, tmp_path, verbose=verbose)

        results_path = output_dir / "gp60_typing_results.tsv"
        extracted_path = output_dir / "extracted_gp60_sequences.fasta"
        failed_path = output_dir / "failed_or_low_confidence_samples.tsv"

        extracted_records: List[Tuple[str, str, Optional[str]]] = []
        result_rows = []

        for idx, genome in enumerate(genomes, start=1):
            sample_id = sanitize_id(genome.name)
            log(f"[{idx}/{len(genomes)}] Processing sample: {sample_id}", verbose)
            blast_tsv = blast_dir / f"{sample_id}.gp60_refs_vs_genome.tsv"
            run_blast_refs_against_genome(ref_fasta, genome, blast_tsv, threads, verbose=verbose)
            hits = parse_blast_tsv(blast_tsv)
            log(f"[{idx}/{len(genomes)}] BLAST hits found: {len(hits)}", verbose)
            best, close = pick_best_hits(hits, max_delta=ambiguity_delta)

            if best is None:
                log(f"[{idx}/{len(genomes)}] No gp60 hit found for {sample_id}.", verbose)
                result_rows.append({
                    "sample_id": sample_id, "input_file": str(genome), "likely_species": "NA",
                    "gp60_family": "NA", "gp60_subtype": "NA", "best_reference": "NA",
                    "accession": "NA", "identity_pct": "0.000", "reference_coverage_pct": "0.000",
                    "alignment_length_bp": "0", "contig": "NA", "contig_start": "NA",
                    "contig_end": "NA", "strand": "NA", "status": "not_found",
                    "confidence": "none", "notes": "no gp60 hit",
                })
                continue

            best_meta = meta.get(best.qseqid)
            if best_meta is None:
                die(f"Best reference {best.qseqid} not found in DB metadata.")
            log(
                f"[{idx}/{len(genomes)}] Best hit: {best.qseqid} "
                f"species={best_meta.species} family={best_meta.gp60_family} "
                f"subtype={best_meta.gp60_subtype} identity={best.pident:.3f}% "
                f"qcov={best.qcov:.3f}%",
                verbose,
            )

            status, confidence, notes = classify_status(best, close, min_qcov, min_pident)
            close_species, close_families, close_subtypes, close_accessions = collapse_close_hit_labels(close, meta)

            if len(set(close_subtypes.split("|"))) > 1 and status == "complete":
                confidence = "medium"
                notes = (notes + "; " if notes != "-" else "") + f"close subtype-level tie: {close_subtypes}"
                reported_subtype = "ambiguous:" + close_subtypes
            else:
                reported_subtype = best_meta.gp60_subtype

            if len(set(close_families.split("|"))) > 1 and status == "complete":
                status = "ambiguous"
                confidence = "low"
                notes = (notes + "; " if notes != "-" else "") + f"close family-level tie: {close_families}"

            contigs = load_contigs(genome)
            extracted_seq = extract_region_from_hit(best, contigs, pad=0)
            contig_start = min(best.sstart, best.send)
            contig_end = max(best.sstart, best.send)

            if extracted_seq:
                ex_id = f"{sample_id}|best_ref={best.qseqid}|species={best_meta.species}|family={best_meta.gp60_family}|subtype={reported_subtype}"
                ex_desc = f"identity={best.pident:.3f} qcov={best.qcov:.3f} status={status}"
                extracted_records.append((ex_id, extracted_seq, ex_desc))

            result_rows.append({
                "sample_id": sample_id, "input_file": str(genome),
                "likely_species": best_meta.species,
                "gp60_family": best_meta.gp60_family if status != "ambiguous" else "ambiguous:" + close_families,
                "gp60_subtype": reported_subtype,
                "best_reference": best.qseqid,
                "accession": best_meta.accession,
                "identity_pct": f"{best.pident:.3f}",
                "reference_coverage_pct": f"{best.qcov:.3f}",
                "alignment_length_bp": str(best.length),
                "contig": best.sseqid,
                "contig_start": str(contig_start),
                "contig_end": str(contig_end),
                "strand": best.strand,
                "status": status,
                "confidence": confidence,
                "notes": notes,
            })

        fieldnames = [
            "sample_id", "input_file", "likely_species", "gp60_family", "gp60_subtype",
            "best_reference", "accession", "identity_pct", "reference_coverage_pct",
            "alignment_length_bp", "contig", "contig_start", "contig_end", "strand",
            "status", "confidence", "notes",
        ]
        log("Writing output tables and extracted gp60 FASTA.", verbose)
        with results_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(result_rows)

        write_fasta(extracted_records, extracted_path)

        with failed_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in result_rows:
                if row["confidence"] != "high":
                    writer.writerow(row)

    print(f"Wrote: {results_path}")
    print(f"Wrote: {extracted_path}")
    print(f"Wrote: {failed_path}")
    print(f"Wrote BLAST files in: {blast_dir}")


def type_extracted_gp60(
    query_fasta: Path,
    db_path: Path,
    output_dir: Path,
    threads: int,
    min_qcov: float,
    min_pident: float,
    ambiguity_delta: float,
    verbose: bool = False,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = output_dir / "blast_results"
    blast_dir.mkdir(exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        ref_fasta, meta = load_single_json_db(db_path, tmp_path, verbose=verbose)
        db_prefix = tmp_path / "gp60_ref_db"
        make_blast_db(ref_fasta, db_prefix, verbose=verbose)
        blast_tsv = blast_dir / "queries_vs_gp60_refs.tsv"
        run_blast_query_against_db(query_fasta, db_prefix, blast_tsv, threads, verbose=verbose)

    hits = parse_blast_tsv(blast_tsv)
    by_query: Dict[str, List[BlastHit]] = defaultdict(list)
    for h in hits:
        by_query[h.qseqid].append(h)

    query_ids = [rec.id for rec in parse_fasta(query_fasta)]
    log(f"Found {len(query_ids)} query sequence(s) in {query_fasta}.", verbose)
    rows = []
    for idx, qid in enumerate(query_ids, start=1):
        log(f"[{idx}/{len(query_ids)}] Typing query: {qid}", verbose)
        qhits = by_query.get(qid, [])
        best, close = pick_best_hits(qhits, max_delta=ambiguity_delta)

        if best is None:
            log(f"[{idx}/{len(query_ids)}] No gp60 reference hit found.", verbose)
            rows.append({
                "sample_id": qid, "likely_species": "NA", "gp60_family": "NA",
                "gp60_subtype": "NA", "best_reference": "NA", "accession": "NA",
                "identity_pct": "0.000", "query_coverage_pct": "0.000",
                "reference_coverage_pct": "0.000", "alignment_length_bp": "0",
                "status": "not_found", "confidence": "none", "notes": "no gp60 reference hit",
            })
            continue

        best_meta = meta.get(best.sseqid)
        if best_meta is None:
            die(f"Best reference {best.sseqid} not found in DB metadata.")
        log(
            f"[{idx}/{len(query_ids)}] Best hit: {best.sseqid} "
            f"species={best_meta.species} family={best_meta.gp60_family} "
            f"subtype={best_meta.gp60_subtype} identity={best.pident:.3f}% "
            f"query_cov={best.qcov:.3f}% ref_cov={best.scov:.3f}%",
            verbose,
        )

        status = "complete"
        confidence = "high"
        notes = []
        if best.scov < min_qcov:
            status = "partial"
            confidence = "low" if best.scov < 50 else "medium"
            notes.append(f"reference coverage below threshold: {best.scov:.1f}%")
        if best.pident < min_pident:
            status = "low_identity"
            confidence = "low"
            notes.append(f"identity below threshold: {best.pident:.2f}%")

        close_species = sorted({meta[h.sseqid].species for h in close if h.sseqid in meta})
        close_families = sorted({meta[h.sseqid].gp60_family for h in close if h.sseqid in meta})
        close_subtypes = sorted({meta[h.sseqid].gp60_subtype for h in close if h.sseqid in meta})
        reported_subtype = best_meta.gp60_subtype

        if len(set(close_subtypes)) > 1 and status == "complete":
            confidence = "medium"
            reported_subtype = "ambiguous:" + "|".join(close_subtypes)
            notes.append("close subtype-level tie: " + "|".join(close_subtypes))
        if len(set(close_families)) > 1 and status == "complete":
            status = "ambiguous"
            confidence = "low"
            notes.append("close family-level tie: " + "|".join(close_families))

        rows.append({
            "sample_id": qid,
            "likely_species": best_meta.species if len(close_species) <= 1 else "ambiguous:" + "|".join(close_species),
            "gp60_family": best_meta.gp60_family if len(close_families) <= 1 else "ambiguous:" + "|".join(close_families),
            "gp60_subtype": reported_subtype,
            "best_reference": best.sseqid,
            "accession": best_meta.accession,
            "identity_pct": f"{best.pident:.3f}",
            "query_coverage_pct": f"{best.qcov:.3f}",
            "reference_coverage_pct": f"{best.scov:.3f}",
            "alignment_length_bp": str(best.length),
            "status": status,
            "confidence": confidence,
            "notes": "; ".join(notes) if notes else "-",
        })

    out_path = output_dir / "gp60_typing_results.tsv"
    fieldnames = [
        "sample_id", "likely_species", "gp60_family", "gp60_subtype",
        "best_reference", "accession", "identity_pct", "query_coverage_pct",
        "reference_coverage_pct", "alignment_length_bp", "status", "confidence", "notes",
    ]
    log("Writing query-mode typing table.", verbose)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote: {out_path}")
    print(f"Wrote BLAST file in: {blast_dir}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Cryptosporidium gp60 typer using a single JSON database file with embedded references and metadata."
    )
    sub = p.add_subparsers(dest="mode", required=True)

    g = sub.add_parser("genome", help="Type gp60 from genome/assembly FASTA files.")
    g.add_argument("-d", "--input-dir", required=True, type=Path, help="Directory with genome FASTA files.")
    g.add_argument("--db", required=True, type=Path, help="Single gp60 JSON DB file.")
    g.add_argument("-o", "--output-dir", required=True, type=Path, help="Output directory.")
    g.add_argument("-t", "--threads", default=1, type=int, help="BLAST threads. Default: 1.")
    g.add_argument("--min-qcov", default=85.0, type=float, help="Minimum reference coverage percent for complete call. Default: 85.")
    g.add_argument("--min-pident", default=95.0, type=float, help="Minimum percent identity for confident call. Default: 95.")
    g.add_argument("--ambiguity-delta", default=1.0, type=float, help="Score delta for close competing references. Default: 1.0.")
    g.add_argument("-v", "--verbose", action="store_true", help="Print progress messages to stderr.")

    q = sub.add_parser("query", help="Type already extracted gp60 FASTA sequences against the DB.")
    q.add_argument("-i", "--input-fasta", required=True, type=Path, help="FASTA with extracted gp60/query sequences.")
    q.add_argument("--db", required=True, type=Path, help="Single gp60 JSON DB file.")
    q.add_argument("-o", "--output-dir", required=True, type=Path, help="Output directory.")
    q.add_argument("-t", "--threads", default=1, type=int, help="BLAST threads. Default: 1.")
    q.add_argument("--min-qcov", default=85.0, type=float, help="Minimum reference coverage percent for complete call. Default: 85.")
    q.add_argument("--min-pident", default=95.0, type=float, help="Minimum percent identity for confident call. Default: 95.")
    q.add_argument("--ambiguity-delta", default=1.0, type=float, help="Score delta for close competing references. Default: 1.0.")
    q.add_argument("-v", "--verbose", action="store_true", help="Print progress messages to stderr.")

    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.mode == "genome":
        type_assemblies(
            input_dir=args.input_dir,
            db_path=args.db,
            output_dir=args.output_dir,
            threads=args.threads,
            min_qcov=args.min_qcov,
            min_pident=args.min_pident,
            ambiguity_delta=args.ambiguity_delta,
            verbose=args.verbose,
        )
    elif args.mode == "query":
        type_extracted_gp60(
            query_fasta=args.input_fasta,
            db_path=args.db,
            output_dir=args.output_dir,
            threads=args.threads,
            min_qcov=args.min_qcov,
            min_pident=args.min_pident,
            ambiguity_delta=args.ambiguity_delta,
            verbose=args.verbose,
        )


if __name__ == "__main__":
    main()
