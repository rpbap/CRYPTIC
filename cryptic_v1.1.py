#!/usr/bin/env python3
# CRYPTIC v1.1.8 draft: Cryptosporidium gp60 Typing and Identification Classifier
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
        print(f"[CRYPTIC] {msg}", file=sys.stderr, flush=True)


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



def strip_fastq_extensions(path: Path) -> str:
    name = path.name
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def is_fastq_file(path: Path) -> bool:
    return path.name.lower().endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))


def detect_illumina_pair_token(path: Path) -> Optional[Tuple[str, str]]:
    base = strip_fastq_extensions(path)
    patterns = [
        r"^(?P<sample>.+?)(?P<sep>[_\.-])R(?P<read>[12])(?:_001)?$",
        r"^(?P<sample>.+?)(?P<sep>[_\.-])(?P<read>[12])(?:_001)?$",
    ]
    for pat in patterns:
        m = re.match(pat, base, flags=re.IGNORECASE)
        if m:
            return sanitize_id(m.group("sample")), m.group("read")
    return None


def discover_illumina_pairs(input_dir: Path) -> Tuple[List[Tuple[str, Path, Path]], List[Dict[str, str]]]:
    fastqs = [p for p in sorted(input_dir.iterdir()) if p.is_file() and is_fastq_file(p)]
    grouped: Dict[str, Dict[str, Path]] = defaultdict(dict)
    report_rows: List[Dict[str, str]] = []

    for fq in fastqs:
        parsed = detect_illumina_pair_token(fq)
        if parsed is None:
            report_rows.append({"sample_id": sanitize_id(strip_fastq_extensions(fq)), "R1": "NA", "R2": "NA", "status": "unrecognized_name"})
            continue
        sample, read = parsed
        if read in grouped[sample]:
            report_rows.append({"sample_id": sample, "R1": "NA", "R2": "NA", "status": f"duplicate_R{read}:{fq.name}"})
            continue
        grouped[sample][read] = fq

    pairs: List[Tuple[str, Path, Path]] = []
    for sample in sorted(grouped):
        r1 = grouped[sample].get("1")
        r2 = grouped[sample].get("2")
        if r1 and r2:
            pairs.append((sample, r1, r2))
            report_rows.append({"sample_id": sample, "R1": str(r1), "R2": str(r2), "status": "paired"})
        elif r1 and not r2:
            report_rows.append({"sample_id": sample, "R1": str(r1), "R2": "NA", "status": "missing_R2"})
        elif r2 and not r1:
            report_rows.append({"sample_id": sample, "R1": "NA", "R2": str(r2), "status": "missing_R1"})
    return pairs, report_rows


def discover_ont_reads(input_dir: Path) -> List[Tuple[str, Path]]:
    fastqs = [p for p in sorted(input_dir.iterdir()) if p.is_file() and is_fastq_file(p)]
    return [(sanitize_id(strip_fastq_extensions(p)), p) for p in fastqs]


def write_pairing_report(rows: List[Dict[str, str]], output_path: Path) -> None:
    fields = ["sample_id", "R1", "R2", "status"]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_json_db_to_ref_fasta_and_meta(db_path: Path, work_dir: Path, verbose: bool = False) -> Tuple[Path, Dict[str, RefMeta], Dict[str, str]]:
    ref_fasta, meta = load_single_json_db(db_path, work_dir, verbose=verbose)
    ref_sequences = {rec.id: rec.seq for rec in parse_fasta(ref_fasta)}
    return ref_fasta, meta, ref_sequences


def run_minimap2_to_sorted_bam(ref_fasta: Path, bam_path: Path, preset: str, threads: int, read1: Path, read2: Optional[Path], verbose: bool = False) -> None:
    require_executable("minimap2")
    require_executable("samtools")
    log(f"Mapping reads with minimap2 preset '{preset}'.", verbose)
    if read2 is None:
        minimap_cmd = ["minimap2", "-ax", preset, "-t", str(threads), str(ref_fasta), str(read1)]
    else:
        minimap_cmd = ["minimap2", "-ax", preset, "-t", str(threads), str(ref_fasta), str(read1), str(read2)]
    sort_cmd = ["samtools", "sort", "-@", str(max(1, threads - 1)), "-o", str(bam_path), "-"]
    p1 = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p1.stdout:
        p1.stdout.close()
    _, sort_err = p2.communicate()
    map_err = p1.stderr.read() if p1.stderr else b""
    ret1 = p1.wait()
    ret2 = p2.returncode
    if ret1 != 0:
        die("minimap2 failed:\n" + map_err.decode(errors="replace"))
    if ret2 != 0:
        die("samtools sort failed:\n" + sort_err.decode(errors="replace"))
    subprocess.run(["samtools", "index", str(bam_path)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def parse_samtools_idxstats(bam_path: Path) -> Dict[str, Dict[str, int]]:
    require_executable("samtools")
    result = subprocess.run(["samtools", "idxstats", str(bam_path)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stats: Dict[str, Dict[str, int]] = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        ref, length, mapped, unmapped = line.split("\t")[:4]
        if ref == "*":
            continue
        stats[ref] = {"length": int(length), "mapped": int(mapped), "unmapped": int(unmapped)}
    return stats


def count_fastq_reads(path: Optional[Path]) -> Optional[int]:
    if path is None:
        return None
    import gzip
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    try:
        with opener(path, "rt", errors="replace") as handle:
            return sum(1 for _ in handle) // 4
    except Exception:
        return None


def pileup_bases_to_counts_with_ref(bases: str, ref_base: str) -> Dict[str, int]:
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    i = 0
    ref_base = ref_base.upper()
    while i < len(bases):
        ch = bases[i]
        if ch == "^":
            i += 2
            continue
        if ch == "$":
            i += 1
            continue
        if ch in "+-":
            i += 1
            num = []
            while i < len(bases) and bases[i].isdigit():
                num.append(bases[i])
                i += 1
            indel_len = int("".join(num)) if num else 0
            i += indel_len
            continue
        if ch in ".,":
            if ref_base in counts:
                counts[ref_base] += 1
            i += 1
            continue
        base = ch.upper()
        if base in counts:
            counts[base] += 1
        i += 1
    return counts


def build_reference_guided_consensus(bam_path: Path, ref_fasta: Path, ref_id: str, ref_seq: str, min_depth: int, min_base_fraction: float, verbose: bool = False) -> Tuple[str, D
ict[str, float]]:
    require_executable("samtools")
    cmd = ["samtools", "mpileup", "-aa", "-A", "-d", "1000000", "-f", str(ref_fasta), "-r", ref_id, str(bam_path)]
    log(f"Building reference-guided consensus for {ref_id}.", verbose)
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    consensus = list("N" * len(ref_seq))
    depths: List[int] = []
    covered_positions = 0
    ambiguous_positions = 0
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        pos = int(parts[1])
        ref_base = parts[2].upper()
        depth = int(parts[3])
        bases = parts[4]
        if pos < 1 or pos > len(consensus):
            continue
        depths.append(depth)
        if depth >= min_depth:
            covered_positions += 1
            counts = pileup_bases_to_counts_with_ref(bases, ref_base)
            total = sum(counts.values())
            if total == 0:
                ambiguous_positions += 1
                continue
            best_base, best_count = max(counts.items(), key=lambda x: x[1])
            if best_count / total >= min_base_fraction:
                consensus[pos - 1] = best_base
            else:
                ambiguous_positions += 1
    non_n = sum(1 for b in consensus if b in "ACGT")
    metrics = {
        "mean_depth": (sum(depths) / len(depths)) if depths else 0.0,
        "breadth_coverage_pct": (100.0 * non_n / len(consensus)) if consensus else 0.0,
        "consensus_n_count": sum(1 for b in consensus if b == "N"),
        "covered_positions_min_depth": covered_positions,
        "ambiguous_positions": ambiguous_positions,
    }
    return "".join(consensus), metrics


def classify_consensus_sequence(sample_id: str, consensus_seq: str, db_path: Path, threads: int, min_qcov: float, min_pident: float, ambiguity_delta: float) -> Dict[str, str]:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        query_fasta = tmp_path / f"{sample_id}.consensus.fasta"
        write_fasta([(sample_id, consensus_seq, "read_derived_gp60_consensus")], query_fasta)
        ref_fasta, meta = load_single_json_db(db_path, tmp_path, verbose=False)
        db_prefix = tmp_path / "gp60_ref_db"
        make_blast_db(ref_fasta, db_prefix, verbose=False)
        blast_tsv = tmp_path / "consensus_vs_gp60_refs.tsv"
        run_blast_query_against_db(query_fasta, db_prefix, blast_tsv, threads, verbose=False)
        hits = parse_blast_tsv(blast_tsv)
        best, close = pick_best_hits(hits, max_delta=ambiguity_delta)
        if best is None:
            return {"likely_species": "NA", "gp60_family": "NA", "gp60_subtype": "NA", "best_reference": "NA", "accession": "NA", "identity_pct": "0.000", "query_coverage_pct": 
"0.000", "reference_coverage_pct": "0.000", "alignment_length_bp": "0", "status": "not_found", "confidence": "none", "notes": "consensus had no gp60 reference hit"}
        best_meta = meta.get(best.sseqid)
        if best_meta is None:
            die(f"Best reference {best.sseqid} not found in DB metadata.")
        status = "complete"
        confidence = "high"
        notes: List[str] = []
        if best.scov < min_qcov:
            status = "partial"
            confidence = "low" if best.scov < 50 else "medium"
            notes.append(f"reference coverage below threshold: {best.scov:.1f}%")
        if best.pident < min_pident:
            status = "low_identity"
            confidence = "low"
            notes.append(f"identity below threshold: {best.pident:.2f}")
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
        return {
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
        }


def process_read_sample(
    sample_id: str,
    read_type: str,
    read1: Path,
    read2: Optional[Path],
    db_path: Path,
    ref_fasta: Path,
    meta: Dict[str, RefMeta],
    ref_sequences: Dict[str, str],
    output_dir: Path,
    threads: int,
    min_qcov: float,
    min_pident: float,
    ambiguity_delta: float,
    min_depth: int,
    min_base_fraction: float,
    min_mapped_reads: int,
    candidate_limit: int = 25,
    min_candidate_breadth: float = 20.0,
    verbose: bool = False,
) -> Tuple[Dict[str, str], Optional[Tuple[str, str, str]]]:
    """
    Process one reads sample.

    v1.1.1 fix:
    Do NOT choose the read-derived reference by raw mapped-read count alone.
    gp60 references share conserved regions, and multi-mapped reads can make an
    unrelated reference appear to have many reads over a tiny fraction of the gene.

    Instead:
      1. Map reads to all gp60 references.
      2. Take top references by mapped reads as candidates.
      3. Build a simple consensus for each candidate.
      4. Classify each candidate consensus back against the DB.
      5. Select the best candidate by breadth, reference coverage, identity, and depth.
    """
    bam_dir = output_dir / "bam"
    bam_dir.mkdir(exist_ok=True)
    bam_path = bam_dir / f"{sample_id}.gp60.sorted.bam"

    preset = "sr" if read_type == "illumina" else "map-ont"
    run_minimap2_to_sorted_bam(ref_fasta, bam_path, preset, threads, read1, read2, verbose=verbose)

    idxstats = parse_samtools_idxstats(bam_path)
    total_mapped = sum(v["mapped"] for v in idxstats.values())

    candidate_refs = [
        (ref_id, payload["mapped"])
        for ref_id, payload in idxstats.items()
        if payload["mapped"] >= min_mapped_reads and ref_id in ref_sequences
    ]
    candidate_refs.sort(key=lambda x: x[1], reverse=True)
    candidate_refs = candidate_refs[:candidate_limit]

    total_reads_1 = count_fastq_reads(read1)
    total_reads_2 = count_fastq_reads(read2) if read2 else None
    total_reads = "NA"
    if total_reads_1 is not None:
        total_reads = str(total_reads_1 + (total_reads_2 or 0))

    if not candidate_refs:
        row = {
            "sample_id": sample_id,
            "read_type": read_type,
            "R1_or_reads": str(read1),
            "R2": str(read2) if read2 else "NA",
            "total_reads": total_reads,
            "mapped_gp60_reads": str(total_mapped),
            "diagnostic_top_raw_read_reference": "NA",
            "diagnostic_top_raw_read_reference_reads": "0",
            "mean_gp60_depth": "0.000",
            "gp60_breadth_coverage_pct": "0.000",
            "consensus_n_count": "NA",
            "likely_species": "NA",
            "gp60_family": "NA",
            "gp60_subtype": "NA",
            "best_reference": "NA",
            "accession": "NA",
            "identity_pct": "0.000",
            "query_coverage_pct": "0.000",
            "reference_coverage_pct": "0.000",
            "alignment_length_bp": "0",
            "status": "not_found",
            "confidence": "none",
            "notes": f"no reference had at least {min_mapped_reads} mapped reads",
        }
        return row, None

    log(
        f"Evaluating {len(candidate_refs)} candidate reference(s) by consensus breadth/identity, "
        f"not raw read count alone.",
        verbose,
    )

    candidate_rows = []
    candidate_consensus_records = []

    for rank, (candidate_ref, candidate_mapped) in enumerate(candidate_refs, start=1):
        ref_seq = ref_sequences.get(candidate_ref)
        if not ref_seq:
            continue

        log(
            f"Candidate {rank}/{len(candidate_refs)}: {candidate_ref} "
            f"({candidate_mapped} mapped reads by idxstats)",
            verbose,
        )

        consensus, metrics = build_reference_guided_consensus(
            bam_path=bam_path,
            ref_fasta=ref_fasta,
            ref_id=candidate_ref,
            ref_seq=ref_seq,
            min_depth=min_depth,
            min_base_fraction=min_base_fraction,
            verbose=False,
        )

        class_row = classify_consensus_sequence(
            sample_id=sample_id,
            consensus_seq=consensus,
            db_path=db_path,
            threads=threads,
            min_qcov=min_qcov,
            min_pident=min_pident,
            ambiguity_delta=ambiguity_delta,
        )

        try:
            identity = float(class_row.get("identity_pct", 0))
        except Exception:
            identity = 0.0
        try:
            ref_cov = float(class_row.get("reference_coverage_pct", 0))
        except Exception:
            ref_cov = 0.0

        breadth = float(metrics["breadth_coverage_pct"])
        mean_depth = float(metrics["mean_depth"])

        # Selection score emphasizes breadth and ref coverage first.
        # This avoids selecting a conserved 6-7% fragment with many reads.
        candidate_score = (
            breadth * 2.0
            + ref_cov * 2.0
            + identity * 0.5
            + min(mean_depth, 200.0) * 0.05
            + min(candidate_mapped, 10000) * 0.0001
        )

        candidate_rows.append({
            "candidate_ref": candidate_ref,
            "candidate_mapped": candidate_mapped,
            "consensus": consensus,
            "metrics": metrics,
            "class_row": class_row,
            "selection_score": candidate_score,
            "identity": identity,
            "ref_cov": ref_cov,
            "breadth": breadth,
            "mean_depth": mean_depth,
        })

    if not candidate_rows:
        row = {
            "sample_id": sample_id,
            "read_type": read_type,
            "R1_or_reads": str(read1),
            "R2": str(read2) if read2 else "NA",
            "total_reads": total_reads,
            "mapped_gp60_reads": str(total_mapped),
            "diagnostic_top_raw_read_reference": candidate_refs[0][0] if candidate_refs else "NA",
            "diagnostic_top_raw_read_reference_reads": str(candidate_refs[0][1]) if candidate_refs else "0",
            "mean_gp60_depth": "0.000",
            "gp60_breadth_coverage_pct": "0.000",
            "consensus_n_count": "NA",
            "likely_species": "NA",
            "gp60_family": "NA",
            "gp60_subtype": "NA",
            "best_reference": "NA",
            "accession": "NA",
            "identity_pct": "0.000",
            "query_coverage_pct": "0.000",
            "reference_coverage_pct": "0.000",
            "alignment_length_bp": "0",
            "status": "not_found",
            "confidence": "none",
            "notes": "candidate consensus generation failed",
        }
        return row, None

    # Prefer candidates with enough breadth. If none pass, still report the best,
    # but keep it low-confidence/partial.
    passing = [c for c in candidate_rows if c["breadth"] >= min_candidate_breadth]
    selection_pool = passing if passing else candidate_rows
    selected = sorted(selection_pool, key=lambda c: c["selection_score"], reverse=True)[0]

    best_ref_by_reads, best_ref_mapped_reads = candidate_refs[0]
    selected_ref = selected["candidate_ref"]
    selected_mapped = selected["candidate_mapped"]
    metrics = selected["metrics"]
    class_row = selected["class_row"]
    consensus = selected["consensus"]

    log(
        f"Selected candidate: {selected_ref} "
        f"mapped_reads={selected_mapped} breadth={selected['breadth']:.3f}% "
        f"ref_cov={selected['ref_cov']:.3f}% identity={selected['identity']:.3f}%",
        verbose,
    )

    notes = class_row["notes"]

    if selected["breadth"] < min_qcov:
        extra = f"consensus breadth below threshold: {selected['breadth']:.1f}%"
        notes = extra if notes == "-" else notes + "; " + extra
        if class_row["confidence"] == "high":
            class_row["confidence"] = "medium"
        if class_row["status"] == "complete":
            class_row["status"] = "partial"

    if selected["breadth"] < min_candidate_breadth:
        extra = (
            f"very low breadth across all candidates; best raw-read reference was "
            f"{best_ref_by_reads} with {best_ref_mapped_reads} mapped reads"
        )
        notes = extra if notes == "-" else notes + "; " + extra
        class_row["confidence"] = "low"
        if class_row["status"] == "complete":
            class_row["status"] = "partial"

    row = {
        "sample_id": sample_id,
        "read_type": read_type,
        "R1_or_reads": str(read1),
        "R2": str(read2) if read2 else "NA",
        "total_reads": total_reads,
        "mapped_gp60_reads": str(total_mapped),
        "diagnostic_top_raw_read_reference": best_ref_by_reads,
        "diagnostic_top_raw_read_reference_reads": str(best_ref_mapped_reads),
        "mean_gp60_depth": f"{metrics['mean_depth']:.3f}",
        "gp60_breadth_coverage_pct": f"{metrics['breadth_coverage_pct']:.3f}",
        "consensus_n_count": str(int(metrics["consensus_n_count"])),
        **class_row,
        "notes": notes,
    }

    consensus_record = (
        f"{sample_id}|read_type={read_type}|selected_ref={selected_ref}|top_raw_read_ref={best_ref_by_reads}",
        consensus,
        f"mean_depth={metrics['mean_depth']:.3f} breadth={metrics['breadth_coverage_pct']:.3f} mapped_reads={selected_mapped}",
    )

    # Write candidate diagnostics for this sample.
    cand_dir = output_dir / "read_candidate_diagnostics"
    cand_dir.mkdir(exist_ok=True)
    cand_path = cand_dir / f"{sample_id}.candidate_refs.tsv"
    cand_fields = [
        "sample_id", "candidate_ref", "candidate_mapped_reads", "selection_score",
        "mean_depth", "breadth_coverage_pct", "consensus_n_count",
        "classified_species", "classified_family", "classified_subtype",
        "classified_best_reference", "identity_pct", "reference_coverage_pct",
        "status", "confidence", "notes",
    ]
    with cand_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cand_fields, delimiter="\t")
        writer.writeheader()
        for c in sorted(candidate_rows, key=lambda x: x["selection_score"], reverse=True):
            cr = c["class_row"]
            writer.writerow({
                "sample_id": sample_id,
                "candidate_ref": c["candidate_ref"],
                "candidate_mapped_reads": c["candidate_mapped"],
                "selection_score": f"{c['selection_score']:.6f}",
                "mean_depth": f"{c['mean_depth']:.3f}",
                "breadth_coverage_pct": f"{c['breadth']:.3f}",
                "consensus_n_count": str(int(c["metrics"]["consensus_n_count"])),
                "classified_species": cr.get("likely_species", "NA"),
                "classified_family": cr.get("gp60_family", "NA"),
                "classified_subtype": cr.get("gp60_subtype", "NA"),
                "classified_best_reference": cr.get("best_reference", "NA"),
                "identity_pct": cr.get("identity_pct", "0.000"),
                "reference_coverage_pct": cr.get("reference_coverage_pct", "0.000"),
                "status": cr.get("status", "NA"),
                "confidence": cr.get("confidence", "NA"),
                "notes": cr.get("notes", "-"),
            })

    return row, consensus_record


def type_reads(db_path: Path, output_dir: Path, threads: int, min_qcov: float, min_pident: float, ambiguity_delta: float, min_depth: int, min_base_fraction: float, min_mapped_re
ads: int, candidate_limit: int = 25, min_candidate_breadth: float = 20.0, verbose: bool = False, illumina_r1: Optional[Path] = None, illumina_r2: Optional[Path] = None, illumina
_dir: Optional[Path] = None, ont_reads: Optional[Path] = None, ont_dir: Optional[Path] = None, skip_unpaired: bool = False) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    input_modes = sum(x is not None for x in [illumina_r1, illumina_dir, ont_reads, ont_dir])
    if input_modes != 1:
        die("Provide exactly one reads input mode: -1/-2, --illumina-dir, --ont, or --ont-dir.")
    samples = []
    pairing_rows: List[Dict[str, str]] = []
    if illumina_r1 is not None:
        if illumina_r2 is None:
            die("Paired-end Illumina mode requires both -1 and -2.")
        sample_id = sanitize_id(strip_fastq_extensions(illumina_r1))
        parsed = detect_illumina_pair_token(illumina_r1)
        if parsed:
            sample_id = parsed[0]
        samples.append((sample_id, "illumina", illumina_r1, illumina_r2))
        pairing_rows.append({"sample_id": sample_id, "R1": str(illumina_r1), "R2": str(illumina_r2), "status": "paired"})
    elif illumina_dir is not None:
        pairs, pairing_rows = discover_illumina_pairs(illumina_dir)
        bad = [r for r in pairing_rows if r["status"] != "paired"]
        if bad and not skip_unpaired:
            report_path = output_dir / "read_input_report.tsv"
            write_pairing_report(pairing_rows, report_path)
            die(f"Found unpaired or unrecognized FASTQ files. See {report_path}. Use --skip-unpaired to continue with complete pairs.")
        samples = [(sample, "illumina", r1, r2) for sample, r1, r2 in pairs]
    elif ont_reads is not None:
        sample_id = sanitize_id(strip_fastq_extensions(ont_reads))
        samples.append((sample_id, "ont", ont_reads, None))
        pairing_rows.append({"sample_id": sample_id, "R1": str(ont_reads), "R2": "NA", "status": "single_ont"})
    elif ont_dir is not None:
        ont_samples = discover_ont_reads(ont_dir)
        samples = [(sample, "ont", reads, None) for sample, reads in ont_samples]
        pairing_rows = [{"sample_id": sample, "R1": str(reads), "R2": "NA", "status": "single_ont"} for sample, reads in ont_samples]
    if not samples:
        die("No read samples found.")
    pairing_report = output_dir / "read_input_report.tsv"
    write_pairing_report(pairing_rows, pairing_report)
    log(f"Wrote read input report: {pairing_report}", verbose)
    log(f"Found {len(samples)} read sample(s) to process.", verbose)
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        ref_fasta, meta, ref_sequences = write_json_db_to_ref_fasta_and_meta(db_path, tmp_path, verbose=verbose)
        result_rows = []
        consensus_records = []
        for idx, (sample_id, read_type, r1, r2) in enumerate(samples, start=1):
            log(f"[{idx}/{len(samples)}] Processing read sample: {sample_id} ({read_type})", verbose)
            row, consensus_record = process_read_sample(
                sample_id=sample_id,
                read_type=read_type,
                read1=r1,
                read2=r2,
                db_path=db_path,
                ref_fasta=ref_fasta,
                meta=meta,
                ref_sequences=ref_sequences,
                output_dir=output_dir,
                threads=threads,
                min_qcov=min_qcov,
                min_pident=min_pident,
                ambiguity_delta=ambiguity_delta,
                min_depth=min_depth,
                min_base_fraction=min_base_fraction,
                min_mapped_reads=min_mapped_reads,
                candidate_limit=candidate_limit,
                min_candidate_breadth=min_candidate_breadth,
                verbose=verbose,
            )
            result_rows.append(row)
            if consensus_record is not None:
                consensus_records.append(consensus_record)
    out_tsv = output_dir / "gp60_read_typing_results.tsv"
    out_consensus = output_dir / "gp60_read_consensus_sequences.fasta"
    failed_tsv = output_dir / "failed_or_low_confidence_read_samples.tsv"
    fields = ["sample_id", "read_type", "R1_or_reads", "R2", "total_reads", "mapped_gp60_reads", "diagnostic_top_raw_read_reference", "diagnostic_top_raw_read_reference_reads", 
"mean_gp60_depth", "gp60_breadth_coverage_pct", "consensus_n_count", "likely_species", "gp60_family", "gp60_subtype", "best_reference", "accession", "identity_pct", "query_cover
age_pct", "reference_coverage_pct", "alignment_length_bp", "status", "confidence", "notes"]
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(result_rows)
    write_fasta(consensus_records, out_consensus)
    with failed_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in result_rows:
            if row["confidence"] != "high":
                writer.writerow(row)
    print(f"Wrote: {out_tsv}")
    print(f"Wrote: {out_consensus}")
    print(f"Wrote: {failed_tsv}")
    print(f"Wrote: {pairing_report}")
    print(f"Wrote BAM files in: {output_dir / 'bam'}")
    print(f"Wrote candidate diagnostics in: {output_dir / 'read_candidate_diagnostics'}")

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


    r = sub.add_parser("reads", help="Type gp60 directly from FASTQ reads using minimap2 and samtools.")
    read_group = r.add_mutually_exclusive_group(required=True)
    read_group.add_argument("-1", "--r1", type=Path, help="Illumina R1 FASTQ/FASTQ.GZ for one paired-end sample. Requires -2/--r2.")
    read_group.add_argument("--illumina-dir", type=Path, help="Directory with Illumina paired-end FASTQ files to auto-pair.")
    read_group.add_argument("--ont", type=Path, help="ONT FASTQ/FASTQ.GZ for one sample.")
    read_group.add_argument("--ont-dir", type=Path, help="Directory with ONT FASTQ/FASTQ.GZ files.")
    r.add_argument("-2", "--r2", type=Path, help="Illumina R2 FASTQ/FASTQ.GZ for one paired-end sample.")
    r.add_argument("--db", required=True, type=Path, help="Single gp60 JSON DB file.")
    r.add_argument("-o", "--output-dir", required=True, type=Path, help="Output directory.")
    r.add_argument("-t", "--threads", default=1, type=int, help="Threads for minimap2/samtools. Default: 1.")
    r.add_argument("--min-qcov", default=85.0, type=float, help="Minimum reference coverage percent for complete call. Default: 85.")
    r.add_argument("--min-pident", default=95.0, type=float, help="Minimum percent identity for confident call. Default: 95.")
    r.add_argument("--ambiguity-delta", default=1.0, type=float, help="Score delta for close competing references. Default: 1.0.")
    r.add_argument("--min-depth", default=3, type=int, help="Minimum read depth to call a consensus base. Default: 3.")
    r.add_argument("--min-base-fraction", default=0.70, type=float, help="Minimum fraction for majority consensus base. Default: 0.70.")
    r.add_argument("--min-mapped-reads", default=3, type=int, help="Minimum mapped gp60 reads needed to attempt consensus. Default: 3.")
    r.add_argument("--candidate-limit", default=25, type=int, help="Number of top read-mapped references to evaluate by consensus. Default: 25.")
    r.add_argument("--min-candidate-breadth", default=20.0, type=float, help="Minimum consensus breadth percent for a candidate to be preferred. Default: 20.")
    r.add_argument("--skip-unpaired", action="store_true", help="In --illumina-dir mode, continue with complete pairs and skip unpaired/unrecognized files.")
    r.add_argument("-v", "--verbose", action="store_true", help="Print progress messages to stderr.")

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
    elif args.mode == "reads":
        type_reads(
            db_path=args.db,
            output_dir=args.output_dir,
            threads=args.threads,
            min_qcov=args.min_qcov,
            min_pident=args.min_pident,
            ambiguity_delta=args.ambiguity_delta,
            min_depth=args.min_depth,
            min_base_fraction=args.min_base_fraction,
            min_mapped_reads=args.min_mapped_reads,
            candidate_limit=args.candidate_limit,
            min_candidate_breadth=args.min_candidate_breadth,
            verbose=args.verbose,
            illumina_r1=args.r1,
            illumina_r2=args.r2,
            illumina_dir=args.illumina_dir,
            ont_reads=args.ont,
            ont_dir=args.ont_dir,
            skip_unpaired=args.skip_unpaired,
        )


if __name__ == "__main__":
    main()
