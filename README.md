# CRYPTIC v1.0

**CRYPTIC** is a command-line tool for **Cryptosporidium gp60 typing and identification** from genome assemblies or already extracted gp60 sequences.

**CRYPTIC** stands for:

> **Cryptosporidium gp60 Typing and Identification Classifier**

The tool uses a single embedded JSON reference database containing gp60 reference sequences and metadata, including species, gp60 family, subtype, variant, accession, sequence length, and original reference header.

---

## Features

- Types **Cryptosporidium gp60** from genome assemblies.
- Supports already extracted gp60 FASTA sequences.
- Uses one self-contained JSON database file.
- Reports likely species, gp60 family, subtype, best reference, identity, coverage, confidence, and status.
- Extracts the best gp60 sequence from each genome.
- Flags partial, low-confidence, ambiguous, and not-found calls.
- Provides optional verbose progress logging.
- Does not require Biopython.

---

## Current version

```text
CRYPTIC v1.0
```

This first release supports:

```text
assembly/genome FASTA input
extracted gp60 FASTA input
single JSON database input
BLAST-based reference matching
```

FASTQ read-based typing is planned for a future version.

---

## Requirements

CRYPTIC requires Python 3 and BLAST+.

Check that BLAST+ is available:

```bash
blastn -version
makeblastdb -version
```

If BLAST+ is not installed, install it with conda:

```bash
conda install -c bioconda blast
```

No additional Python packages are required.

---

## Files needed

The simplest setup requires only two files:

```text
cryptic.py
gp60_db.json
```

Where:

```text
cryptic.py      CRYPTIC command-line script
gp60_db.json    single embedded gp60 reference database
```

---

## Quick start

Run CRYPTIC on genome assemblies:

```bash
python cryptic.py genome \
  -d assemblies \
  --db gp60_db.json \
  -o cryptic_results \
  -t 8 \
  --verbose
```

Run CRYPTIC on already extracted gp60 sequences:

```bash
python cryptic.py query \
  -i extracted_gp60s.fasta \
  --db gp60_db.json \
  -o cryptic_query_results \
  -t 8 \
  --verbose
```

---

## Input modes

CRYPTIC has two modes.

### 1. Genome mode

Genome mode searches the gp60 database against each genome assembly FASTA file.

```bash
python cryptic.py genome \
  -d assemblies \
  --db gp60_db.json \
  -o cryptic_results \
  -t 8
```

Input directory example:

```text
assemblies/
├── sample1.fasta
├── sample2.fasta
├── sample3.fa
└── sample4.fna
```

Accepted assembly extensions:

```text
.fa
.fna
.fasta
.fas
```

### 2. Query mode

Query mode classifies already extracted gp60 sequences, amplicons, or consensus sequences.

```bash
python cryptic.py query \
  -i extracted_gp60s.fasta \
  --db gp60_db.json \
  -o cryptic_query_results \
  -t 8
```

---

## Output files

### Genome mode output

```text
cryptic_results/
├── gp60_typing_results.tsv
├── extracted_gp60_sequences.fasta
├── failed_or_low_confidence_samples.tsv
└── blast_results/
```

### Query mode output

```text
cryptic_query_results/
├── gp60_typing_results.tsv
└── blast_results/
```

---

## Main output table

The main output file is:

```text
gp60_typing_results.tsv
```

Genome mode columns:

```text
sample_id
input_file
likely_species
gp60_family
gp60_subtype
best_reference
accession
identity_pct
reference_coverage_pct
alignment_length_bp
contig
contig_start
contig_end
strand
status
confidence
notes
```

Query mode columns:

```text
sample_id
likely_species
gp60_family
gp60_subtype
best_reference
accession
identity_pct
query_coverage_pct
reference_coverage_pct
alignment_length_bp
status
confidence
notes
```

---

## Interpreting results

### status

```text
complete       gp60 was detected with sufficient reference coverage and identity
partial        gp60 was detected but reference coverage was below threshold
low_identity   gp60-like sequence was detected but identity was below threshold
ambiguous      close competing references disagree at the family level
not_found      no gp60 hit was detected
```

### confidence

```text
high      strong gp60 call
medium    usable call but should be reviewed
low       weak, partial, or ambiguous call
none      no gp60 detected
```

---

## Default thresholds

By default, CRYPTIC uses:

```text
--min-qcov 85
--min-pident 95
```

This means a confident complete call requires at least:

```text
85% reference coverage
95% nucleotide identity
```

For stricter subtype-level typing, use:

```bash
python cryptic.py genome \
  -d assemblies \
  --db gp60_db.json \
  -o cryptic_results_strict \
  -t 8 \
  --min-qcov 95 \
  --min-pident 98 \
  --verbose
```

Recommended interpretation:

```text
Discovery / screening:
  --min-qcov 85 --min-pident 95

Stricter subtype calling:
  --min-qcov 95 --min-pident 98
```

---

## Verbose mode

Use `--verbose` to show progress.

```bash
python cryptic.py genome \
  -d assemblies \
  --db gp60_db.json \
  -o cryptic_results \
  -t 8 \
  --verbose
```

Example verbose output:

```text
[gp60typer] Found 10 genome FASTA file(s) in assemblies.
[gp60typer] Loading DB JSON: gp60_db.json
[gp60typer] Loaded 233 gp60 reference records from DB.
[gp60typer] [1/10] Processing sample: sample1
[gp60typer] [1/10] BLAST hits found: 24
[gp60typer] [1/10] Best hit: C_parvum_gp60_IIc_a_AF164491 species=C_parvum family=IIc subtype=IIc_a identity=99.800% qcov=98.500%
```

Save the progress log:

```bash
python cryptic.py genome \
  -d assemblies \
  --db gp60_db.json \
  -o cryptic_results \
  -t 8 \
  --verbose 2> cryptic_run.log
```

---

## Database format

CRYPTIC uses one JSON database file:

```text
gp60_db.json
```

Each record contains:

```text
ref_id
sequence
description
species
species_raw
gene
gp60_family
gp60_variant
gp60_subtype
accession
length_bp
original_header
```

Example record structure:

```json
{
  "ref_id": "C_parvum_gp60_IIc_a_AF164491",
  "sequence": "TGTCTCCGCTGT...",
  "species": "C_parvum",
  "species_raw": "C.parvum",
  "gene": "gp60",
  "gp60_family": "IIc",
  "gp60_variant": "a",
  "gp60_subtype": "IIc_a",
  "accession": "AF164491",
  "length_bp": 1020,
  "original_header": "C.parvum|IIc|a(AF164491)"
}
```

Internally, CRYPTIC temporarily writes the embedded reference sequences to FASTA because BLAST requires FASTA input. The user only needs to provide the single JSON database file.

---

## Example result

Example `gp60_typing_results.tsv` row (subtypes are still under work in progress):

```text
sample_id    likely_species    gp60_family    gp60_subtype    best_reference                  identity_pct    reference_coverage_pct    status      confidence
sample1      C_parvum          IIc            IIc_a           C_parvum_gp60_IIc_a_AF164491    99.800          98.500                    complete    high
```

---

## Important notes

`likely_species` is inferred from the best gp60 reference match. This should be interpreted as **gp60-based species inference**, not definitive species identification.

For final species assignment, confirm with an independent species marker or genome-wide method when possible.

Subtype-level calls are most reliable when the gp60 region is nearly complete and has high identity to a reference.

---

## Citation

No formal citation is available yet.

If you use CRYPTIC before publication, cite the GitHub repository and include the version used:

```text
CRYPTIC v1.0, Cryptosporidium gp60 Typing and Identification Classifier.
```

---

## Author

Developed by Rodrigo P. Baptista

---

## Version history

### v1.0

Initial release.

Features:

```text
single JSON gp60 database
genome assembly mode
extracted gp60 query mode
BLAST-based gp60 typing
confidence and ambiguity reporting
verbose progress logging
```
