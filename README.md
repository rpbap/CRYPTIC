# CRYPTIC v1.1

**CRYPTIC**: **Cryptosporidium gp60 Typing and Identification Classifier**
<p align="center">
<img width="300" height="300" alt="CRYPTIC_color" src="https://github.com/user-attachments/assets/f8b4a5c6-87d4-4d36-84ad-e93ddad82001" />
</p>

CRYPTIC is a command-line tool for reference-based **Cryptosporidium gp60 typing** from genome assemblies, extracted gp60 sequences, and sequencing reads.

Version **v1.1** supports:

- genome/assembly FASTA typing
- extracted gp60 FASTA typing
- Illumina paired-end FASTQ typing
- ONT FASTQ typing
- embedded single-file JSON gp60 database
- batch processing
- Illumina FASTQ pair auto-detection
- confidence/status reporting
- extracted/read-derived gp60 consensus output

---

## Overview

CRYPTIC classifies gp60 sequences using a curated embedded gp60 reference database. It reports the likely gp60-associated species, gp60 family, subtype, best reference, identity, coverage, call status, and confidence.

The tool currently uses **gp60-based typing**. Therefore, the `likely_species` column should be interpreted as the species assignment supported by the gp60 reference match, not as a full genome-based species identification.

---

## Main modes

CRYPTIC has three main modes:

```text
genome   Type gp60 directly from genome assemblies or contigs
query    Type already extracted gp60 sequences, amplicons, or consensus FASTA
reads    Type gp60 directly from Illumina or ONT FASTQ reads
```

---

## Requirements

### Required for genome and query modes

```text
Python >= 3.8
BLAST+
```

Required BLAST+ executables:

```text
blastn
makeblastdb
```

### Additional requirements for reads mode

```text
minimap2
samtools
```

Reads mode uses minimap2 for mapping and samtools for BAM processing, depth estimation, and reference-guided consensus generation.

### Conda installation

```bash
conda create -n cryptic -c conda-forge -c bioconda python=3.10 blast minimap2 samtools
conda activate cryptic
```

---

## Repository structure

Suggested repository layout:

```text
CRYPTIC/
├── cryptic.py
├── gp60_db.json
├── README.md
└── example_data/
```

For your current GitHub release, you can keep the latest tested script as:

```text
cryptic_v1.1.py
```

or rename it to the simpler command name:

```bash
mv cryptic_v1.1.py cryptic.py
```

All examples below use `cryptic.py`. Replace it with `cryptic_v1.1.py` if you keep the versioned file name.

---

## Quick start

### 1. Genome/assembly mode

Use this mode when you have assembled genomes, contigs, or scaffolds in FASTA format.

```bash
python cryptic.py genome \
  -d assemblies/ \
  --db gp60_db.json \
  -o cryptic_genome_results \
  -t 8 \
  -v
```

Input directory example:

```text
assemblies/
├── sample1.fasta
├── sample2.fa
└── sample3.fna
```

Supported assembly extensions include:

```text
.fasta
.fa
.fna
.fas
```

Main outputs:

```text
cryptic_genome_results/
├── gp60_typing_results.tsv
├── extracted_gp60_sequences.fasta
├── failed_or_low_confidence_samples.tsv
└── blast_results/
```

---

### 2. Extracted gp60/query mode

Use this mode when you already have gp60 sequences, amplicons, or consensus sequences in FASTA format.

```bash
python cryptic.py query \
  -i extracted_gp60s.fasta \
  --db gp60_db.json \
  -o cryptic_query_results \
  -t 8 \
  -v
```

Main outputs:

```text
cryptic_query_results/
├── gp60_query_typing_results.tsv
└── failed_or_low_confidence_queries.tsv
```

---

### 3. Illumina paired-end reads mode

Use this mode when you have paired-end Illumina FASTQ files.

#### Single sample

```bash
python cryptic.py reads \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --db gp60_db.json \
  -o cryptic_reads_results \
  -t 8 \
  -v
```

#### Batch folder with automatic pair detection

```bash
python cryptic.py reads \
  --illumina-dir reads/Illumina/ \
  --db gp60_db.json \
  -o cryptic_reads_results \
  -t 8 \
  -v
```

Supported Illumina pair naming examples:

```text
sample_R1.fastq.gz        sample_R2.fastq.gz
sample_R1_001.fastq.gz    sample_R2_001.fastq.gz
sample_1.fastq.gz         sample_2.fastq.gz
sample.1.fq.gz            sample.2.fq.gz
sample-1.fastq            sample-2.fastq
sample_1.fq               sample_2.fq
```

Supported FASTQ extensions:

```text
.fastq
.fq
.fastq.gz
.fq.gz
```

If some files are unpaired or have unrecognized names, CRYPTIC will stop and write a read input report. To continue with complete pairs only, use:

```bash
python cryptic.py reads \
  --illumina-dir reads/Illumina/ \
  --db gp60_db.json \
  -o cryptic_reads_results \
  -t 8 \
  --skip-unpaired \
  -v
```

---

### 4. ONT reads mode

#### Single ONT FASTQ

```bash
python cryptic.py reads \
  --ont sample.fastq.gz \
  --db gp60_db.json \
  -o cryptic_ont_results \
  -t 8 \
  -v
```

#### ONT batch folder

```bash
python cryptic.py reads \
  --ont-dir reads/ONT/ \
  --db gp60_db.json \
  -o cryptic_ont_results \
  -t 8 \
  -v
```

---

## Reads-mode workflow

Reads mode performs a conservative reference-guided gp60 typing workflow:

```text
FASTQ reads
  -> minimap2 mapping to embedded gp60 DB
  -> sorted/indexed BAM with samtools
  -> candidate gp60 references selected from mapped reads
  -> reference-guided consensus for top candidates
  -> consensus classification against gp60 DB
  -> final gp60 call with status/confidence
```

The final read-mode call is **not** selected by raw mapped-read count alone. CRYPTIC evaluates candidate references using breadth coverage, reference coverage, identity, and depth. This avoids false calls from reads mapping only to short conserved gp60 regions.

---

## Reads-mode outputs

```text
cryptic_reads_results/
├── gp60_read_typing_results.tsv
├── gp60_read_consensus_sequences.fasta
├── failed_or_low_confidence_read_samples.tsv
├── read_input_report.tsv
├── bam/
└── read_candidate_diagnostics/
```

### `gp60_read_typing_results.tsv`

This is the main reads-mode result table.

Important interpretation columns:

```text
sample_id
likely_species
gp60_family
gp60_subtype
best_reference
accession
identity_pct
reference_coverage_pct
gp60_breadth_coverage_pct
status
confidence
notes
```

Diagnostic columns are included for troubleshooting:

```text
diagnostic_top_raw_read_reference
diagnostic_top_raw_read_reference_reads
```

These diagnostic columns show the reference with the highest raw mapped-read count before candidate evaluation. They are **not** the final biological typing call.

### `read_candidate_diagnostics/`

This folder contains one candidate-reference diagnostics table per sample. It is useful for troubleshooting ambiguous or low-confidence calls.

Example file:

```text
read_candidate_diagnostics/sample1.candidate_refs.tsv
```

---

## Output interpretation

### Key columns

| Column | Meaning |
|---|---|
| `sample_id` | Sample name derived from input file name |
| `likely_species` | gp60-supported likely species |
| `gp60_family` | gp60 family, for example `IIc` |
| `gp60_subtype` | gp60 subtype/variant when resolvable |
| `best_reference` | Best matching gp60 reference record |
| `accession` | Accession associated with the best reference |
| `identity_pct` | Percent identity to the best reference |
| `reference_coverage_pct` | Percent of the reference covered by the alignment |
| `gp60_breadth_coverage_pct` | For reads mode, percent of selected gp60 reference covered by consensus bases |
| `status` | Call status |
| `confidence` | Confidence category |
| `notes` | Additional warnings or ambiguity notes |

### Status values

| Status | Meaning |
|---|---|
| `complete` | High-coverage gp60 match meeting default thresholds |
| `partial` | gp60 detected but coverage/breadth below complete-call threshold |
| `ambiguous` | Close competing references prevent a unique family or subtype call |
| `low_identity` | gp60-like hit detected but identity is below threshold |
| `not_found` | No acceptable gp60 hit detected |

### Confidence values

| Confidence | Meaning |
|---|---|
| `high` | Strong identity and coverage support |
| `medium` | gp60 detected but with subtype ambiguity or moderate limitations |
| `low` | Weak, partial, low-identity, or ambiguous evidence |
| `none` | No gp60 call |

---

## Thresholds

Default thresholds:

```text
--min-qcov 85
--min-pident 95
```

For stricter subtype-level calls, especially for publication-quality reporting:

```bash
python cryptic.py genome \
  -d assemblies/ \
  --db gp60_db.json \
  -o cryptic_genome_strict \
  -t 8 \
  --min-qcov 95 \
  --min-pident 98 \
  -v
```

For reads mode, additional parameters are available:

```text
--min-depth              Minimum read depth to call a consensus base. Default: 3
--min-base-fraction      Minimum majority base fraction for consensus. Default: 0.70
--min-mapped-reads       Minimum mapped gp60 reads to attempt consensus. Default: 3
--candidate-limit        Number of top mapped references to evaluate. Default: 25
--min-candidate-breadth  Minimum candidate breadth to prefer a candidate. Default: 20
```

Example stricter reads-mode run:

```bash
python cryptic.py reads \
  --illumina-dir reads/Illumina/ \
  --db gp60_db.json \
  -o cryptic_reads_strict \
  -t 8 \
  --min-qcov 95 \
  --min-pident 98 \
  --min-depth 5 \
  --candidate-limit 40 \
  -v
```

---

## Example result

Example reads-mode output:

```text
sample_id    likely_species    gp60_family    gp60_subtype          best_reference                    identity_pct    reference_coverage_pct    status      confidence
SRR15694494  C_parvum          IIc            IIc_d                 C_parvum_gp60_IIc_d_AF440636      100.000         99.884                    complete    high
SRR15694501  C_parvum          IIc            ambiguous:IIc_a|IIc_j C_parvum_gp60_IIc_a_AF164491      100.000         99.887                    complete    medium
SRR15694513  C_parvum          IIc            ambiguous:IIc_a|IIc_j C_parvum_gp60_IIc_a_AF164491      100.000         99.662                    complete    medium
```

Subtype ambiguity means that more than one closely related subtype reference is nearly indistinguishable under the current scoring thresholds. In this case, the family-level call may still be reliable while exact subtype assignment should be interpreted conservatively.

---

## Database

CRYPTIC uses a single embedded JSON database:

```text
gp60_db.json
```

Each record contains the reference sequence and metadata:

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

The script materializes this JSON database into a temporary FASTA internally when BLAST or minimap2 requires FASTA input. Users only need to provide the JSON database.

---

## Important limitations

1. **gp60 is a typing marker, not a whole-genome species classifier.**  
   The `likely_species` result is based on gp60 similarity.

2. **Reads mode is reference-guided.**  
   It is useful when assemblies are unavailable, but exact subtype calls should be interpreted carefully when coverage is low or when subtype references are highly similar.

3. **ONT reads may require more conservative interpretation. (still under tests and not validated)**  
   Homopolymers, repeats, and local errors can affect exact subtype-level calls.

4. **Subtype ambiguity can be biologically meaningful. (Still under development)**  
   Closely related gp60 subtype references may not be distinguishable from short or highly similar sequences.
   Subtype feature is still in process for nomenclatures to better be represented (like adding the specific mutation), for this current version we suggest     avoiding using them.

6. **Reference database completeness matters.**  
   CRYPTIC can only classify against the references present in `gp60_db.json`.
   We are keeping the database as uptaded as possible, but some genomes might be underrepresented and give you low confidence calls

---

## Recommended reporting language

For high-confidence calls:

```text
CRYPTIC assigned sample X to Cryptosporidium parvum gp60 family IIc, subtype IIc_d, with 100% identity and 99.9% reference coverage to reference AF440636.
```

For ambiguous subtype calls:

```text
CRYPTIC assigned sample X to Cryptosporidium parvum gp60 family IIc, with an ambiguous subtype-level call between IIc_a and IIc_j due to nearly identical reference matches.
```

For partial/low-confidence calls:

```text
CRYPTIC detected a partial gp60 signal in sample X, but coverage and/or identity were below thresholds for a confident subtype assignment.
```

---

## Troubleshooting

### `blastn` or `makeblastdb` not found

Install BLAST+:

```bash
conda install -c bioconda blast
```

### `minimap2` or `samtools` not found

Install reads-mode dependencies:

```bash
conda install -c bioconda minimap2 samtools
```

### Illumina files are not paired

Check:

```text
read_input_report.tsv
```

Use `--skip-unpaired` if you want CRYPTIC to continue with complete pairs only.

### Reads mode reports low coverage

Check:

```text
read_candidate_diagnostics/
gp60_read_consensus_sequences.fasta
bam/
```

Low breadth may indicate that gp60 is poorly covered, absent, fragmented, or that only conserved regions are represented in the reads.

---

## Citation

If you use CRYPTIC, please cite the associated repository and manuscript/preprint when available.

Suggested software citation placeholder:

```text
Baptista RP. CRYPTIC: Cryptosporidium gp60 Typing and Identification Classifier. Version 1.1.
```

---

## Author

Developed by **Rodrigo P. Baptista**

---

## Version history

### v1.1

- Added FASTQ reads mode
- Added Illumina paired-end single-sample mode
- Added Illumina batch auto-pairing
- Added ONT single-file and batch-folder support
- Added read-derived gp60 consensus output
- Added candidate-reference diagnostics for reads mode
- Improved reads-mode candidate selection using breadth, coverage, identity, and depth rather than raw mapped-read count alone
- Cleaned final reads-mode output by keeping raw-read reference information as diagnostic columns only

### v1.0

- Initial genome assembly mode
- Initial extracted gp60/query FASTA mode
- Single embedded JSON database support
- BLAST-based gp60 reference classification
- Status and confidence reporting
