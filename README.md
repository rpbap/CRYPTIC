# GP60typer

## Description

This tool uses Blastn to retrieve *Cryptosporidium* GP60 sequences from the genomes and compare to the reference database provided. 
Best hits would be kept and the equivalent GP60 type will be stored in typing_log.csv
The tool is not perfect, but with tests against reported gp60 types from Genomes from GenBank we got 98% match to the reported *C. parvum* gp60 types.

**Note for users**
* Refdb is currently only having *C. parvum* GP60 types, others will be added in new releases
* The --tree flag is disabled since it was malfunctioning don't use it until new release

## Usage

```
python typer.py -d <directory containing all genomes> -r Refdb.fasta -o <output_directory>
```
