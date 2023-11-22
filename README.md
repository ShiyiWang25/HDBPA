# Toolkits for Nanopore sequencing data

## HDBPA (Hamming-distance-based phylogenetic analysis) (Last UpdateL 2023-11-22)
### Description
This _Python_ script is designed to trace the evolution of the HIV population on the single virus level.
It utilizes Hamming distance to calculate the similarity between sequences detected at different time points and could handle thousands of long reads in seconds.h
Hamming distance is calculated on the non-synonymous mutation level in the current version. The calculation on the nucleotide level will be added soon.

### Features
1. HDBPA requires sequencing data collected at two different time points (e.g., pre-treatment and post-treatment)
2. HDBPA takes SAM files as input.
3. One reference sequence in the FASTA file format is required.
4. The phylogenetic analysis outcome will be stored in one CSV file.
5. A PNG figure will show the clonal expansion of the top 10 major subpopulations in the viral swarm.

### Installation
```
python -m venv env
source env/bin/activate
python -m pip install git+https://github.com/ShiyiWang25/HDBPA.git
```
### Command Line ARGS:

| Options | Description | Default |
| --- | --- | --- |
| `--SAM1`  | Input aligned reads in the 1st sample. | |
| `--SAM2`  | Input aligned reads in the 2nd sample. | |
| `-r`  | Import the reference sequence. | |
| `--mode`| Calculate the Hamming distance on the nucleotide level  or non-synonymous amino acid level | {NT, AA} |

### Example:
```
python3 -m HDBPA --SAM1 ./materials/DN.sam --SAM2 ./materials/VF01.sam -r ./materials/hxb2.fa
```
Expected outputs:


