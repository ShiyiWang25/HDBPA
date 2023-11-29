# Toolkits for Nanopore sequencing data

## HDBPA (Hamming-distance-based phylogenetic analysis) 
(Last Update: 2023-11-22)
### Description
This _Python_ script is designed to trace the evolution of the HIV population on the single virus level.

It utilizes Hamming distance to calculate the similarity between sequences detected at different time points and could handle thousands of long reads in seconds.

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
python3 -m HDBPA --SAM1 ./TestData/DN.sam --SAM2 ./TestData/VF01.sam -r ./TestData/hxb2.fa
```
**Expected outputs:**
<img align="right" src="https://github.com/ShiyiWang25/HDBPA/blob/main/Figures/HDBPA_plot.png" width=20% height=20%>

Two significant subpopulations formed during this treatment. 

They are derived from the MRCAs with mutational patterns A518 and A1594, respectively. 

More details can be found in the CSV file in the Figures folder.




## Further application: analysis across multiple time points
The script for analysis across multiple time points will be updated soon. 

Here is an example of viral evolution across various time points and therapies by manually combining evolution between every two timepoints and plotting using a Fish plot.

<img align="center" src="https://github.com/ShiyiWang25/HDBPA/blob/main/Figures/Example.png" width=80% height=80%>


