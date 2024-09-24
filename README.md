# moonshot : Metagenome assembly pipeline 

## I-Snakemake Pipeline for Genome Assembly, Annotation, and Pangenome Analysis
This pipeline automates the process of read filtering, genome assembly, annotation, and pangenome analysis. It is designed to process multiple samples and produce reproducible results in a structured directory.

**Note**: Due to the computational intensity of genome assembly, it is recommended to run this pipeline on a high-performance computing cluster (HPC). Running the pipeline on a local machine might require substantial resources (RAM, CPU threads), particularly when handling large datasets.

## Requirements
- Python 3.x
- Snakemake
- Bowtie2
- Samtools
- SPAdes
- Prokka
- Panaroo
- IQ-TREE

Before running the pipeline on a cluster, make sure to load the necessary modules. You can use the following commands to load the required software on your environment:
```bash
module load Python/3.8.1 snakemake bowtie2 samtools SPAdes panaroo IQ-TREE
```
**Note**:Depending on the modules already installed in the environment, the following modules must be installed or loaded using the appropriate `module load` commands:
```bash
module load graalvm minced hmmer bedtools aragorn blast+ infernal prodigal ncbitools barrnap signalp cd-hit
```

## Directory Structure
For each sample listed in `list_samples.txt`, the pipeline will create a structured directory with subdirectories for the outputs of each tool.
```
results/
  ├── sample1/
  │   ├── bowtie2/
  │   ├── spades/
  │   └── prokka/
  ├── sample2/
  │   ├── bowtie2/
  │   ├── spades/
  │   └── prokka/
  ...
panaroo/
iqtree/
```
## Input Files
- `list_samples.txt`: A list of sample names, one per line.
- `sample_correspondance.txt`: A tab-delimited file mapping sample names from `list_samples.txt` to their corresponding FASTA filenames.

Example `sample_correspondance.txt`:
```
IMVM_0001	A4138_0001
IMVM_0002	A4138_0002
```

- `CustomDB` : This database will be used to filter the reads, ensuring that only reads mapping to this database are included in the metagenome-assembled genomes (MAGs). This database typically consists of reference genomes, and can be expanded to include environmental or clinical strains. The FASTA files should be downloaded, and the database must be constructed using Bowtie2 with the following command:
```bash
bowtie2-build <reference_fasta> <CustomDB_name>
```
Exemple: 
`Ecoli_customDB_Accession#`: E. coli genome accesssion numbers used to build a customDB used to assemble E. coli MAGs from metagenomic data of vaginal samples. This database includes both reference and clinical strains.

##  Launch the Snakemake Pipeline

```bash
sbatch -p <partition> -c <number_of_cores> snakemake -p -j <number_of_jobs> -s <snakemake_file> 
```
`<partition>`: The partition or resource group you want to use.
`<number_of_cores>`: The number of cores you wish to allocate.
`<number_of_jobs>`: The maximum number of jobs to run in parallel.
`<snakemake_file>`: The name of the Snakemake file (usually Snakefile).

## II-Quality Assessment of the generated MAGs
At the end of the assembly process, general quality metrics of the generated MAGs can be calculated using two scripts: `run_stat.py` and `assembly_stat.py`.

### assembly_stat.py
This script computes statistics for a given contigs FASTA file. It parses the contigs, calculates lengths, and determines the N50 and N90, which are helpful for assessing the quality of the assembled genomes.

### run_stat.py
This script iterates through sample folders containing contigs.fasta files generated during the assembly process. It invokes assembly_stat.py for each file, collects various assembly statistics, and writes the results to a CSV file. Key statistics include:
- Total length of contigs greater than 0 bp and 500 bp
- Number of contigs greater than 0 bp and 500 bp
- Largest and smallest contig sizes
- Average contig size
- N50 and N90 values

## Acknowledgements and contributions 
We thank **Thi Ngoc Anh Vu**, who initially developed this pipeline under the supervision of **Sean P. Kennedy**. The work on assembling E. coli MAGs from vaginal samples of pregnant women was carried out by **Nassim Boutouchent**.

