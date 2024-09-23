import os
import csv
import subprocess

# Path to the assembly_stat.py script
assembly_stat_script = "/Path/to/assembly_stat.py"

# Directory containing sample folders with contigs.fasta files
sample_dir = "/Path/to/directory/fasta/files"

# Path to the output CSV file
output_csv = "/Path/to/stat_output_genome_assembly.csv"

# Initialize the CSV header
csv_header = [
    "Sample",
    "Total length of contigs >0bp",
    "Total length of contigs >500bp",
    "Number of contigs >0bp",
    "Number of contigs >500bp",
    "Largest contig",
    "Smallest contig",
    "Average contig size",
    "N50",
    "N90"
]

with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(csv_header)

# Iterate through sample folders
for folder in os.listdir(sample_dir):
    contigs_file = os.path.join(sample_dir, folder, f'spades/{folder}.contigs.fasta')

    # Run the assembly_stat.py script for the contigs.fasta file
    try:
        result = subprocess.check_output(["python3", assembly_stat_script, contigs_file], text=True)
        result_lines = result.strip().split("\n")
        result_dict = {line.split(": ")[0]: line.split(": ")[1] for line in result_lines}

        if not "A4138" in folder:
            result_dict["Sample"] = f'A4138_{folder}'
        else: 
            result_dict["Sample"] = folder

        with open(output_csv, "a", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_header)
            writer.writerow(result_dict)
    except subprocess.CalledProcessError:
        print(f"Error processing {contigs_file}")
