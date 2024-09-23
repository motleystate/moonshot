import os

# load file containing names of samples 
list_samples_file = "list_samples.txt"

with open(list_samples_file, 'r') as file:
    SAMPLES = [line.strip() for line in file]

#create directory for each sample's results
for sam in SAMPLES:
    sample_dir = os.path.join("results", sam)
    subdirectories = ["bowtie2", "spades","prokka"]  
    os.makedirs(sample_dir, exist_ok=True)
    for subdir in subdirectories:
        subdir_path = os.path.join(sample_dir, subdir)
        os.makedirs(subdir_path, exist_ok=True)

dirlist = ["panaroo","iqtree"]
for name in dirlist:
    os.makedirs(name, exist_ok=True)

# load a file containing mappings from the names in “list_samples.txt” to the names in the fasta file  
sample_corres = "sample_correspondance.txt"

sample_dict = {}
with open(sample_corres, "r") as sample_file:
    for line in sample_file:
        cols = line.strip().split("\t")
        sample_dict[cols[0]] = cols[1]

def getrawreads(sample):
    prefix = sample_dict[sample]
    raw_read_file = [f'/Path/to/fasta/files/{prefix}_1.fastq.gz', f'/Path/to/fasta/files/{prefix}_2.fastq.gz']
    return raw_read_file

rule all:
    input:
        "pipeline_end.txt"

rule aligned_via_bowtie2: 
    input:
        rawreads = lambda wildcards: getrawreads(wildcards.sample)
    params:
        btindex = "/Path/to/CustomDB/for/filtering/Bowtie2_DB/species",
        log = "results/{sample}/bowtie2/{sample}_bowtie2.log"
    output:
        samfile = "results/{sample}/bowtie2/{sample}.sam"
    shell:
        """
        bowtie2 -p 4 \
        -x {params.btindex} -1 {input.rawreads[0]} -2 {input.rawreads[1]} \
        -S {output.samfile} 2> {params.log} 
        """

rule filtered_bam:
    input: "results/{sample}/bowtie2/{sample}.sam"
    output: "results/{sample}/bowtie2/{sample}.filtered.bam"
    shell:
    # exclude reads where neither pair is mapped with -F 12 (UNMAP, MUNMAP)
    # -f include, -F exclude
    # use 'samtools flags 12' to check. Can use any number
        """
        samtools view -SbF 12 {input} | \
        samtools sort -o {output}
        """

rule extract_mapped_reads:
    input: "results/{sample}/bowtie2/{sample}.filtered.bam"
    output:
        read1 = "results/{sample}/bowtie2/{sample}.aligned.1.fastq.gz",
        read2 = "results/{sample}/bowtie2/{sample}.aligned.2.fastq.gz"
    shell:
    ### samtools collate – shuffles and groups reads together by their names
    ### -u : Write uncompressed BAM output
    ### -O : Output to stdout.
    ### samtools fastq - converts a SAM/BAM/CRAM file to FASTA or FASTQ
    ### -n By default, either '/1' or '/2' is added to the end of read names where the corresponding READ1 or READ2 FLAG bit is set. 
    ### Using -n causes read names to be left as they are.
    ### -s write singleton reads to FILE
    ### -0 Write reads where the READ1 and READ2 FLAG bits set are either both set or both unset to FILE instead of outputting them.
        """
        samtools collate -Ou {input} | \
        samtools fastq -1 {output.read1} -2 {output.read2} -0 /dev/null -s /dev/null -n
        """

rule assembly_via_spades:
    input: 
        read1 = rules.extract_mapped_reads.output.read1,
        read2 = rules.extract_mapped_reads.output.read2
    output: 
        contigs = "results/{sample}/spades/{sample}.contigs.fasta",
        sample_dir = directory("results/{sample}/spades/")
    threads:
        8
    shell:
        """
        spades.py --careful --cov-cutoff auto --threads {threads} \
        -1 {input.read1} -2 {input.read2} -o {output.sample_dir} &&
        mv "results/{wildcards.sample}/spades/contigs.fasta" "results/{wildcards.sample}/spades/{wildcards.sample}.contigs.fasta"
        """

rule annotation_via_prokka:
    input:
        contigs = "results/{sample}/spades/{sample}.contigs.fasta"
    output:
        result_file = "results/{sample}/prokka/{sample}.gff",
        result_dir = directory("results/{sample}/prokka/")
    shell:
        """
        prokka --force --kingdom Bacteria --prefix {wildcards.sample} --locustag {wildcards.sample}_ --outdir {output.result_dir} {input.contigs}
        """

rule pangenome_panaroo:
    input: expand("results/{sample}/prokka/{sample}.gff", sample = SAMPLES)
    output:
        outdir = directory("panaroo"),      
        core_aln = "panaroo/core_gene_alignment.aln",
        pre_abs_table = expand("panaroo/gene_presence_absence.{ext}", ext = ["csv","Rtab"])
    threads:
        12
    shell:
        """
        panaroo -i {input} --clean-mode moderate -o {output.outdir} -t {threads} -a core --aligner mafft --core_threshold 0.50

        """

rule iqtree:
    input: rules.pangenome_panaroo.output.core_aln
    output: "panaroo/core_gene_alignment.aln.treefile"
    threads:
        12
    shell:
        """
        iqtree -s {input} -m MFP -bb 1000 -nt {threads} -ntmax {threads} -st DNA
        """

rule run_completed:
    input:
        "panaroo/core_gene_alignment.aln.treefile"
    output:
        end_report = "pipeline_end.txt"
    shell:
        """
        touch {output.end_report}
        """
