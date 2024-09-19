import random
import gzip
import os
from tqdm import tqdm
import argparse
from argparse import Namespace

# Simulation parameters
NUM_CELLS = 50           # Number of cells
NUM_UMIS = 100            # UMIs per cell
READ_LENGTH = 100         # Length of the DNA sequence
BARCODE_LENGTH = 10       # Length of the cell barcode
UMI_LENGTH = 8            # Length of the UMI
NUM_READS_PER_CELL = 100  # Number of reads per cell
CHROMOSOME_COUNT = 10     # Number of chromosomes in the simulated genome
GENES_PER_CHROM = 100     # Number of genes per chromosome
MIN_EXONS_PER_GENE = 1    # Minimum number of exons per gene
MAX_EXONS_PER_GENE = 5    # Maximum number of exons per gene
EXON_LENGTH = 200         # Fixed exon length for simplicity

# Generate a random DNA sequence (for the genome)
def generate_random_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

# Generate a reference genome with specified size (number of base pairs)
def generate_reference_genome(genome_size_bp):
    genome = {}
    for chrom in range(1, CHROMOSOME_COUNT + 1):
        chromosome_seq = generate_random_sequence(genome_size_bp // CHROMOSOME_COUNT)
        genome[f"chr{chrom}"] = chromosome_seq
    return genome

# Save reference genome to a FASTA file
def save_reference_genome(genome, output_file):
    with open(output_file, 'w') as f:
        for chrom, seq in genome.items():
            f.write(f">{chrom}\n")
            # Split sequence into lines of 80 characters (FASTA format standard)
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

# Generate a GTF file with simulated gene annotations
def generate_gtf(genome, gtf_output_file):
    with open(gtf_output_file, 'w') as gtf_file:
        for chrom, chrom_seq in genome.items():
            chrom_length = len(chrom_seq)
            for gene_id in range(1, GENES_PER_CHROM + 1):
                # Randomly place the gene on the chromosome
                gene_start = random.randint(0, chrom_length - EXON_LENGTH * MAX_EXONS_PER_GENE)
                gene_end = gene_start + random.randint(MIN_EXONS_PER_GENE, MAX_EXONS_PER_GENE) * EXON_LENGTH
                
                # Write gene entry
                gene_name = f"gene{gene_id}"
                gtf_file.write(f"{chrom}\tSim\tgene\t{gene_start+1}\t{gene_end}\t.\t+\t.\tgene_id \"{gene_name}\";\n")
                
                # Add exons for this gene
                exon_count = random.randint(MIN_EXONS_PER_GENE, MAX_EXONS_PER_GENE)
                exon_start = gene_start
                for exon_id in range(1, exon_count + 1):
                    exon_end = exon_start + EXON_LENGTH
                    gtf_file.write(f"{chrom}\tSim\texon\t{exon_start+1}\t{exon_end}\t.\t+\t.\tgene_id \"{gene_name}\"; exon_number \"{exon_id}\";\n")
                    exon_start = exon_end + 1

# Simulate reads based on a reference genome
def simulate_read_from_genome(genome, read_length):
    chrom = random.choice(list(genome.keys()))
    chrom_seq = genome[chrom]
    start_pos = random.randint(0, len(chrom_seq) - read_length)
    dna_seq = chrom_seq[start_pos:start_pos + read_length]
    return chrom, start_pos, dna_seq

# Generate unique barcodes and UMIs
def generate_barcode():
    return generate_random_sequence(BARCODE_LENGTH)

def generate_umi():
    return generate_random_sequence(UMI_LENGTH)

# Generate DNA read (R1) and cell barcode/UMI (R2) from reference genome
def generate_fastq_read(cell_barcode, umi, genome):
    
    r1_read:str
    r2_read:str
    quality:str
    
    chrom, start_pos, dna_seq = simulate_read_from_genome(genome, READ_LENGTH)
    quality = ''.join(['?'] * READ_LENGTH)
    
    # FASTQ format (R1 - DNA) with gene-tag XT
    r1_read = f"@{cell_barcode}:{umi}:R1:chrom={chrom}:start={start_pos} XT:Z:{umi}\n{dna_seq}\n+\n{quality}\n"
    
    # FASTQ format (R2 - Barcode + UMI)
    barcode_umi = cell_barcode + umi
    r2_read = f"@{cell_barcode}:{umi}:R2\n{barcode_umi}\n+\n{''.join(['?'] * len(barcode_umi))}\n"
    
    return r1_read, r2_read

# Generate FASTQ file for R1 and R2 based on the reference genome
def simulate_fastq_data(genome, num_cells, num_umis, num_reads_per_cell, r1_output, r2_output, whitelist_output) -> None:
    barcodes_set = set()  # Set to store unique barcodes for whitelist

    with gzip.open(r1_output, 'wt') as r1_file, gzip.open(r2_output, 'wt') as r2_file:
        for cell_id in tqdm(range(num_cells)):
            cell_barcode = generate_barcode()
            barcodes_set.add(cell_barcode)  # Store unique barcodes
            for umi_id in range(num_umis):
                umi = generate_umi()
                for read_id in range(num_reads_per_cell):
                    r1_read, r2_read = generate_fastq_read(cell_barcode, umi, genome)
                    r1_file.write(r1_read)
                    r2_file.write(r2_read)
    
    # Write unique barcodes to whitelist.txt
    with open(whitelist_output, 'w') as whitelist_file:
        for barcode in barcodes_set:
            whitelist_file.write(f"{barcode}\n")
            
def main(genome_dir:str, samples_dir:str, genome_size_bp:int=10_000_000, num_of_samples:int=3) -> None:
    
    # Create directories if they don't exist
    os.makedirs(genome_dir, exist_ok=True)
    os.makedirs(samples_dir, exist_ok=True)

    # Generate reference genome
    reference_genome = generate_reference_genome(genome_size_bp)
    
    # Save the reference genome to a FASTA file in the reference_genome directory
    genome_output_file = os.path.join(genome_dir, "simulated_genome.fa")
    save_reference_genome(reference_genome, genome_output_file)
    print(f"Reference genome saved to {genome_output_file}")

    # Generate GTF file for the reference genome
    gtf_output_file = os.path.join(genome_dir, "simulated_genes.gtf")
    generate_gtf(reference_genome, gtf_output_file)
    print(f"GTF file generated: {gtf_output_file}")

    # Simulate FASTQ data using the reference genome
    sample_names = [f'tissue{indx}' for indx in range(1, num_of_samples + 1)]  # List of sample names
    for sample_name in sample_names:
        # Create a directory for the sample
        sample_dir = os.path.join(samples_dir, sample_name)
        os.makedirs(sample_dir, exist_ok=True)

        # Output parameters
        r1_output = os.path.join(sample_dir, f"sample_test_1_{sample_name}_L001_R1_001.fastq.gz")
        r2_output = os.path.join(sample_dir, f"sample_test_1_{sample_name}_L001_R2_001.fastq.gz")
        whitelist_output = os.path.join(sample_dir, f"whitelist_sample_test_1_{sample_name}.txt")

        # Run simulation
        simulate_fastq_data(reference_genome, NUM_CELLS, NUM_UMIS, NUM_READS_PER_CELL, r1_output, r2_output, whitelist_output)

        print(f"FASTQ files generated in {sample_dir}: {r1_output}, {r2_output}")
        print(f"Whitelist file generated in {sample_dir}: {whitelist_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process genome and sample directories.")
    
    parser.add_argument('--genome_dir', type=str, required=True, help='Path to the genome directory')
    parser.add_argument('--samples_dir', type=str, required=True, help='Path to the samples directory')
    parser.add_argument('--genome_size_bp', type=int, required=True, help='Size of the genome in base pairs')
    parser.add_argument('--num_of_samples', type=int, required=True, help='Number_of_samples')
    
    args:Namespace = parser.parse_args()
    
    main(args.genome_dir, args.samples_dir, args.genome_size_bp, args.num_of_samples)
    

