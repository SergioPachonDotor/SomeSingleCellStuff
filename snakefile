SAMPLES = ["tissue1", "tissue2", "tissue3"]

rule all:
    input:
        expand("output/second_analysis/{sample}/deduplicated.bam", sample=SAMPLES)

rule simulate_data:
    output:
        "samples/tissue1/sample_test_1_tissue1_L001_R1_001.fastq.gz",
    shell:
        """
        python3 ./scripts/01_single_cell_data_simulator.py \
            --genome_dir ./reference_genome \
            --samples_dir ./samples \
            --genome_size_bp 5000000 \
            --num_of_samples {len(SAMPLES)}
        python3 ./scripts/02_output_dirs_creator.py \
            --samples_dir ./samples \
            --output_dir ./output
        """

rule genome_index:
    output:
        "reference_genome/Genome"
    params:
        genomeDir="./reference_genome/",
        genomeFastaFiles="./reference_genome/simulated_genome.fa",
        sjdbGTFfile="./reference_genome/simulated_genes.gtf"
    shell:
        """
        STAR --runThreadN 8 \
             --runMode genomeGenerate \
             --genomeDir {params.genomeDir} \
             --genomeFastaFiles {params.genomeFastaFiles} \
             --sjdbGTFfile {params.sjdbGTFfile} \
             --genomeSAindexNbases 10
        """

rule align_reads:
    input:
        genome_index="reference_genome/Genome",
        r1="samples/{sample}/sample_test_1_{sample}_L001_R1_001.fastq.gz",
        r2="samples/{sample}/sample_test_1_{sample}_L001_R2_001.fastq.gz",
        whitelist="samples/{sample}/whitelist_sample_test_1_{sample}.txt"
    output:
        "output/first_analysis/{sample}/STAR_output_Aligned.sortedByCoord.out.bam"
    params:
        out_prefix="output/first_analysis/{sample}/STAR_output_"
    shell:
        """
        STAR --genomeDir ./reference_genome/ \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.out_prefix} \
            --runThreadN 8 \
            --soloType Droplet \
            --soloCBwhitelist {input.whitelist} \
            --soloCBlen 10 \
            --soloUMIlen 8 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM CR UR \
            --soloFeatures Gene
        """

rule filter_bam:
    input:
        "output/first_analysis/{sample}/STAR_output_Aligned.sortedByCoord.out.bam"
    output:
        "output/first_analysis/{sample}/filtered.bam"
    shell:
        """
        samtools view -b -F 4 -q 20 {input} > {output}
        samtools index {output}
        """

rule deduplication:
    input:
        "output/first_analysis/{sample}/filtered.bam"
    output:
        "output/second_analysis/{sample}/deduplicated.bam"
    params:
        stats="output/second_analysis/{sample}/dedup_stats"
    shell:
        """
        umi_tools dedup -I {input} -S {output} --output-stats={params.stats}
        """
