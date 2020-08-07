
FASTA = "resources/gencode.vM21.transcripts.fa.gz"
INDEX = "resources/gencode.idx"
GENCODE_MAP = "resources/gene_to_transcript_map_gencode.rds"

IDS, = glob_wildcards("resources/raw/{id}_R1.fastq.gz")

rule create_kallisto_index:
    input: FASTA
    output: INDEX
    conda:
        "env/kallisto_and_sleuth.yml"
    shell:
        "kallisto index -i {output} {input}"

rule kallisto_alignment_and_quantification:
    input:
        "resources/raw/{sample}_R1.fastq.gz",
        INDEX
    output:
        directory("results/kallisto_aligned/{sample}")
    conda:
        "../env/kallisto_and_sleuth.yml"
    cache:
        True
    shell:
        "kallisto quant --single "
        "-l 200 "
        "-s 20 "
        "-i {INDEX} "
        "-b 100 "
        "-t 8 "
        "-o {output[0]} "
        "{input[0]}"

rule DE_analysis_with_transcript_aggregation:
    input:
        "resources/{strain}_metadata.csv",
        GENCODE_MAP,
        expand("results/kallisto_aligned/{sample}/abundance.tsv", sample = IDS)
    output:
        "results/DE/{strain}_DE_list.txt",
        "results/DE/{strain}_DE_gene_table.rds",
	"results/DE/{strain}_sleuth_object.rds",
	"results/DE/{strain}_phenodata.rds"
    conda:
        "../env/kallisto_and_sleuth.yml"
    script:
        "../scripts/DEscript.R"

rule create_gene_count_matrices_for_GENIE3:
    input:
        GENCODE_MAP,
        "results/DE/{strain}_phenodata.rds"
    output:
        "results/mats/{strain}_respmat.rds",
        "results/mats/{strain}_nonrespmat.rds"
    conda:
        "../env/tximport_and_dplyr.yml"
    script:
        "../scripts/create_count_matrices.R"
