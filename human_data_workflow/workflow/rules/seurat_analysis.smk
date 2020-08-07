
rule make_seurat_pbmc:
    input:
        "resources/GSE130157_14546.14554.14562.RawCounts.txt",
        "resources/GSE130157_15424R.15435R.RawCounts.txt",
        "resources/biomart_map_12072020.rds"
    container:
        "docker://wlc27/singler_legacy:v2"
    benchmark:
        "benchmarks/data_download.txt"
    output:
        "results/pbmc.rds",
        "results/fulldf.rds"
    script:
        "../scripts/create_pbmc_object_v2.R"

rule stuff_of_label:
    input:
        "results/pbmc.rds",
        "results/fulldf.rds",
        "resources/hpca_real.rda"
    container:
        "docker://wlc27/singler_legacy:v2"
    output:
        "results/labels_clusters.rds",
        "results/picture_labels.jpeg",
        "results/singler_result.rds",
        "results/pbmc_human_label_seurat_object.rds",
    script:
        "../scripts/single_r_labelling.R"

rule projection_experiment:
    input:
        "resources/ab1_combined_identified.rds",
        "results/pbmc.rds",
    container:
        "docker://wlc27/singler_legacy:v2"
    benchmark:
        "benchmarks/projection.txt"
    output:
        "results/predictions.rds",
        "results/prediction_projection.jpeg",
        "results/metadata_pbmc.rds",
    script:
        "../scripts/projection_experiment.R"

rule merge_metadata_and_prepare_counts:
    input:
        "results/pbmc.rds",
        "results/fulldf.rds",
        "results/predictions.rds",
        "results/labels_clusters.rds"
        "resources/metadata_human.csv",
    container:
        "docker://wlc27/singler_legacy:v2" 
    output:
        "results/metadata_for_scde.rds",
        "results/metadata_for_scde.tsv",
        "results/counts_monocytes.rds",
        "results/monocytes.tsv"
    script:
        "../scripts/merge_results_with_seurat_and_generate_counts.R"

rule benchmark_summary:
    input:
        "benchmarks/predictions.txt",
        "benchmarks/data_download.txt",
        "benchmarks/projection.txt",
    output:
        "benchmarks/seurat_workflow_summary.txt"
    shell:
        "cat {input} > {output}"
        
