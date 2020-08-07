
rule make_tsv_files_of_counts_and_pheno:
    input:
        "testdata/testdata_ab1.rds",
        #"results/cell_labelling_ab1/ab1_combined_identified.rds",
    output:
        "results/scenic/countmatrix.tsv",
        "results/scenic/metadata.tsv"
    conda:
        "../envs/R_env.yaml"
    script:
        "../scripts/scenic_generate_count_matrix.R"

rule grn_make:
    input:
        "results/scenic/countmatrix.tsv",
        rules.get_network_analysis_databases.output.networkdir,
    output:
        "results/scenic/adjecencies.p"
    resources:
        mem_mb=96000
    threads: 12
    conda:
        "../envs/arboreto_env.yaml"
    script:
        "../scripts/scenic_grn_construction.py"

rule regulon_processing_from_grn:
    input:
        "results/scenic/adjecencies.p",
        "results/scenic/countmatrix.tsv",
    output:
        "results/scenic/df_MODULES_PRUNED.p",
        "results/scenic/REGULONS_all.p"
    conda:
        "../envs/pyscenic.yaml"
    script:
        "../scripts/scenic_regulon_processing.py"

rule auc_cell_to_score_regulons_on_matrix:
    input:
        "results/scenic/REGULONS_all.p",
        "results/scenic/countmatrix.tsv"
    output:
        "results/scenic/auc_binarized.csv"
    conda:
        "../envs/pyscenic_latest.yaml"
    script:
        "../scripts/scenic_auc_processing.py"

rule clean_auc_output:
    input:
        "results/scenic/auc_binarized.csv",
    output:
        "results/scenic/binarised.rds",
    conda:
        "../envs/R_env.yaml"
    script:
        "../scripts/scenic_preprocessing_TFnames.R"
        
rule plot_regulon_heatmaps:
    input:
        "results/scenic/binarised.rds",
        "results/scenic/metadata.tsv"
    output:
        "results/scenic/TF_clustermap_small.pdf",
        "results/scenic/TF_clustermap_full.pdf"
    conda:
        "../envs/R_env.yaml"
    script:
        "../scripts/scenic_delta_script.R"
