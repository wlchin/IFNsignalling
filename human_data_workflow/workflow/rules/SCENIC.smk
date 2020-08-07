
rule grn_make:
    input:
        "results/monocytes.tsv",
        "resources/tfs.txt"
    output:
        "results/scenic/adjecencies.p"
    resources:
        mem_mb=96000
    threads: 12
    container:
        "docker://aertslab/pyscenic:0.9.19"
    script:
        "../scripts/scenic_grn_construction.py"

rule regulon_processing_from_grn:
    input:
        "results/scenic/adjecencies.p",
        "results/monocytes.tsv",
    output:
        "results/scenic/df_MODULES_PRUNED.p",
        "results/scenic/REGULONS_all.p"
    container:
        "docker://aertslab/pyscenic:0.9.19"
    script:
        "../scripts/scenic_regulon_processing.py"

rule auc_cell_to_score_regulons_on_matrix:
    input:
        "results/scenic/REGULONS_all.p",
        "results/monocytes.tsv"
    output:
        "results/scenic/auc_binarized.csv"
    container:
        "docker://aertslab/pyscenic:0.9.19"
    script:
        "../scripts/scenic_auc_processing.py"

rule clean_auc_output:
    input:
        "results/scenic/auc_binarized.csv",
        "results/metadata_for_scde.rds",
    output:
        "results/scenic/binarised_stuff_pbmc.rds"
        "results/binarized_mono.tsv"
    script:
        "../scripts/analysed_binarised.R"
