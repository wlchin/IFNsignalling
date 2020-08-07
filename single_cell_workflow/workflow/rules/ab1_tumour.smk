
rule single_r_label_ab1:
    input:
        "results/combined_seurat_objects/ab1_combined.rds"
    conda:
        "../envs/bioconductor.yaml"
    output:
        "results/cell_labelling_ab1/ab1_single_r_labels.rds",
        "results/cell_labelling_ab1/ab1_labels.jpeg",
    script:
        "../scripts/singler_cell_labelling.R"
        
rule tumour_predictions_ab1:
    input:
        "results/combined_seurat_objects/ab1_combined.rds",
        "results/combined_seurat_objects/ab1_tagged_combined.rds"
    container:
        "docker://wlc27/singler_legacy:v2"
    output:
        "results/cell_labelling_ab1/ab1_metadata_predictions.rds",
        "results/cell_labelling_ab1/ab1_prediction_scores.jpeg"
    script:
        "../scripts/identify_HAtagged_cells.R"

rule analyse_metadata_ab1:
    input:
        "results/combined_seurat_objects/ab1_combined.rds",
        "results/cell_labelling_ab1/ab1_single_r_labels.rds",
        "results/cell_labelling_ab1/ab1_metadata_predictions.rds"
    conda:
        "../envs/bioconductor.yaml"
    output:
        "results/cell_labelling_ab1/ab1_combined_identified.rds",
        "results/cell_labelling_ab1/tumour_ident_and_singleR_ab1.jpeg"
    script:
        "../scripts/process_metadata.R"

rule infer_cnv_ab1:
    input:
        "results/cell_labelling_ab1/ab1_combined_identified.rds",
    output:
        directory("results/cell_labelling_ab1/ab1_infer_cnv"),
        "results/cell_labelling_ab1/infer_cnv_ab1_object.rds"
    conda:
        "../envs/bioconductor.yaml"
    script:
        "../scripts/infercnv_analysis.R"

rule DE_analysis_ab1:
    input:
        "results/cell_labelling_ab1/ab1_combined_identified.rds",
    conda:
        "../envs/bioconductor.yaml"
    output:
        "results/cell_labelling_ab1/DE_object_ab1.rds",
        "results/cell_labelling_ab1/metascape_table_ab1.csv"
    script:
        "../scripts/DE_list.R"
