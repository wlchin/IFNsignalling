
rule single_r_label_renca:
    input:
        "results/combined_seurat_objects/renca_combined.rds"
    conda:
        "../envs/bioconductor.yaml"
    output:
        "results/cell_labelling_renca/renca_single_r_labels.rds",
        "results/cell_labelling_renca/renca_labels.jpeg",
    script:
        "../scripts/singler_cell_labelling.R"

rule tumour_predictions_renca:
    input:
        "results/combined_seurat_objects/renca_combined.rds",
        "results/combined_seurat_objects/Park.rds"
    conda:
        "../envs/bioconductor.yaml"
    output:
        "results/cell_labelling_renca/renca_metadata_predictions.rds"
    script:
        "../scripts/detect_tumour_kidney_clusters.R"

rule analyse_metadata_renca:
    input:
        "results/combined_seurat_objects/renca_combined.rds",
        "results/cell_labelling_renca/renca_single_r_labels.rds",
        "results/cell_labelling_renca/renca_metadata_predictions.rds"
    conda:
        "../envs/bioconductor.yaml"
    output:
        "results/cell_labelling_renca/renca_combined_identified.rds",
        "results/cell_labelling_renca/tumour_ident_and_singleR_renca.jpeg"
    script:
        "../scripts/process_metadata.R"

rule infer_cnv_renca:
    input:
        "results/cell_labelling_renca/renca_combined_identified.rds",
    conda:
        "../envs/bioconductor.yaml"
    output:
        directory("results/cell_labelling_renca/renca_infer_cnv"),
        "results/cell_labelling_renca/infer_cnv_renca_object.rds"
    script:
        "../scripts/infercnv_analysis.R"

rule DE_analysis_renca:
    input:
        "results/cell_labelling_renca/renca_combined_identified.rds",
    conda:
        "../envs/bioconductor.yaml"
    output:
        "results/cell_labelling_renca/DE_object_renca.rds",
        "results/cell_labelling_renca/metascape_table_renca.csv"
    script:
        "../scripts/DE_list.R"
