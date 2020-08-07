

rule create_individual_seurat_objects:
    input:
        rules.get_count_matrices.output.countdir
    params:
        folder="{samp_name}",
        fpath = "resources/count_data/{samp_name}/filtered_feature_bc_matrix/"
    container:
        "docker://wlc27/singler_legacy:v2"
    output:
        "results/seurat_objects/{samp_name}.rds"
    script:
        "../scripts/create_seurat_object.R"

rule merge_ab1_objects:
    input:
        expand("results/seurat_objects/{samp}.rds", samp = AB1)
    output:
        "results/combined_seurat_objects/ab1_combined.rds"
    container:
        "docker://wlc27/singler_legacy:v2"
    script:
        "../scripts/sc_transform.R"

rule merge_renca_objects:
    input:
        expand("results/seurat_objects/{samp}.rds", samp = Renca)
    output:
        "results/combined_seurat_objects/renca_combined.rds"
    container:
        "docker://wlc27/singler_legacy:v2"
    script:
        "../scripts/sc_transform.R"

rule merge_all_ab1_and_renca:
    input:
        expand("results/seurat_objects/{samp}.rds", samp = allsamples)
    container:
        "docker://wlc27/singler_legacy:v2"
    output:
        "results/combined_seurat_objects/renca_and_ab1_combined.rds"
    script:
        "../scripts/sc_transform.R"

rule merge_ab1_tag_objects:
    input:
        expand("results/seurat_objects/{samp}.rds", samp = tagged)
    output:
        "results/combined_seurat_objects/ab1_tagged_combined.rds",
    container:
        "docker://wlc27/singler_legacy:v2"
    script:
        "../scripts/sc_transform.R"

rule create_park_et_al_object:
    input:
        "resources/count_data/GSE107585_datamatrix.txt"
    output:
        "results/combined_seurat_objects/Park.rds",
    container:
        "docker://wlc27/singler_legacy:v2"
    script:
        "../scripts/create_park_dataset.R"

