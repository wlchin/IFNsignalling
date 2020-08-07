        
rule subset_matrix_and_get_groups:
    input:
        "testdata/testdata_ab1_seurat314.rds",
        #"results/cell_labelling_ab1/ab1_combined_identified.rds",
    output:
        "results/scde/clean_matrix.rds",
        "results/scde/groups.rds"
    script:
        "../scripts/scde_filter_cells_by_library_size.R"

rule create_knn:
    input:
        "results/scde/clean_matrix.rds",
        "results/scde/groups.rds"
    cache: True
    output:
        "results/scde/knn.rds"
    script:
        "../scripts/scde_get_knn.R"

rule calculate_variance:
    input:
        "results/scde/clean_matrix.rds",
        "results/scde/knn.rds"
    output:
        "results/scde/varinfo.rds",
    script:
        "../scripts/scde_calculate_variance.R"

rule create_go_env:
    input:
        "results/scde/clean_matrix.rds",
        "resources/metadata/fast_isg.rds"
    output:
        "results/scde/go_env.rds",
        "results/scde/go_env_custom.rds"
    script:
        "../scripts/scde_goterms.R"
        
rule check_enrichment:
    input:
        "results/scde/go_env.rds",
        "results/scde/go_env_custom.rds",
        "results/scde/varinfo.rds",
        "testdata/testdata_ab1.rds",
        #"results/ab1_combined_identified.rds",
    output:
        "results/scde/enrichment.csv"
    script:
        "../scripts/scde_create_dataframe_for_seaborn.R"

rule plot_figures:
    input:
        "results/scde/enrichment.csv"
    output:
        "results/scde/typeIproduction_UMAP.png",
        "results/scde/typeIIproduction_UMAP.png",
        "results/scde/fastISG_UMAP.png",
        "results/scde/typeIproduction_violin.png",
        "results/scde/typeIIproduction_violin.png",
        "results/scde/fastISG_violin.png"
    notebook:
        "../notebooks/plot_scde_results.ipynb"
