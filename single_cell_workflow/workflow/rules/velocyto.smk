

rule loom_file_processing_ab1:
    input:
        rules.get_loomfiles.output.loomdir,
    params:
        loomfilenames = expand("resources/loomfiles/{samp}.loom", samp = AB1)
    conda:
        "../envs/velocyto.yaml"
    output:
        "results/velocyto/combined_loomfiles_ab1.loom"
    script:
        "../scripts/velocyto_combine_files.py"

rule loom_file_processing_renca:
    input:
        rules.get_loomfiles.output.loomdir,
    params:
        loomfilenames = expand("resources/loomfiles/{samp}.loom", samp = Renca)
    output:
        "results/velocyto/combined_loomfiles_renca.loom"
    conda:
        "../envs/velocyto.yaml"
    script:
        "../scripts/velocyto_combine_files.py"

rule process_cellmetadata_from_seurat_object:
    input:
        "results/test_combined_seurat_objects/{samples}_combined.rds"
    output:
        "results/velocyto/{samples}.csv"
    script:
        "../scripts/velocyto_create_cell_index.R"

rule momentum_and_embedding_individual_strains:
    input:
        "results/velocyto/combined_loomfiles_{sample}.loom",
        "results/velocyto/{sample}.csv"
    threads: 12
    conda:
        "../envs/velocyto.yaml"
    output:
        "results/velocyto/all_combined_{sample}.hdf5"
    script:
        "../scripts/velocyto_embedding_calculations.py"

rule create_velocity_plots_ab1:
    input:
        "results/velocyto/all_combined_ab1.hdf5",
    conda:
        "../envs/velocyto.yaml"
    output:
        "results/velocyto/small_velocity.png",
        "results/velocyto/quiver_velocity.png",
        "results/velocyto/Irf1_velocity.png"
    script:
        "../scripts/velocyto_AB1_diagrams.py"

rule interferon_analysis_on_ab1:
    input:
        "results/velocyto/all_combined_ab1.hdf5",
        "resources/metadata/fastgenes.tsv",
        "results/velocyto/ab1.csv"
    conda:
        "../envs/velocyto.yaml"
    output:
        "results/velocyto/ab1_clustermap.png",
        "results/velocyto/ab1_isg_velocities.png"
    script:
        "../scripts/velocyto_interferon_analysis.py"

rule identify_ab1_and_renca_monocytes:
    input:
        "results/test_combined_seurat_objects/ab1_combined.rds",
        "results/test_combined_seurat_objects/renca_and_ab1_combined.rds"
    conda:
        "../envs/velocyto.yaml"
    output:
        "results/velocyto/cluster_mono.csv",
        "results/velocyto/Ly6c.jpeg",
    script:
        "../scripts/velocyto_detect_monocytes_UMAP.R"

rule momentum_and_embedding_calc_on_common_UMAP:
    input:
        "results/velocyto/combined_loomfiles_{sample}.loom",
        "results/velocyto/renca_and_ab1.csv"
    threads: 12
    conda:
        "../envs/velocyto.yaml"
    output:
        "results/velocyto/common_embedding_combined_{sample}.hdf5"
    script:
        "../scripts/velocyto_embedding_calculations.py"

rule momentum_on_ab1_and_renca:
    input:
        "results/velocyto/common_embedding_combined_ab1.hdf5",
        "results/velocyto/common_embedding_combined_renca.hdf5",
        "results/velocyto/renca_and_ab1.csv",
        "results/velocyto/cluster_mono.csv"
    conda:
        "../envs/velocyto.yaml"
    output:
        "results/velocyto/facet_by_response.png",
        "results/velocyto/facet_by_strain.png"
    script:
        "../scripts/velocyto_combined_ab1_renca_analysis.py"
