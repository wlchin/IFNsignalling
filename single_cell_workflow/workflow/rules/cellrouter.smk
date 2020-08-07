LIBDIR = "CellRouter"
        
rule create_cellrouter_object:
    input:
        "results/cell_labelling_ab1/ab1_combined_identified.rds",
    output:
        "results/cellrouter/cellrouter_obj.rda"
    conda:
        "../envs/R_env.yaml"
    script:
        "../scripts/cellrouter_create_object.R"

rule colour_with_markers:
    input:
        "results/ab1_combined_identified.rds",
    output:
        "results/cellrouter/markers.png",
        "results/cellrouter/markers2.png"
    conda:
        "../envs/R_envt.yaml"
    script:
        "../scripts/cellrouter_marker_colours.R"

rule diffusion_pseudotime:
    input:
        "results/cellrouter/cellrouter_obj.rda",
    output:
        "results/cellrouter/cellrouter_obj_with_diffusion.rda",
        temp("cell_edge_weighted_network.txt"),
        temp(directory("cluster_12.cluster_3")),
        temp("graph_clusters.gml"),
        temp("graph_subpopulations.gml"),
        temp("Cells_FlowNetwork_paths.ser")
    params:
        libdir = "CellRouter"
    container:
        "docker://wlc27/cellrouter_env" 
    script:
        "../scripts/cellrouter_pseudotime.R"

rule enrichment_of_diffusion_components:
    input:
        "results/cellrouter/cellrouter_obj_with_diffusion.rda"
    output:
        expand("results/cellrouter/cellrouter_component_{ind}.txt", ind = [1,2,3,4,5]),
        "results/cellrouter/diffusion_comp_list.rds",
        "results/cellrouter/dynamics_diffusion.pdf"
    conda:
        "../envs/R_env.yaml"
    script:
        "../scripts/cellrouter_extract_diffusion_components.R"

