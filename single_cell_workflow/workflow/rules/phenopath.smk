               
rule make_counts:
    input:
        "results/cell_labelling_ab1/ab1_combined_identified.rds",
    output:
        "results/phenopath/expression_matrix.rds",
        "results/phenopath/phenovector.rds",
    script:
        "../scripts/phenopath_get_highvar.R"

rule run_phenopath:
    input:
        "results/phenopath/expression_matrix.rds",
        "results/phenopath/phenovector.rds"
    output:
        "results/phenopath/run_output.rds"
    script:
        "../scripts/phenopath_run.R"
        
rule process_phenopath_run:
    input:
        "results/phenopath/run_output.rds",
    output:
        "results/phenopath/interaction_output.txt",
        "results/phenopath/top_interactions.jpeg",
        "results/phenopath/chiplot.jpeg",
        "results/phenopath/trajectory_info.rds"
    script:
        "../scripts/phenopath_interaction_analysis.R"
        
rule extract_significant_genes:
    input:
        "results/phenopath/interaction_output.txt",
    output:
        "results/phenopath/responder_genes.txt",
        "results/phenopath/nonresponder_genes.txt"
    script:
        "../scripts/phenopath_significant_genes.R"

