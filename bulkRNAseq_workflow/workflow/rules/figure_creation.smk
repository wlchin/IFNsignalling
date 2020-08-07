
rule get_scaled_counts_matrices_for_fig1_fig2_fig4:
    input:
        "results/DE/{sample}_phenodata.rds",
	"results/mats/{sample}_respmat.rds",
        "results/mats/{sample}_nonrespmat.rds"
    output:
        "results/figures/{sample}_scaled_counts.rds",
        "results/figures/{sample}_scaled_individual_counts.rds",
    script:
        "../scripts/generate_scaled_averages_per_timepoint.R"
        
rule get_figure1:
    input:
        "results/figures/renca_ranks.rds",
        "results/figures/ab1_ranks.rds",
        "results/figures/ab1_scaled_individual_counts.rds",
        "results/figures/renca_scaled_individual_counts.rds"
    output:
        "results/figures/figure1.pdf",
    conda:
        "../env/visualisation.yml"
    script:
        "../scripts/fig1_heatmap.R"

rule get_figure2c:
    input:
        "results/figures/ab1_ISG_edgeweights.rds",
        "results/figures/ab1_scaled_counts.rds",
        "results/figures/renca_ISG_edgeweights.rds",
        "results/figures/renca_scaled_counts.rds"
    output:
        "results/figures/isg_figure2.pdf",
    conda:
        "../env/visualisation.yml"
    script:
        "../scripts/fig2_ISG_heatmap.R"

rule get_figure2_igraph:
    input:
        "results/figures/{sample}_resp_filt.rds",
	"results/figures/{sample}_scaled_counts.rds"
    output:
        "results/figures/vertices_{sample}.rds",
        "results/figures/network_figure2c_{sample}.pdf"
    conda:
        "../env/visualisation.yml"
    script:
        "../scripts/fig2_igraph_plots.R"

rule get_hiveplot:
    input:
        "results/figures/{sample}_resp_filt.rds",
        "workflow/scripts/mod_HPD.R"
    output:
        "results/figures/hiveplot_{sample}.pdf",
    conda:
        "../env/visualisation.yml"
    script:
        "../scripts/fig2_hiveplot.R"


