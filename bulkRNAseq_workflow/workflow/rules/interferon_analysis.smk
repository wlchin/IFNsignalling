
rule get_direct_network_with_kde_thresholds_for_fig2:
    input:
        "results/networks/ab1_direct_network_with_weights.rds",
	"results/networks/renca_direct_network_with_weights.rds"
    conda:
        "../env/network_creation.yml"
    output:
        "results/figures/cutoff.pdf",
	"results/figures/ab1_resp_filt.rds",
	"results/figures/renca_resp_filt.rds"
    script:
        "../scripts/filtering_and_thresholding.R"

rule get_direct_TF_to_ISG_edges_with_weights_for_fig2:
    input:
        "results/networks/ab1_direct_network_with_weights.rds",
        "results/networks/renca_direct_network_with_weights.rds",
        "results/networks/ab1_resp_GENIE3_full_network.rds",
        "results/networks/renca_resp_GENIE3_full_network.rds",
    conda:
        "../env/network_creation.yml"
    output:
        "results/figures/renca_ISG_edgeweights.rds",
        "results/figures/ab1_ISG_edgeweights.rds"
    script:
        "../scripts/merge_ISGs_from_two_networks.R"

rule get_regulator_rankings_for_fig1:
    input:
        "results/networks/{test}_resp_GENIE3_full_network.rds"
    conda:
        "../env/network_creation.yml"
    output:
        "results/figures/{test}_ranks.rds",
    script:
        "../scripts/regulator_rankings.R"
