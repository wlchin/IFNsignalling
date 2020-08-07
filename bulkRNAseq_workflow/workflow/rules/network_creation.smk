
rule create_GENIE3_network:
    input:
        "results/mats/{sample}_respmat.rds",
    cache:
        True
    conda:
        "../env/network_creation.yml"
    threads: 8
    output:
        "results/networks/{sample}_resp_GENIE3_full_network.rds",
    script:
        "../scripts/Genie3run.R"

rule get_TSS_site_for_each_DE_gene_using_UCSC_reference:
    input:
        "results/DE/{sample}_DE_list.txt", 
        "resources/tss.bed"
    conda:
        "../env/network_creation.yml"
    output:
        "results/bedfiles/{sample}_TSS_sites.bed"
    shell:
        "grep -w -f {input} > {output}"

rule get_direct_interactions_between_each_DEgene_and_TFs:
    input:
        jasparTFBSbed = "resources/filteredJASPAR.bed",
        TSSbedfile = "results/bedfiles/{sample}_TSS_sites.bed"
    conda:
        "../env/network_creation.yml"
    output:
        "results/bedfiles/Intersection_bedfile_{sample}.bed"
    shell:
        "bedtools window -a {input.TSSbedfile} -b {input.jasparTFBSbed} -l 400 -r 300 -sw  > {output}"

rule get_direct_network_from_JASPAR_and_DE_data:
    input:
        "results/bedfiles/Intersection_bedfile_{sample}.bed"
    conda:
        "../env/network_creation.yml"
    output:
        "results/networks/{sample}_direct_network.rds"
    script:
        "../scripts/create_direct_network_from_bed.R"

rule transfer_GENIE3_weights_to_direct_network:
    input:
        "results/networks/{sample}_resp_GENIE3_full_network.rds",
        "results/networks/{sample}_direct_network.rds"
    conda:
        "../env/network_creation.yml"
    output:
        "results/networks/{sample}_direct_network_with_weights.rds"
    script:
        "../scripts/populate_direct_network_with_weights.R"
