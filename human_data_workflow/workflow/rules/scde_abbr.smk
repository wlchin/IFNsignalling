rule create_counts_scde:
    input:
         "results/counts_monocytes.rds"
    output:
        "results/filtered_monocytes_scde.rds",
        "results/knn_mono.rds",
        "results/varinfo_mono.rds"
    container:
        "docker://wlc27/scde_with_db:v2"
    script:
        "scripts/run_scde_analysis.R"

rule create_go_env:
    input:
        "results/filtered_monocytes_scde.rds"
        "resources/fast_isg.rds"
    output:
        "results/scde/go_env.rds",
        "results/scde/go_env_custom.rds"
    container:
        "docker://wlc27/scde_with_db:v2"
    script:
        "../scripts/go_terms.R"

rule preparedf:
    input:
        "results/scde/go_env.rds",
        "results/scde/go_env_custom.rds",
        "resources/varinfo_mono.rds",
        "results/metadata_for_scde.rds"
    output:
        "results/scde/df_comb.rds",
        "results/scde.tsv"
    container:
        "docker://wlc27/scde_with_db:v2"
    script:
        "../scripts/df_from_go_env.R"
