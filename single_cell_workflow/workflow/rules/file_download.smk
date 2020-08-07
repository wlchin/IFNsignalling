

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_count_matrices:
    input:
        HTTP.remote("https://cloudstor.aarnet.edu.au/plus/s/xFJji4G6w2Juh8A/download", keep_local=True)
    cache: True
    output:
        countdir = directory("resources/count_data")
    run:
        outputName = "count_data.tar.gz"
        destName = "count_data"
        shell("mv {input} {outputName}; tar -xzvf {outputName}; mv {destName} {output.countdir}")

rule get_loomfiles:
    input:
        HTTP.remote("https://cloudstor.aarnet.edu.au/plus/s/jLlIJeEo0HK3Hsj/download", keep_local=True)
    cache: True
    output:
        loomdir = directory("resources/loomfiles")
    run:
        outputName = "loomfiles.tar.gz"
        destName = "loomfiles"
        shell("mv {input} {outputName}; tar -xzvf {outputName}; mv {destName} {output.loomdir}")

rule get_network_analysis_databases:
    input:
        HTTP.remote("https://cloudstor.aarnet.edu.au/plus/s/G7vyld8gurfAHLy/download", keep_local=True)
    cache: True
    output:
        networkdir = directory("resources/network_analysis")
    run:
        outputName = "network_analysis.tar.gz"
        destName = "network_analysis"
        shell("mv {input} {outputName}; tar -xzvf {outputName}; mv {destName} {output.networkdir}")

