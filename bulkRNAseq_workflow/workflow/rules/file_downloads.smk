
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_jaspar_bedfile:
    input:
        HTTP.remote("cloudstor.aarnet.edu.au/plus/s/WOii3MUHN8tFqWe/download", keep_local=True)
    output:
        "resources/filteredJASPAR.bed"
    cache: True
    run:
        outputName = "filteredJASPAR.bed"
        shell("mv {input} {outputName}; mv {outputName} {output}")

rule get_index_for_kallisto:
    output:
        "resources/gencode.vM21.transcripts.fa.gz"
    shell:
        "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.transcripts.fa.gz"
