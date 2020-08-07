
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

rule get_count_matrices:
    input:
        FTP.remote("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130157/suppl/GSE130157_14546.14554.14562.RawCounts.txt.gz", keep_local=True)
    output:
        "resources/GSE130157_14546.14554.14562.RawCounts.txt"
    run:
        outputName = "GSE130157_14546.14554.14562.RawCounts.txt.gz"
        destName = "GSE130157_14546.14554.14562.RawCounts.txt"
        shell("mv {input} {outputName}; gunzip {outputName}; mv {destName} {output}")

rule get_count_matrices2:
    input:
        FTP.remote("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130157/suppl/GSE130157_15424R.15435R.RawCounts.txt.gz", keep_local=True)
    output:
        "resources/GSE130157_15424R.15435R.RawCounts.txt"
    run:
        outputName = "GSE130157_15424R.15435R.RawCounts.txt.gz"
        destName = "GSE130157_15424R.15435R.RawCounts.txt"
        shell("mv {input} {outputName}; gunzip {outputName}; mv {destName} {output}")

rule get_cistargetDB:
    input:
        HTTP.remote("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-5kb-7species.mc9nr.feather", keep_local=True)
    output:
        "resources/hg19-tss-centered-5kb-10species.mc9nr.feather"
    run:
        outputName = "hg19-tss-centered-5kb-10species.mc9nr.feather"
        shell("mv {input} {outputName}; mv {outputName} {output}")

rule get_tbl:
    input:
        HTTP.remote("https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl", keep_local=True)
    output:
        "resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
    run:
        outputName = "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
        shell("mv {input} {outputName}; mv {outputName} {output}")
