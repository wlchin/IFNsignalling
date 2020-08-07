
import loompy

print(snakemake.params.loomfilenames)

loompy.combine(snakemake.params.loomfilenames, snakemake.output[0], key="Accession")
