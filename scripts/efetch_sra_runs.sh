#!/bin/bash

metadata_dir=$1
outdir=${metadata_dir}/.sra_runs

# Get SRA runs for each ID
mkdir ${outdir}
split -l 5000 sra_ids.txt ${outdir}/sra_ids_
for f in ${outdir}/sra_ids*;
do epost -db sra -input $f -format acc | \
        efetch -format docsum | \
        xtract -pattern DocumentSummary/* \
               -group Biosample -ret '\t' -element Biosample \
               -group Runs -sep ',' -element Runs/Run@acc \
               >> ${metadata_dir}/sra_runs.tsv
done
rm ${outdir}
