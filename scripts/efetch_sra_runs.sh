#!/bin/bash

metadata_dir=$1
sra_ids=${metadata_dir}sra_ids.txt
sra_ids_split=${metadata_dir}_sra_ids_*
sra_runs=${metadata_dir}sra_runs.tsv

find $metadata_dir -name "*sra_ids*" -or -name "sra_runs.tsv" -delete

# Get SRA runs for each ID
split -l 5000 ${sra_ids} ${metadata_dir}_sra_ids_
for f in $sra_ids_split;
do epost -db sra -input $f -format acc | \
        efetch -format docsum | \
        xtract -pattern DocumentSummary/* \
               -group Biosample -ret '\t' -element Biosample \
               -group Runs -sep ',' -element Runs/Run@acc \
               >> $sra_runs
done
