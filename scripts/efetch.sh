#!/bin/bash

metadata_dir=$1
biosample_xml="$metadata_dir/biosample.xml"
biosample_xtract="$metadata_dir/biosample.txt"

# Get SRA runs for each ID
for f in ${metadata_dir}/sra_ids/sra*txt;
do epost -db sra -input $f -format acc | \
        efetch -format docsum | \
        xtract -pattern DocumentSummary/* \
               -group Biosample -ret '\t' -element Biosample \
               -group Runs -sep ',' -element Runs/Run@acc \
               >> ${metadata_dir}/sra_runs.txt
done

# This data is retrieved using Biopython
# esearch -db biosample  -query 'bacteria[orgn] AND biosample_assembly[filter]' | \
#     efetch -format docsum > "$biosample_xml"

# This data is retrieved using Biopython
# xtract -input "$biosample_xml" -pattern DocumentSummary/* \
#     -group SampleData/BioSample \
#         -if Attribute@harmonized_name \
#             -lbl "accession" -element @accession \
#     -group Attribute  \
#         -if Attribute@harmonized_name \
#             -def "missing" -element @harmonized_name Attribute > "$biosample_xtract"
