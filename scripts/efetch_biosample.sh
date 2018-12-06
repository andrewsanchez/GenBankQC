#!/bin/bash

This data is retrieved using Biopython
biosample_xml="$metadata_dir/biosample.xml"
biosample_xtract="$metadata_dir/biosample.txt"
esearch -db biosample  -query 'bacteria[orgn] AND biosample_assembly[filter]' | \
    efetch -format docsum > "$biosample_xml"

This data is retrieved using Biopython
xtract -input "$biosample_xml" -pattern DocumentSummary/* \
    -group SampleData/BioSample \
        -if Attribute@harmonized_name \
            -lbl "accession" -element @accession \
    -group Attribute  \
        -if Attribute@harmonized_name \
            -def "missing" -element @harmonized_name Attribute > "$biosample_xtract"
