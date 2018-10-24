#!/bin/bash

genbank=$1
biosample_xml="$genbank/metadata/biosample.xml"
biosample_xtract="$genbank/metadata/biosample.txt"

mkdir "$genbank/metadata"
esearch -db biosample  -query 'bacteria[orgn] AND biosample_assembly[filter]' | \
    efetch -format docsum > "$biosample_xml"

xtract -input "$biosample_xml" -pattern DocumentSummary/* \
    -group SampleData/BioSample \
        -if Attribute@harmonized_name \
            -lbl "accession" -element @accession \
    -group Attribute  \
        -if Attribute@harmonized_name \
            -def "missing" -element @harmonized_name Attribute > "$biosample_xtract"
