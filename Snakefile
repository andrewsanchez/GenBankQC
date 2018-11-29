rule biosample:
     input:
          {dir}
     output:
          {dir}/metadata/assembly_summary.txt
          {dir}/metadata/biosample.csv
          {dir}/metadata/sra_ids/sra_ids_1.csv
          {dir}/metadata/sra_ids/sra_ids_2.csv
          ... etc ...
     script:
          genbankqc metadata {dir}