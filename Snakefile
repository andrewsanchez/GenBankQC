from pathlib import Path

metadata_dir = "metadata"
sra_ids = [f for f in Path(metadata_dir).joinpath("sra_ids").glob("sra*txt")]

rule biosample:
    input:
        "{metadata_dir}"
    output:
         "{metadata_dir}/biosample.csv",
         "{metadata_dir}/assembly_summary.csv",
    script:
         "genbankqc metadata {metadata_dir}"

rule runs:
    input:
        "{sra_ids}"
    output:
        "{metadata_dir}/sra_runs.txt"
    shell:
        "sh efetch.sh {input}"

rule install:
    conda:
        "requirements/conda.yaml"
    shell:
        """
        pip install -r requirements/pip.txt &&
        pip install -r requirements/dev.txt &&
        pip install --no-deps -e .
        """
