.. image:: https://travis-ci.org/andrewsanchez/genbank-qc.svg?branch=master

================================
           Genbank-QC
================================

A set of tools for curating your genomes from the .. _National Center for Biotechnology Information's:https://www.ncbi.nlm.nih.gov/ public database.

- Assess the integrity of your FASTA collection
  - Labelling/annotation independent quality control using
    +  statistical quality metrics
  +  .. _genome and metagenome distance estimation using MinHash:http://mash.readthedocs.io/en/latest/

====================
    Installation
====================

While I have not yet created a conda recipe for genbank-qc, the recommended installation method is to use conda.
Genbank-qc uses the .. _ETE Toolkit:http://etetoolkit.org/ which has conflicted with other genbank-qc dependencies due to the version of python it needs to use.
Anaconda makes it easy to create a virtual environment with a specific version of python.

If you don't yet have a working conda installation, please download and install.. _Miniconda:https://conda.io/miniconda.html

    conda create --name genbank-qc python=3.4
    pip install genbank-qc
