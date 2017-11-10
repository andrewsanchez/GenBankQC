.. image:: https://travis-ci.org/andrewsanchez/genbank-qc.svg?branch=master

=============================================
           Genbank Quality Control
=============================================

Automatic quality control of genomes downloaded from public repositories such as the `National Center for Biotechnology Information's`_ GenBank database.

.. _National Center for Biotechnology Information's: https://www.ncbi.nlm.nih.gov/ 

- Assess the integrity of your FASTA collection

- Flag potential outliers to exlcude them from polluting your pipelines

- Labelling/annotation-independent quality control

  -  statistical quality metrics

  - `Genome distance estimation using MASH <http://mash.readthedocs.io/en/latest/>`_
  

====================
    Installation
====================

.. _ETE Toolkit: http://etetoolkit.org/ 

If you don't yet have a functional conda environment, please download and install `Miniconda`_.

.. _Miniconda: https://conda.io/miniconda.html

.. code:: bash

    conda create -n genbank-qc -c etetoolkit -c biocore pip ete3 scikit-bio
    source activate genbank-qc
    pip install genbank-qc
