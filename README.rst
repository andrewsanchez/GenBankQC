.. image:: https://travis-ci.org/andrewsanchez/genbank-qc.svg?branch=master

================================
           Genbank-QC
================================

A set of tools for curating your genomes from the `National Center for Biotechnology Information's`_ public database.

.. _National Center for Biotechnology Information's: https://www.ncbi.nlm.nih.gov/ 

- Assess the integrity of your FASTA collection

- Labelling/annotation-independent quality control

  -  statistical quality metrics

  - `Genome and metagenome distance estimation using MinHash <http://mash.readthedocs.io/en/latest/>`_
  

====================
    Installation
====================

.. _ETE Toolkit: http://etetoolkit.org/ 

If you don't yet have a working conda installation, please download and install `Miniconda`_.

.. _Miniconda: https://conda.io/miniconda.html

.. code:: bash

    conda create -n genbank-qc -c etetoolkit -c biocore pip ete3 scikit-bio
    source activate genbank-qc
    pip install genbank-qc
