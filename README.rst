.. image:: https://travis-ci.org/andrewsanchez/genbank-qc.svg?branch=master

================================
           Genbank-QC
================================

A set of tools for curating your genomes from the `National Center for Biotechnology Information's`_ public database.

.. _National Center for Biotechnology Information's: https://www.ncbi.nlm.nih.gov/ 

- Assess the integrity of your FASTA collection

- Labelling/annotation-independent quality control using

  -  statistical quality metrics

  - `Genome and metagenome distance estimation using MinHash <http://mash.readthedocs.io/en/latest/>`_
  

====================
    Installation
====================

While I have not yet created a conda recipe for genbank-qc, the recommended installation method is to use conda.  Genbank-qc uses the `ETE Toolkit`_ which can conflict with other dependencies due to the version of python it requires.  Anaconda makes it easy to create a virtual environment with a specific version of python.

.. _ETE Toolkit: http://etetoolkit.org/ 

If you don't yet have a working conda installation, please download and install `Miniconda`_.

.. _Miniconda: https://conda.io/miniconda.html

.. code:: bash

    conda create --name genbank-qc python=3.4
    pip install genbank-qc
