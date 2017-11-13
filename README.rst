.. image:: https://travis-ci.org/andrewsanchez/genbank-qc.svg?branch=master

=============================================
           GenBank Quality Control
=============================================

Please see http://genbank-qc.readthedocs.io/en/latest/ for complete documentation.

GenBankQC is an effort to address the quality control problem for public databases such as the National Center for Biotechnology Information's `GenBank`_.  The goal is to offer a simple, efficient, and automated solution for assessing the quality of your genomes.

Features
--------

- Flag potential outliers to exclude them from polluting your pipelines

- Labelling/annotation-independent quality control based on:

  -  Simple metrics

  - Genome distance estimation using `MASH`_

The genbankqc work-flow consists of the following steps:

#. Generate statistics for each genome based on the following metrics:

   * Number of unknown bases
   * Number of contigs
   * Assembly size
   * Average `MASH`_ distance compared to other genomes

#. Flag potential outliers based on these statistics:

   * Flag genomes containing more than a certain number of unknown bases.

   * Flag genomes outside of a range based on the median absolute deviation.

     * Applies to number of contigs and assembly size

   * Flag genomes whose `MASH`_ distance is greater than the upper end of the median absolute deviation.

#. Visualize the results with a color coded tree

Usage
-----

.. code:: bash

          genbankqc /path/to/genomes

          open /path/to/genomes/Escherichia_coli/qc/200_3.0_3.0_3.0/tree.svg


Installation
------------

If you don't yet have a functional conda environment, please download and install `Miniconda`_.

.. code::

    conda create -n genbankqc -c etetoolkit -c biocore pip ete3 scikit-bio

    source activate genbankqc

    pip install genbankqc


.. _GenBank: https://www.ncbi.nlm.nih.gov/genbank/
.. _ETE Toolkit: http://etetoolkit.org/ 
.. _Miniconda: https://conda.io/miniconda.html
.. _MASH: http://mash.readthedocs.io/en/latest/
