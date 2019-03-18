.. image:: https://api.travis-ci.org/andrewsanchez/GenBankQC.svg?branch=master
.. image:: https://codecov.io/gh/andrewsanchez/GenBankQC/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/andrewsanchez/GenBankQC

=============================================
           GenBank Quality Control
=============================================

Complete documentation lives at `genbankqc.readthedocs.io`_.  It is a work in progress.

GenBankQC is an effort to address the quality control problem for public databases such as the National Center for Biotechnology Information's `GenBank`_.  The goal is to offer a simple, efficient, and automated solution for assessing the quality of your genomes.

Note
----

    Please note that GenbankQC is currently in alpha.  As a proof of concept for a specific use case, it currently has limitations that users should be aware of.  If there is interest, we will address the issues to make it more convenient to use.  Please see `caveats <#caveats>`__ for more details.


Features
--------

- Labelling/annotation-independent quality control based on:

  -  Simple metrics

  - Genome distance estimation using `MASH`_

- Flag potential outliers to exclude them from polluting your pipelines

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

::

    genbankqc /path/to/genomes
    open /path/to/genomes/Escherichia_coli/qc/200_3.0_3.0_3.0/tree.svg


Installation
------------

If you don't yet have a functional conda environment, please download and install `Miniconda`_.

.. code::

    conda create -n genbankqc -c etetoolkit -c biocore pip ete3 scikit-bio

    source activate genbankqc

    pip install genbankqc


.. _caveats:

Caveats
--------

There are some arbitrary, hard-coded limitations regarding file names.  This is because the project originally began as a part of the NCBI Tool Kit (`NCBITK`_) which we use for downloading genomes from NCBI.  NCBITK generates a specific directory structure and file naming scheme which GenbankQC currently expects.

If you'd like to use GenBankQC without using NCBITK, all that is required is that your file names match the python regular expression ``re.compile('.*(GCA_\d+\.\d.*)(.fasta)')``.  You can quickly test this by following my example at `pythex.org`_.

.. _pythex.org: https://pythex.org/?regex=.*(GCA_%5Cd%2B%5C.%5Cd.*)(.fasta)&test_string=GCA_002415405.1_Acinetobacter_nosocomialis_UBA5139_Scaffold.fasta&ignorecase=0&multiline=0&dotall=0&verbose=0

.. _NCBITK:  https://github.com/andrewsanchez/NCBITK
.. _GenBank: https://www.ncbi.nlm.nih.gov/genbank/
.. _ETE Toolkit: http://etetoolkit.org/
.. _Miniconda: https://conda.io/miniconda.html
.. _MASH: http://mash.readthedocs.io/en/latest/
.. _genbankqc.readthedocs.io: http://genbankqc.readthedocs.io/en/latest/

.. image:: https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square
           :target: https://yangsu.github.io/pull-request-tutorial/
