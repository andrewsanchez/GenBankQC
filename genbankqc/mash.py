"""
Functions for calling MASH with ``subprocess``.

`MASH`_ is a tool for fast genome and metagenome distance estimation using
`MinHash`_.  Make sure you have a MASH executable in your ``$PATH``.  Please
download the appropriate executable for your machine `here`_.

.. _MASH: https://mash.readthedocs.io/en/latest/
.. _here: https://github.com/marbl/Mash/releases
.. _MinHash: https://en.wikipedia.org/wiki/MinHash
"""

import os
import re
import subprocess
import pandas as pd


def generate_sketch_command(genome_path):
    """
    Return MASH sketch command corresponding to the full path of
    ``genome_path``.
    """
    basename = os.path.splitext(genome_path)[0]
    sketch_file = "{}.msh".format(basename)
    sketch_cmd = "mash sketch '{}' -o '{}'".format(genome_path, sketch_file)
    return sketch_cmd


def sketch_genome(genome_path):
    """
    :param genome_path: the full path to a FASTA

    Produce a MASH sketch file for ``genome``
    """
    sketch_cmd = generate_sketch_command(genome_path)
    subprocess.Popen(
        sketch_cmd, shell="True", stdout=subprocess.DEVNULL).wait()


def sketch_dir(directory):
    """
    Produce MASH sketch files for each FASTA in ``directory``
    """

    genome_paths = find_all_genome_paths(directory)
    for genome_path in genome_paths:
        sketch_genome(genome_path)


def paste(species_dir):
    """
    :param species_dir: path a directory containg FASTAs for
    a species or group of related genomes.

    Generate a master sketch file representing all genomes in
    ``species_dir``
    """
    paste_file = os.path.join(species_dir, 'all.msh')
    remove_old(paste_file)
    sketches = os.path.join(species_dir, "GCA*msh")
    paste_cmd = "mash paste {} {}".format(paste_file, sketches)
    subprocess.Popen(paste_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    return paste_file


def dist(species_dir):
    """Generate a MASH distance matrix, ``dmx.txt``

    :param species_dir: genome directory
    :returns: Pandas DataFrame representation of MASH distance matrix
    :rtype: pd.DataFrame

    """
    dmx_path = os.path.join(species_dir, 'dmx.txt')
    paste_file = os.path.join(species_dir, 'all.msh')
    remove_old(dmx_path)
    dist_cmd = "mash dist -t '{}' '{}' > '{}'".format(paste_file, paste_file,
                                                      dmx_path)
    subprocess.Popen(dist_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    dmx = format_dmx(dmx_path)
    return dmx


def mash(species_dir):
    """Convenience function to run the complete MASH workflow:
    * Create MASH sketch files
    * Create master sketch file, ``all.msh`` for all genomes
    * Generate distance matrix
    """
    sketch_dir(species_dir)
    paste(species_dir)
    dmx = dist(species_dir)
    return dmx


def format_dmx(dmx_path):
    """
    Set indices and headers of distance matrix to genome accession ID's
    instead of complete path to genome for readability.
    """
    dmx = pd.read_csv(dmx_path, index_col=0, sep="\t")
    genome_ids = []
    for i in dmx.index:
        genome_id = re.match('.*(GCA_\d+\.\d.*)(.fasta)', i).group(1)
        genome_ids.append(genome_id)
    dmx.index = genome_ids
    dmx.columns = genome_ids
    dmx.to_csv(dmx_path, sep="\t")
    return dmx


def remove_old(f):
    if os.path.isfile(f):
        os.remove(f)


def find_all_genome_paths(directory):
    """
    Return full paths to all FASTAs in ``directory``.
    """
    genome_paths = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.endswith('fasta'):
                genome_path = os.path.join(root, f)
                genome_paths.append(genome_path)
    return genome_paths


def find_all_sketches(genbank):
    """
    Return full paths to all MASH sketch files in ``directory``.
    """
    sketches = []
    for root, dirs, files in os.walk(genbank):
        for f in files:
            if re.match('GCA.*msh', f):
                sketch = os.path.join(root, f)
                sketches.append(sketch)
    return sketches
