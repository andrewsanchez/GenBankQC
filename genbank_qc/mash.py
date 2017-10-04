#!/usr/bin/env python

import os
import re
import subprocess
import pandas as pd


def generate_sketch_command(genome_path):
    """
    Return sketch command for genome_path.
    sketch_cmd is based on full path to genome_path.
    """

    basename = os.path.splitext(genome_path)[0]
    sketch_file = "{}.msh".format(basename)
    sketch_cmd = "mash sketch '{}' -o '{}'".format(genome_path, sketch_file)
    return sketch_cmd


def sketch_genome(genome_path):
    """
    Produce a sketch file for genome, where genome is the full path to a FASTA.
    """

    sketch_cmd = generate_sketch_command(genome_path)
    subprocess.Popen(
        sketch_cmd, shell="True", stdout=subprocess.DEVNULL).wait()


def sketch_dir(directory):
    """
    Produce sketch files for each FASTA in dir.
    """

    genome_paths = find_all_genome_paths(directory)
    for genome_path in genome_paths:
        sketch_genome(genome_path)


def paste(species_dir):
    """
    Generate a master sketch file representing all genomes in species_dir/*msh
    """
    paste_file = os.path.join(species_dir, 'all.msh')
    remove_old(paste_file)
    sketches = os.path.join(species_dir, "GCA*msh")
    paste_cmd = "mash paste {} {}".format(paste_file, sketches)
    subprocess.Popen(paste_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    return paste_file


def dist(species_dir):

    dmx_path = os.path.join(species_dir, 'dmx.txt')
    paste_file = os.path.join(species_dir, 'all.msh')
    remove_old(dmx_path)
    dist_cmd = "mash dist -t '{}' '{}' > '{}'".format(paste_file, paste_file,
                                                      dmx_path)
    subprocess.Popen(dist_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    dmx = format_dmx(dmx_path)
    return dmx


def mash(species_dir):
    sketch_dir(species_dir)
    paste(species_dir)
    dmx = dist(species_dir)
    return dmx


def format_dmx(dmx_path):
    """
    Set indices and headers to the accession ID's
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
    Return full paths to all FASTAs in `directory`.
    """
    genome_paths = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.endswith('fasta'):
                genome_path = os.path.join(root, f)
                genome_paths.append(genome_path)
    return genome_paths


def find_all_sketches(genbank):
    sketches = []
    for root, dirs, files in os.walk(genbank):
        for f in files:
            if re.match('GCA.*msh', f):
                sketch = os.path.join(root, f)
                sketches.append(sketch)
    return sketches
