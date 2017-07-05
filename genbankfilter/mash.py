#!/usr/bin/env python

import os
import re
import argparse
import glob
import subprocess
import pandas as pd

def generate_sketch_command(genome):

    """
    Return sketch command for genome.
    sketch_cmd is based on full path to genome.
    """

    basename = os.path.splitext(genome)[0]
    sketch_file = "{}.msh".format(basename)
    sketch_cmd = "mash sketch '{}' -o '{}'".format(genome, sketch_file)
    return sketch_cmd


def sketch_genome(genome):

    """
    Produce a sketch file for genome, where genome is the full path to a FASTA.
    """

    sketch_cmd = generate_sketch_command(genome)
    subprocess.Popen(sketch_cmd, shell="True", stdout=subprocess.DEVNULL).wait()

def paste(genbank_mirror, assembly_summary, species):

    species_dir = os.path.join(genbank_mirror, species)
    all_msh = os.path.join(species_dir, 'all.msh')
    remove_old(all_msh)
    all_sketchs = os.path.join(species_dir,'sketches.txt')
    paste_cmd = "mash paste -l '{}' '{}'".format(all_msh, all_sketchs)
    subprocess.Popen('ls {}/*msh >> {}'.format(species_dir, all_sketchs), shell='True').wait()
def sketch_species(species_dir):

    """
    Produce sketch files for each FASTA in species_dir.
    """

    fastas = (f for f in os.listdir(species_dir) if f.endswith('fasta'))
    for f in fastas:
        genome_path = os.path.join(species_dir, f)
        sketch_cmd = generate_sketch_command(genome)
        subprocess.Popen(sketch_cmd, shell="True", stdout=subprocess.DEVNULL).wait()

    subprocess.Popen(paste_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    return all_msh

def dist(genbank_mirror, assembly_summary, species):

    species_dir = os.path.join(genbank_mirror, species)
    dst_mx_path = os.path.join(species_dir, 'dst_mx.txt')
    all_msh = os.path.join(species_dir,'all.msh')
    remove_old(dst_mx_path)
    dist_cmd = "mash dist -t '{}' '{}' > '{}'".format(all_msh, all_msh, dst_mx_path)
    subprocess.Popen(dist_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    dst_mx = format_dst_mx(dst_mx_path)
    return dst_mx

def format_dst_mx(dst_mx_path):

    """
    Set indices and headers to the accession ID's
    """

    dst_mx = pd.read_csv(dst_mx_path, index_col=0, sep="\t")
    new_index = []
    new_columns = []
    for i in dst_mx.index:
        name = i.split("/")[-1].strip(".fasta")
        new_index.append(name)
        new_columns.append(name)

    dst_mx.index = new_index
    dst_mx.columns = new_columns
    dst_mx.to_csv(dst_mx_path, sep="\t")
    return dst_mx


def remove_old(f):
    if os.path.isfile(f):
        os.remove(f)
