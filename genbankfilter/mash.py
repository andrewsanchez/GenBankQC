#!/usr/bin/env python

import os
import re
import argparse
import glob
import subprocess
import pandas as pd

def sketch(genbank_mirror, assembly_summary, genome):

    # TODO: Using os.walk will avoid having to look up species_dir using assembly_summary
    species_dir = assembly_summary.scientific_name.loc[genome]
    fasta = os.path.join(genbank_mirror, species_dir, "{}*fasta".format(genome))
    fasta = glob.glob(fasta)
    try:
        fasta = fasta[0]
    except IndexError:
        pass
    sketch_file = os.path.join(genbank_mirror, species_dir, "{}.msh".format(genome))
    sketch_cmd = "mash sketch '{}' -o '{}'".format(fasta, sketch_file)
    subprocess.Popen(sketch_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    return sketch_file

def paste(genbank_mirror, assembly_summary, species):

    species_dir = os.path.join(genbank_mirror, species)
    all_msh = os.path.join(species_dir, 'all.msh')
    remove_old(all_msh)
    all_sketchs = os.path.join(species_dir,'sketches.txt')
    paste_cmd = "mash paste -l '{}' '{}'".format(all_msh, all_sketchs)
    subprocess.Popen('ls {}/*msh >> {}'.format(species_dir, all_sketchs), shell='True').wait()
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
