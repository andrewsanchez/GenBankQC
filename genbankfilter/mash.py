#!/usr/bin/env python

import os
import re
import argparse
import glob
import subprocess

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
    dst_mx = os.path.join(species_dir, 'all_dist.msh')
    all_msh = os.path.join(species_dir,'all.msh')
    remove_old(dst_mx)
    dist_cmd = "mash dist -t '{}' '{}' > '{}'".format(all_msh, all_msh, dst_mx)
    subprocess.Popen(dist_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    return dst_mx

def remove_old(f):
    if os.path.isfile(f):
        os.remove(f)
