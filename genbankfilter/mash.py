#!/usr/bin/env python

import os
import re
import argparse
import glob
import subprocess

def write_sketch_commands(genbank_mirror, assembly_summary, new_genomes):

    sketch_commands = os.path.join(genbank_mirror, ".info", "sketch_commands.txt")
    if os.path.isfile(sketch_commands):
        os.remove(sketch_commands)

    with open(sketch_commands, 'a') as cmds:
        for genome in new_genomes:
            species_dir = assembly_summary.scientific_name.loc[genome]
            # TODO: This is no good.  Make it more efficient.
            fasta = os.path.join(genbank_mirror, species_dir, "{}*fasta".format(genome))
            fasta = glob.glob(fasta)

            try:
                fasta = fasta[0]
            except IndexError:
                continue

            sketch_dst = os.path.join(genbank_mirror, species_dir, "{}.msh".format(genome))
            cmd = "/common/contrib/bin/mash-Linux64-v1.1.1/mash sketch {} -o {}\n".format(fasta, sketch_dst)
            cmds.write(cmd)

    return sketch_commands # get the line count of this file

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
    master_sketch = os.path.join(species_dir, 'all.msh')
    if os.path.isfile(master_sketch):
        os.remove(master_sketch)
    all_sketchs = os.path.join(species_dir,'sketches.txt')
    paste_cmd = "mash paste -l '{}' '{}'".format(master_sketch, all_sketchs)
    subprocess.Popen('ls {}/*msh >> {}'.format(species_dir, all_sketchs), shell='True').wait()
    subprocess.Popen(paste_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    return master_sketch

def dist(genbank_mirror, assembly_summary, species):

    species_dir = os.path.join(genbank_mirror, species)
    dist_matrix = os.path.join(species_dir, 'all_dist.msh')
    all_msh = os.path.join(species_dir,'all.msh')
    if os.path.isfile(dist_matrix):
        os.remove(dist_matrix)
    dist_cmd = "mash dist -t '{}' '{}' > '{}'".format(all_msh, all_msh, dist_matrix)
    subprocess.Popen(dist_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
    return dist_matrix
