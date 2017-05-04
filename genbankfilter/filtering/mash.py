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

def sketch(genbank_mirror, assembly_summary, missing_sketch_files, logger):

    for genome in missing_sketch_files:
        species_dir = assembly_summary.scientific_name.loc[genome]
        fasta = os.path.join(genbank_mirror, species_dir, "{}*fasta".format(genome))
        fasta = glob.glob(fasta)
        try:
            fasta = fasta[0]
        except IndexError:
            logger.info("IndexError for {} during sketch".format(genome))
            continue
        sketch_dst = os.path.join(genbank_mirror, species_dir, "{}.msh".format(genome))
        sketch_cmd = "mash sketch '{}' -o '{}'".format(fasta, sketch_dst)
        subprocess.Popen(sketch_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
        logger.info("Created sketch file for {}".format(genome))

def paste(genbank_mirror, assembly_summary, species, logger):

    def paste_cmd(species):

        for name in species:
            path = os.path.join(genbank_mirror, name)
            master_sketch = os.path.join(path, 'all.msh')
            if os.path.isfile(master_sketch):
                os.remove(master_sketch)
                logger.info('Old Master sketch file for {} removed'.format(name))
            master_sketch = os.path.join(path, 'all.msh')
            all_sketch_files = os.path.join(path,'*msh')
            paste_cmd = "mash paste '{}' '{}'".format(master_sketch, all_sketch_files)
            subprocess.Popen(paste_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
            logger.info('Master sketch fille for {} created'.format(name))

    paste_cmd(species)

def dist(genbank_mirror, assembly_summary, species, logger):

    def dist_cmd(species):

        for name in species:
            path = os.path.join(genbank_mirror, name)
            dist_matrix = os.path.join(path, 'all_dist.msh')
            all_msh = os.path.join(path,'all.msh')
            if os.path.isfile(dist_matrix):
                logger.info('Old distance matrix for {} removed'.format(name))
                os.remove(dist_matrix)
            dist_cmd = "mash dist -t '{}' '{}' > '{}'".format(all_msh, all_msh, dist_matrix)
            subprocess.Popen(dist_cmd, shell="True", stdout=subprocess.DEVNULL).wait()
            logger.info('Distance matrix for {} created'.format(name))

    dist_cmd(species)
