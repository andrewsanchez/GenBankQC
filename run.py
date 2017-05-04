#!/usr/bin/env python

import os
import argparse
import logging

from genbankfilter import get_resources
from genbankfilter import config
from genbankfilter import curate
from genbankfilter import mash

def setup(genbank_mirror, species_list, update):

    path_vars = config.instantiate_path_vars(genbank_mirror)
    info_dir, slurm, out, logger = path_vars
    assembly_summary = get_resources.get_resources(genbank_mirror, logger, update)
    logger.info('{} genomes in assembly_summary.txt'.format(len(assembly_summary)))
    species = curate.get_species_list(assembly_summary, species_list)

    return path_vars, assembly_summary, species

def assess_genbank(genbank_mirror, assembly_summary, species_list, logger):

    # TODO: Why does this take so long?
    genbank_status = curate.assess_genbank_mirror(genbank_mirror, assembly_summary, species_list)
    local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_status

    logger.info("{} genomes present in local collection.".format(len(local_genomes)))
    logger.info("{} genomes missing from local collection.".format(len(new_genomes)))
    if len(new_genomes) == 0:
        logger.info("Local collection is up to do with latest assembly summary file.")
    logger.info('{} sketch files present in local collection.'.format(len(sketch_files)))
    logger.info('{} sketch files not in local collection.'.format(len(missing_sketch_files)))

    return genbank_status

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-s", "--species", help = 'List of species', nargs='+', default='all')
    parser.add_argument("-p", "--slurm", help = 'Submit jobs in parallel via SLURM pipeline tool.', action="store_true")
    parser.add_argument("--no_update", action='store_false')
    args = parser.parse_args()

    if args.no_update:
        update=False
    else:
        update=True

    genbank_mirror = args.genbank_mirror
    path_vars, assembly_summary, species = setup(genbank_mirror, args.species, update)
    info_dir, slurm, out, logger = path_vars
    genbank_status = assess_genbank(genbank_mirror, assembly_summary, species, logger)
    local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_status
    mash.sketch(genbank_mirror, assembly_summary, missing_sketch_files, logger)
    mash.paste(genbank_mirror, assembly_summary, species, logger)
    mash.dist(genbank_mirror, assembly_summary, species, logger)

if __name__ == "__main__":
    main()
