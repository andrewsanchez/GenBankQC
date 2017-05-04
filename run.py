#!/usr/bin/env python

import os
import argparse
import logging

from genbankfilter import get_resources
from genbankfilter import config
from genbankfilter import curate
from genbankfilter import mash

def setup(genbank_mirror, species, update):

    path_vars = config.instantiate_path_vars(genbank_mirror)
    info_dir, slurm, out, logger = path_vars
    assembly_summary = get_resources.get_resources(genbank_mirror, logger, update)
    species = curate.get_species_list(assembly_summary, species)

    return path_vars, assembly_summary, species

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-s", "--species", help = 'List of species', nargs='+', default='all')
    parser.add_argument("-p", "--parallel", help = 'Submit jobs in parallel via SLURM pipeline tool.', action="store_true")
    parser.add_argument("--no_update", action='store_true')
    args = parser.parse_args()

    if args.no_update:
        update=False
    else:
        update=True

    genbank_mirror = args.genbank_mirror
    path_vars, assembly_summary, species = setup(genbank_mirror, args.species, update)
    info_dir, slurm, out, logger = path_vars
    genbank_status = curate.assess_genbank_mirror(genbank_mirror, assembly_summary, species, logger)
    local_genomes, sketch_files, missing_sketch_files = genbank_status
    mash.sketch(genbank_mirror, assembly_summary, missing_sketch_files, logger)
    # TODO: Generate a new list of species from accession ids in missing_sketch_files so that funs below
    # are only run on species that get new sketch files
    mash.paste(genbank_mirror, assembly_summary, species, logger)
    mash.dist(genbank_mirror, assembly_summary, species, logger)

if __name__ == "__main__":
    main()
