#!/usr/bin/env python

import os, argparse
import pandas as pd
from shutil import rmtree

def make_passed_dir(passed_dir, passed_genomes):
    
    if os.path.isdir(passed_dir):
        rmtree(passed_dir)
    if not os.path.isdir(passed_dir):
        os.mkdir(passed_dir)

def link_passed_genomes(species_dir, passed_dir, passed_genomes):

    passed_genomes = pd.read_csv(passed_genomes, index_col=0)
    for src in passed_genomes.index:
        genome = "{}.fasta".format(src)
        src = os.path.join(species_dir, genome)
        dst = os.path.join(passed_dir, genome)
        os.link(src, dst)

def main():
    parser = argparse.ArgumentParser(description = "Choose which genomes to link into the passed directory.")
    parser.add_argument("passed_genomes",help = "The passed*.csv file containing the genomes you want to use.")
    args = parser.parse_args()

    passed_genomes = args.passed_genomes
    filter_level = passed_genomes.split("/")[-1].strip(".csv")
    species_dir = "/".join(passed_genomes.split("/")[:-1])
    passed_dir = os.path.join(species_dir, filter_level)
    make_passed_dir(passed_dir, passed_genomes)
    link_passed_genomes(species_dir, passed_dir, passed_genomes)

if __name__ == "__main__":
    main()
