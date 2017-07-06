import os
import argparse
import logging
import pandas as pd

from genbankfilter import get_resources
from genbankfilter import config
from genbankfilter import curate
from genbankfilter import mash
from genbankfilter import filter

def mash_stats_and_filter():
    clean_up(species_dir)
    assess_fastas(species_dir)
    dst_mx = mash(species_dir, mash_exe)
    stats = generate_species_stats(fasta_dir, dst_mx)
    filter_med_ad(species_dir, stats, max_ns, c_range, s_range, m_range)

def just_filter():
    stats = os.path.join(os.path.join(species_dir, "info"), "stats.csv")
    stats = pd.read_csv(stats, index_col=0)
    filter_med_ad(species_dir, stats, max_ns, c_range, s_range, m_range)

def stats_are_current():
    genomes_in_dir = 0
    for f in os.listdir(species_dir):
        if f.endswith("fasta"):
            genomes_in_dir += 1
    genomes_in_stats = len(stats.index)
    if genomes_in_dir == genomes_in_stats:
        print("stats.csv is current. Filtering will be based on existing stats.csv")
        return True
    else:
        print("stats.csv is not current.  MASH will be run and stats.csv will be updated.")
        return False

def main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("species_dir", help = "directory containing your FASTA files")
    parser.add_argument("-d", "--directories", help = "The complete path to one or more directories containing\
            FASTA's you want to run the filters on.", action="store_true")
    parser.add_argument("--from_list", help = "Specify a list of one more more directories to fun filters on.", nargs="+")
    parser.add_argument("-p", "--parent_dir", help = "The parent directory containing subdirectories with FASTA's for\
            each of the collections you want to run the filters on.")
    parser.add_argument("-x", "--mash_exe", help = "Path to MASH")
    parser.add_argument("-n", "--max_n_count", help = "Maximum number of N's acceptable", type=int, default=100)
    parser.add_argument("-c", "--c_range", help = "", type=float)
    parser.add_argument("-s", "--s_range", help = "", type=float)
    parser.add_argument("-m", "--m_range", help = "", type=float)
    parser.add_argument("-l", "--filter_level", help = "Value to be used for all filters", type=float)

    args = parser.parse_args()
    species_dir = args.fasta_dir
    # max_ns = args.max_n_count 

    # if args.filter_level:
    #     c_range = args.filter_level
    #     s_range = args.filter_level
    #     m_range = args.filter_level
    # if args.c_range:
    #     c_range = args.c_range
    # if args.s_range:
    #     s_range = args.s_range
    # if args.m_range:
    #     m_range = args.m_range

    if len(os.listdir(species_dir)) <= 5: # pass if there are <= 5 FASTA's
        print("{} contains less than 5 genomes.".format(species_dir))
    # elif os.path.isfile(os.path.join(species_dir, "info", "stats.csv")):
    #     stats = os.path.join(species_dir, "info", "stats.csv")
    #     stats = pd.read_csv(stats, index_col=0)
    #     if stats_are_current():
    #         filter_med_ad(species_dir, stats)
    #     else:
    #         mash_stats_and_filter()
    else:
        mash_stats_and_filter()

if __name__ == '__main__':
    main()
