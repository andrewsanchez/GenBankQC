import os
import argparse
import logging

from genbankfilter import get_resources
from genbankfilter import config
from genbankfilter import curate
from genbankfilter import mash
from genbankfilter import filter

def main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("fasta_dir", help = "directory containing your FASTA files")
    parser.add_argument("-d", "--directories", help = "The complete path to one or more directories containing\
            FASTA's you want to run the filters on.", action="store_true")
    parser.add_argument("-x", "--mash_exe", help = "Path to MASH")
    parser.add_argument("-p", "--parent_dir", help = "The parent directory containing subdirectories with FASTA's for\
            each of the collections you want to run the filters on.")
    parser.add_argument("--from_list", help = "Specify a list of one more more directories to fun filters on.", nargs="+")
    parser.add_argument("-n", "--max_n_count", help = "Maximum number of N's acceptable", type=int, default=100)
    parser.add_argument("-c", "--c_range", help = "", type=float, default=3.5)
    parser.add_argument("-s", "--s_range", help = "", type=float, default=3.5)
    parser.add_argument("-m", "--m_range", help = "", type=float, default=3.5)
    parser.add_argument("-l", "--filter_level", help = "Value to be used for all filters", type=float)
    args = parser.parse_args()

    def mash_stats_and_filter():
        clean_up(fasta_dir)
        assess_fastas(fasta_dir)
        dst_mx = mash(fasta_dir, mash_exe)
        stats = generate_fasta_stats(fasta_dir, dst_mx)
        filter_med_ad(fasta_dir, stats, max_ns, c_range, s_range, m_range)

    def just_filter():
        stats = os.path.join(os.path.join(fasta_dir, "info"), "stats.csv")
        stats = pd.read_csv(stats, index_col=0)
        filter_med_ad(fasta_dir, stats, max_ns, c_range, s_range, m_range)

    def stats_are_current():
        genomes_in_dir = 0
        for f in os.listdir(fasta_dir):
            if f.endswith("fasta"):
                genomes_in_dir += 1
        genomes_in_stats = len(stats.index)
        if genomes_in_dir == genomes_in_stats:
            print("stats.csv is current. Filtering will be based on existing stats.csv")
            return True
        else:
            print("stats.csv is not current.  MASH will be run and stats.csv will be updated.")
            return False

    mash_exe = "/common/contrib/bin/mash-Linux64-v1.1.1/mash"
    fasta_dir = args.fasta_dir
    max_ns = args.max_n_count 

    if args.filter_level:
        c_range = args.filter_level
        s_range = args.filter_level
        m_range = args.filter_level
    else:
        c_range = args.c_range
        s_range = args.s_range
        m_range = args.m_range

    if len(os.listdir(fasta_dir)) <= 5: # pass if there are <= 5 FASTA's
        print("{} contains less than 5 genomes.".format(fasta_dir))
    elif os.path.isfile(os.path.join(fasta_dir, "info", "stats.csv")):
        stats = os.path.join(fasta_dir, "info", "stats.csv")
        stats = pd.read_csv(stats, index_col=0)
        if stats_are_current():
            filter_med_ad(fasta_dir, stats, max_ns, c_range, s_range, m_range)
        else:
            mash_stats_and_filter()
    else:
        mash_stats_and_filter()

if __name__ == '__main__':
    main()
