import os, argparse
import pandas as pd
from collections import namedtuple
from subprocess import Popen
from shutil import rmtree, move
from Bio import SeqIO
from re import findall
from pandas.io.parsers import EmptyDataError

# Move to config
def clean_up(fasta_dir):

    sketch_file = os.path.join(fasta_dir, "all.msh")
    distance_matrix = os.path.join(fasta_dir, "distance_matrix.csv")
    filter_log = os.path.join(fasta_dir, "filter_log.txt")

    info = os.path.join(fasta_dir, "info")
    if not os.path.isdir(info):
        os.mkdir(info)

    files = [sketch_file, distance_matrix, filter_log]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)

def mash(fasta_dir, mash_exe):
    info = os.path.join(fasta_dir, "info")
    sketch_file = os.path.join(info, "all.msh")
    all_fastas = os.path.join(fasta_dir, "*.fasta")
    distance_matrix = os.path.join(info, "distance_matrix.csv")
    sketch_command = "{} sketch -o {} {}".format(mash_exe, sketch_file, all_fastas)
    distance_command = "{} dist -t {} {} > {}".format(mash_exe, sketch_file, all_fastas, distance_matrix)
    #distance_command = "{} dist -p 4 -t {} {} > {}".format(mash_exe, sketch_file, all_fastas, distance_matrix)
    Popen(sketch_command, shell="True").wait()
    Popen(distance_command, shell="True").wait()
    distance_matrix = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    clean_up_matrix(info, distance_matrix)

    return distance_matrix

# Convenience
def pre_process_all(genbank_mirror):

    x = 1
    total_species = len(os.listdir(genbank_mirror))
    for d in os.listdir(genbank_mirror):
        fasta_dir = os.path.join(genbank_mirror, d)
        info_dir = os.path.join(fasta_dir, "info")
        all_dist = os.path.join(genbank_mirror, d, "all_dist.msh")
        if not os.path.isfile(all_dist):
            try:
                distance_matrix = pd.read_csv(all_dist, index_col=0, delimiter="\t")
                clean_up_matrix(info_dir, distance_matrix) # cleans up matrix in place
                print("Formatted matrix for {}".format(d))
                print("{} out of {}".format(x, total_species))
                x += 1
            except FileNotFoundError:
                continue
            except EmptyDataError:
                continue
        else:
            continue

# Convenience
def generate_stats_genbank(genbank_mirror):

    x = 1
    all_genbank_species = os.listdir(genbank_mirror)
    for d in all_genbank_species:
        print("generating stats for {}".format(d))
        fasta_dir = os.path.join(genbank_mirror, d)
        info_dir = os.path.join(fasta_dir, "info")
        distance_matrix_all = os.path.join(info_dir, "distance_matrix_all.csv")
        stats = os.path.join(info_dir, "stats.csv")
        if os.path.isfile(distance_matrix_all) and not os.path.isfile(stats):
            generate_fasta_stats(fasta_dir, pd.read_csv(distance_matrix_all, index_col=0, delimiter="\t"))
            print("Generated stats for {} out of {}".format(x, len(all_genbank_species)))
            x += 1

def generate_stats(species_dir, dst_mx):

    for root, dirs, files, in os.walk(species_dir):
        file_names = []
        contig_totals = []
        assembly_sizes = []
        n_counts = []

        for f in files:
            if f.endswith(".fasta"):
                print("Getting stats for {}".format(f))
                fasta = (os.path.join(root, f))
                name = fasta.split("/")[-1].strip(".fasta")
                file_names.append(name)

                # Read all contigs for current fasta into list
                try:
                    contigs = [ seq.seq for seq in SeqIO.parse(fasta, "fasta") ]
                except UnicodeDecodeError:
                    print("{} threw UnicodeDecodeError".format(f))

                # Append the total number of contigs to contig_totals
                contig_totals.append(len(contigs))

                # Read the length of each contig into a list
                assembly_size = [ len(str(seq)) for seq in contigs ]
                # Append the sum of all contig lengths to lengths
                assembly_sizes.append(sum(assembly_size))

                # Read the N_Count for each contig into a list
                N_Count = [len(findall("[^ATCG]", str(seq))) for seq in contigs]
                # Append the total N_Count to n_counts
                n_counts.append(sum(N_Count))

        SeqDataSet = list(zip(n_counts, contig_totals, assembly_sizes))
        stats = pd.DataFrame(data=SeqDataSet, index=file_names, columns=["N_Count", "Contigs", "Assembly_Size"], dtype="float64")
        mean_distances = dst_mx.mean()
        stats["MASH"] = mean_distances

        return stats

def filter_med_ad(species_dir, stats, max_n_count=100,
                  c_range=3.5, s_range=3.5, m_range=3.5):

    filter_ranges = "{}_{}_{}".format(c_range, s_range, m_range)

    filter_df = pd.DataFrame()

    # Filter based on N's first
    passed_I = stats[stats["N_Count"] <= max_n_count ]
    failed = pd.DataFrame(index=stats.index, columns=stats.columns)
    failed_N = []
    for i in stats.index:
        if i not in passed_I.index:
            failed_N.append(i)
            failed["N_Count"][i] = stats["N_Count"][i]
        else:
            failed["N_Count"][i] = "passed"
    filter_df.set_value(filter_ranges, "N's", len(failed_N))

    # Filter using special function for contigs
    if len(passed_I) > 5:
        filter_contigs_results = filter_contigs(stats, passed_I, filter_ranges, c_range, failed, filter_df)
        passed_II = filter_contigs_results.passed
        failed_contigs = filter_contigs_results.failed

        if len(passed_II) > 5:
            assembly_med_ad = abs(passed_II["Assembly_Size"] - passed_II["Assembly_Size"].median()).mean()# Median absolute deviation
            assembly_dev_ref = assembly_med_ad * s_range
            passed_III = passed_II[abs(passed_II["Assembly_Size"] - passed_II["Assembly_Size"].median()) <= assembly_dev_ref]
            failed_assembly_size = []
            for i in passed_II.index:
                if i not in passed_III.index:
                    failed["Assembly_Size"][i] = stats["Assembly_Size"][i]
                    failed_assembly_size.append(i)
                else:
                    failed["Assembly_Size"][i] = "passed"

            assembly_lower = passed_II["Assembly_Size"].median() - assembly_dev_ref
            assembly_upper = passed_II["Assembly_Size"].median() + assembly_dev_ref
            filter_df.set_value(filter_ranges, "Assembly_Size", len(failed_assembly_size))
            filter_df.set_value(filter_ranges, "Assembly_Range", "{:.0f}-{:.0f}".format(assembly_lower, assembly_upper))

            if len(passed_III) > 5:
                mash_med_ad = abs(passed_III["MASH"] - passed_III["MASH"].median()).mean()# Median absolute deviation
                mash_dev_ref = mash_med_ad * m_range
                passed_final = passed_III[abs(passed_III["MASH"] - passed_III["MASH"].median()) <= mash_dev_ref]
                failed_mash = []
                for i in passed_III.index:
                    if i not in passed_final.index:
                        failed["MASH"][i] = stats["MASH"][i]
                        failed_mash.append(i)

                        failed["MASH"][i] = "passed"

                mash_lower = passed_II["MASH"].median() - mash_dev_ref
                mash_upper = passed_II["MASH"].median() + mash_dev_ref
                filter_df.set_value(filter_ranges, "MASH", len(failed_mash))
                filter_df.set_value(filter_ranges, "MASH_Range", "{:04.3f}-{:04.3f}".format(mash_lower, mash_upper))

            else:
                passed_final = passed_III
                print("Removing genomes outside the range of {}-{} for assembly size resulted in < 5 genomes.\n\
                        Filtering will not commence past this stage.".format(assembly_lower, assembly_upper))
        else:
            passed_final = passed_II
            print("Removing genomes outside the range of acceptable number of contigs resulted in < 5 genomes.\n\
                    Filtering will not commence past this stage.")
    else:
        passed_final = passed_I
        print("Removing genomes with > than {} N's resulted in dataset with < 5 genomes.\n\
                Filtering will not commence past this stage.".format(max_n_count))

    failed.drop(list(passed_final.index), inplace=True)
    filter_df.set_value(filter_ranges, "Filtered", "{}/{}".format(len(failed), len(stats)))
    return filter_df, failed, passed_final

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
        distance_matrix = mash(fasta_dir, mash_exe)
        stats = generate_fasta_stats(fasta_dir, distance_matrix)
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
def filter_contigs(stats, passed_I, filter_ranges, c_range, failed, summary_df):

    contigs = passed_I["Contigs"]
    contigs_above_median = contigs[contigs >= contigs.median()]
    contigs_below_median = contigs[contigs <= contigs.median()]
    contigs_lower = contigs[contigs <= 10] # Save genomes with < 10 contigs to add them back in later.
    contigs = contigs[contigs > 10] # Only look at genomes with > 10 contigs to avoid throwing off the Median AD
    contigs_med_ad = abs(contigs - contigs.median()).mean() # Median absolute deviation
    contigs_dev_ref = contigs_med_ad * c_range
    contigs = contigs[abs(contigs - contigs.median()) <= contigs_dev_ref]
    contigs = pd.concat([contigs, contigs_lower])
    contigs_lower = contigs.median() - contigs_dev_ref
    contigs_upper = contigs.median() + contigs_dev_ref
  # contigs = pd.concat([contigs, contigs_lower, contigs_below_median])

    # Avoid returning empty DataFrame when no genomes are removed above
    if len(contigs) == len(passed_I):
        passed_II = passed_I
        failed_contigs = []
    else:
        failed_contigs = [i for i in passed_I.index if i not in contigs.index]
        passed_II = passed_I.drop(failed_contigs)

        for i in failed_contigs:
            failed["Contigs"][i] = stats["Contigs"][i]
        for i in contigs.index:
            failed["Contigs"][i] = "passed"

    summary_df.set_value(filter_ranges, "Contigs", len(failed_contigs))
    summary_df.set_value(filter_ranges, "Contigs_Range", "{:.0f}-{:.0f}".format(contigs_lower, contigs_upper))

    results = namedtuple("filter_contigs_results", ["passed", "failed"])
    filter_contigs_results = results(passed_II, failed_contigs)

    return filter_contigs_results

def assess_fastas(fasta_dir):

    # Check for empty FASTA's and move them before running MASH
    info = os.path.join(fasta_dir, "info")
    empty = os.path.join(info, "corrupt_fastas")
    for f in os.listdir(fasta_dir):
        if f.endswith("fasta") and os.path.getsize(os.path.join(fasta_dir, f)) == 0:
            print("{} is empty and will be moved to {} before runing MASH.".format(f, empty))
            if not os.path.isdir(empty):
                os.mkdir(empty)
            src = os.path.join(fasta_dir, f)
            dst = os.path.join(empty, f)
            move(src, dst)

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
        distance_matrix = mash(fasta_dir, mash_exe)
        stats = generate_fasta_stats(fasta_dir, distance_matrix)
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

def make_passed_dir(passed_dir):
    
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

if __name__ == '__main__':
    main()
