import os
import pandas as pd
from collections import namedtuple
from Bio import SeqIO
import glob
import re


def generate_stats(species_dir, dst_mx):

    fastas = (f for f in os.listdir(species_dir) if f.endswith('fasta'))
    file_names, contig_totals, assembly_sizes, n_counts = [], [], [], []

    for f in fastas:
        fasta = (os.path.join(species_dir, f))
        name = re.search('(GCA.*)(.fasta)', f).group(1)

        # Get all contigs for current FASTA
        try:
            contigs = [seq.seq for seq in SeqIO.parse(fasta, "fasta")]
        except UnicodeDecodeError:
            print("{} threw UnicodeDecodeError".format(f))
        # Length of each contig
        assembly_size = sum([len(str(seq)) for seq in contigs])
        # N_Count for each contig
        N_Count = [len(re.findall("[^ATCG]", str(seq))) for seq in contigs]

        file_names.append(name)
        assembly_sizes.append(assembly_size)
        contig_totals.append(len(contigs))
        n_counts.append(sum(N_Count))

    SeqDataSet = list(zip(n_counts, contig_totals, assembly_sizes, dst_mx.mean()))
    stats = pd.DataFrame(
        data=SeqDataSet,
        index=file_names,
        columns=["N_Count", "Contigs", "Assembly_Size", "MASH"],
        dtype="float64")

    return stats


def filter_med_ad(species_dir, stats, filter_ranges):

    max_n_count, c_range, s_range, m_range = filter_ranges
    filter_ranges = "{}_{}_{}".format(c_range, s_range, m_range)

    filter_summary = pd.DataFrame()

    # Filter based on N's first
    passed_I = stats[stats["N_Count"] <= max_n_count]
    failed = pd.DataFrame(index=stats.index, columns=stats.columns)
    failed_N = []
    failed_N_count_ixs = [i for i in stats.index if i not in passed_I.index]
    for i in stats.index:
        if i not in passed_I.index:
            failed_N.append(i)
            failed["N_Count"][i] = stats["N_Count"][i]
        else:
            failed["N_Count"][i] = "passed"
    filter_summary.set_value(filter_ranges, "N's", len(failed_N))

    # Filter using special function for contigs
    if len(passed_I) > 5:
        filter_contigs_results = filter_contigs(
            stats, passed_I, filter_ranges, c_range, failed, filter_summary)
        passed_II = filter_contigs_results.passed
        failed_contigs = filter_contigs_results.failed

        if len(passed_II) > 5:
            assembly_med_ad = abs(passed_II["Assembly_Size"]
                                  - passed_II["Assembly_Size"].median()).mean(
                                  )  # Median absolute deviation
            assembly_dev_ref = assembly_med_ad * s_range
            passed_III = passed_II[abs(passed_II["Assembly_Size"] -
                                       passed_II["Assembly_Size"].median()) <=
                                   assembly_dev_ref]
            failed_assembly_size = []
            for i in passed_II.index:
                if i not in passed_III.index:
                    failed["Assembly_Size"][i] = stats["Assembly_Size"][i]
                    failed_assembly_size.append(i)
                else:
                    failed["Assembly_Size"][i] = "passed"

            assembly_lower = passed_II["Assembly_Size"].median(
            ) - assembly_dev_ref
            assembly_upper = passed_II["Assembly_Size"].median(
            ) + assembly_dev_ref
            filter_summary.set_value(filter_ranges, "Assembly_Size",
                                     len(failed_assembly_size))
            filter_summary.set_value(filter_ranges, "Assembly_Range",
                                     "{:.0f}-{:.0f}".format(
                                         assembly_lower, assembly_upper))

            if len(passed_III) > 5:
                mash_med_ad = abs(passed_III["MASH"] - passed_III["MASH"].
                                  median()).mean()  # Median absolute deviation
                mash_dev_ref = mash_med_ad * m_range
                passed_final = passed_III[
                    abs(passed_III["MASH"] - passed_III["MASH"].median()) <=
                    mash_dev_ref]
                failed_mash = []
                for i in passed_III.index:
                    if i not in passed_final.index:
                        failed["MASH"][i] = stats["MASH"][i]
                        failed_mash.append(i)

                        failed["MASH"][i] = "passed"

                mash_lower = passed_II["MASH"].median() - mash_dev_ref
                mash_upper = passed_II["MASH"].median() + mash_dev_ref
                filter_summary.set_value(filter_ranges, "MASH",
                                         len(failed_mash))
                filter_summary.set_value(filter_ranges, "MASH_Range",
                                         "{:04.3f}-{:04.3f}".format(
                                             mash_lower, mash_upper))

            else:
                passed_final = passed_III
                print(
                    "Removing genomes outside the range of {}-{} for assembly size resulted in < 5 genomes.\n\
                        Filtering will not commence past this stage."
                    .format(assembly_lower, assembly_upper))
        else:
            passed_final = passed_II
            print(
                "Removing genomes outside the range of acceptable number of contigs resulted in < 5 genomes.\n\
                    Filtering will not commence past this stage.")
    else:
        passed_final = passed_I
        print(
            "Removing genomes with > than {} N's resulted in dataset with < 5 genomes.\n\
                Filtering will not commence past this stage."
            .format(max_n_count))

    failed.drop(list(passed_final.index), inplace=True)
    filter_summary.set_value(filter_ranges, "Filtered", "{}/{}".format(
        len(failed), len(stats)))
    return filter_summary, failed, passed_final


def stats_and_filter(species_dir, dst_mx, filter_ranges):
    stats = generate_stats(species_dir, dst_mx)
    results = filter_med_ad(species_dir, stats, filter_ranges)
    filter_summary, failed, passed_final = results
    stats.to_csv(os.path.join(species_dir, 'stats.csv'))
    filter_summary.to_csv(os.path.join(species_dir, 'summary.csv'))
    failed.to_csv(os.path.join(species_dir, 'failed.csv'))
    passed_final.to_csv(os.path.join(species_dir, 'passed.csv'))


def filter_contigs(stats, passed_I, filter_ranges, c_range, failed,
                   summary_df):

    contigs = passed_I["Contigs"]
    contigs_above_median = contigs[contigs >= contigs.median()]
    contigs_below_median = contigs[contigs <= contigs.median()]
    contigs_lower = contigs[
        contigs <=
        10]  # Save genomes with < 10 contigs to add them back in later.
    contigs = contigs[
        contigs >
        10]  # Only look at genomes with > 10 contigs to avoid throwing off the Median AD
    contigs_med_ad = abs(contigs -
                         contigs.median()).mean()  # Median absolute deviation
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
    summary_df.set_value(filter_ranges, "Contigs_Range",
                         "{:.0f}-{:.0f}".format(contigs_lower, contigs_upper))

    results = namedtuple("filter_contigs_results", ["passed", "failed"])
    filter_contigs_results = results(passed_II, failed_contigs)

    return filter_contigs_results


def stats_and_filter(species_dir, dst_mx, filter_ranges):
    stats = generate_stats(species_dir, dst_mx)
    stats.to_csv(os.path.join(species_dir, 'stats.csv'))
    results = filter_med_ad(species_dir, stats, filter_ranges)
    filter_summary, failed, passed_final = results
    filter_summary.to_csv(os.path.join(species_dir, 'summary.csv'), index_label='Filter Ranges')
    failed.to_csv(os.path.join(species_dir, 'failed.csv'))
    passed_final.to_csv(os.path.join(species_dir, 'passed.csv'))


def assess_fastas(fasta_dir):

    # Check for empty FASTA's and move them before running MASH
    info = os.path.join(fasta_dir, "info")
    empty = os.path.join(info, "corrupt_fastas")
    for f in os.listdir(fasta_dir):
        if f.endswith("fasta") and os.path.getsize(
                os.path.join(fasta_dir, f)) == 0:
            print("{} is empty and will be moved to {} before runing MASH.".
                  format(f, empty))
            if not os.path.isdir(empty):
                os.mkdir(empty)
            src = os.path.join(fasta_dir, f)
            dst = os.path.join(empty, f)
            shutil.move(src, dst)


def check_passed_dir(species_dir):
    passed_dir = os.path.join(species_dir, "passed")
    if os.path.isdir(passed_dir):
        shutil.rmtree(passed_dir)
    os.mkdir(passed_dir)
    return passed_dir


def min_fastas_check(species_dir):
    """
    Check if speices_dir contains at least 5 FASTAs
    """
    if len(os.listdir(species_dir)) <= 5:
        return False
    return True


def link_passed_genomes(species_dir, passed_final, passed_dir):
    for genome in passed_final.index:
        fasta = "{}.fasta".format(genome)
        # glob_pattern = "{}*fasta".format(genome)
        # fasta = glob.glob(os.path.join(species_dir, glob_pattern))[0]
        # fasta = fasta.split('/')[-1]
        src = os.path.join(species_dir, fasta)
        dst = os.path.join(passed_dir, fasta)
        os.link(src, dst)


# Move to config
def clean_up(species_dir):

    sketch_file = os.path.join(species_dir, "all.msh")
    dst_mx = os.path.join(species_dir, "dst_mx.csv")
    filter_log = os.path.join(species_dir, "filter_log.txt")

    info = os.path.join(species_dir, "info")
    if not os.path.isdir(info):
        os.mkdir(info)

    files = [sketch_file, dst_mx, filter_log]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)


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
                dst_mx = pd.read_csv(all_dist, index_col=0, delimiter="\t")
                clean_up_matrix(info_dir, dst_mx)  # cleans up matrix in place
                print("Formatted matrix for {}".format(d))
                print("{} out of {}".format(x, total_species))
                x += 1
            except FileNotFoundError:
                continue
            except pd.io.parsers.EmptyDataError:
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
        dst_mx_all = os.path.join(info_dir, "dst_mx_all.csv")
        stats = os.path.join(info_dir, "stats.csv")
        if os.path.isfile(dst_mx_all) and not os.path.isfile(stats):
            generate_fasta_stats(fasta_dir,
                                 pd.read_csv(
                                     dst_mx_all, index_col=0, delimiter="\t"))
            print("Generated stats for {} out of {}".format(
                x, len(all_genbank_species)))
            x += 1
