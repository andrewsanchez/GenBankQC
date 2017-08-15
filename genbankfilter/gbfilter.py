import os
import re
import shutil
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from collections import namedtuple


def generate_stats(species_dir, dmx):
    """
    Generate a data frame containing all of the stats for genomes
    in species_dir.
    """

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

    SeqDataSet = list(
        zip(n_counts, contig_totals, assembly_sizes, dmx.mean()))
    stats = pd.DataFrame(
        data=SeqDataSet,
        index=file_names,
        columns=["N_Count", "Contigs", "Assembly_Size", "MASH"],
        dtype="float64")

    return stats


def filter_all(species_dir, stats, tree, filter_ranges):
    """
    This function strings together all of the steps
    involved in filtering your genomes.
    """

    max_n_count, c_range, s_range, m_range = filter_ranges
    summary = {}
    failed = pd.DataFrame(index=stats.index, columns=stats.columns)

    # Filter based on N's first
    passed_N_count, failed_N_count, failed = filter_Ns(stats, summary, failed,
                                                       tree, max_n_count)
    # Filter contigs
    if check_df_len(passed_N_count):
        filter_results = filter_contigs(stats, passed_N_count, c_range, failed,
                                        summary, tree)
        passed = filter_results.passed
    else:
        print("Filtering based on unknown bases resulted in < 5 genomes.  "
              "Filtering will not commence past this stage.")
    for criteria in ["Assembly_Size", "MASH"]:
        if check_df_len(passed):
            filter_results = filter_med_ad(criteria, passed, failed, summary,
                                           tree, s_range)
            passed = filter_results.passed
        else:
            print("Filtering based on {} resulted in < 5 genomes.  "
                  "Filtering will not commence past this stage.".format(
                      criteria))
            break

    failed.drop(list(passed.index), inplace=True)
    write_summary(species_dir, summary, filter_ranges)
    tree.render(os.path.join(species_dir, 'tree.png'))
    tree.render(os.path.join(species_dir, 'tree.svg'))
    tree.render(os.path.join(species_dir, 'tree.pdf'))
    tree.show()
    return failed, passed


def filter_Ns(stats, summary, failed, tree, max_n_count):
    """
    Identify genomes with too many unknown bases.
    """
    passed_N_count = stats[stats["N_Count"] <= max_n_count]
    failed_N_count = stats[stats["N_Count"] >= max_n_count]
    failed["N_Count"][passed_N_count.index] = "+"
    for i in failed_N_count.index:
        failed["N_Count"][i] = stats["N_Count"][i]
    summary["N_Count"] = (max_n_count, len(failed_N_count))
    color_clade(tree, 'N_Count', failed_N_count.index)
    return passed_N_count, failed_N_count, failed


def filter_contigs(stats, passed_N_count, c_range, failed, summary, tree):

    contigs = passed_N_count["Contigs"]
    # Only look at genomes with > 10 contigs to avoid throwing off the
    # Median AD Save genomes with < 10 contigs to add them back in later.
    not_enough_contigs = contigs[contigs <= 10]
    contigs = contigs[contigs > 10]
    contigs_med_ad = abs(contigs -
                         contigs.median()).mean()  # Median absolute deviation
    contigs_dev_ref = contigs_med_ad * c_range
    contigs = contigs[abs(contigs - contigs.median()) <= contigs_dev_ref]
    # Add genomes with < 10 contigs back in
    contigs = pd.concat([contigs, not_enough_contigs])
    lower = contigs.median() - contigs_dev_ref
    upper = contigs.median() + contigs_dev_ref

    # Avoid returning empty DataFrame when no genomes are removed above
    if len(contigs) == len(passed_N_count):
        passed_contigs = passed_N_count
        failed_contigs = []
    else:
        failed_contigs = [
            i for i in passed_N_count.index if i not in contigs.index
        ]
        passed_contigs = passed_N_count.drop(failed_contigs)

        for i in failed_contigs:
            failed["Contigs"][i] = stats["Contigs"][i]
        for i in contigs.index:
            failed["Contigs"][i] = "+"
    color_clade(tree, 'Contigs', failed_contigs)
    range_str = "{:.0f}-{:.0f}".format(lower, upper)
    summary["Contigs"] = (range_str, len(failed_contigs))
    results = namedtuple("filter_contigs_results", ["passed", "failed"])
    filter_contigs_results = results(passed_contigs, failed_contigs)

    return filter_contigs_results


def filter_med_ad(criteria, passed, failed, summary, tree, f_range):
    """
    Filter based on median absolute deviation
    """
    # Get the median absolute deviation
    med_ad = abs(passed[criteria] - passed[criteria].median()).mean()
    dev_ref = med_ad * f_range
    passed = passed[abs(passed[criteria] - passed[criteria].median()) <=
                    dev_ref]
    failed = passed.index[abs(passed[criteria] - passed[criteria].median()) >=
                          dev_ref].tolist()
    lower = passed[criteria].median() - dev_ref
    upper = passed[criteria].median() + dev_ref
    range_str = '{:.0f}-{:.0f}'.format(lower, upper)
    summary[criteria] = (range_str, len(failed))
    results = namedtuple("filter_results", ["passed", "failed"])
    filter_results = results(passed, failed)
    # need a way to dynamically color based on criteria
    color_clade(tree, criteria, failed)
    return filter_results


def check_df_len(df, num=5):
    """
    Verify that df has > than num genomes
    """
    if len(df) > num:
        return True
    else:
        return False


def write_summary(species_dir, summary, filter_ranges):
    """
    Write a summary of the filtering results.
    """
    max_n_count, c_range, s_range, m_range = filter_ranges
    out = 'summary_{}-{}-{}-{}.txt'.format(max_n_count, c_range, s_range,
                                           m_range)
    out = os.path.join(species_dir, out)
    if os.path.isfile(out):
        os.remove(out)
    with open(out, 'a') as f:
        for k, v in summary.items():
            f.write('{}\n'.format(k))
            f.write('Range: {}\n'.format(v[0]))
            f.write('Filtered: {}\n\n'.format(v[1]))


def dmx_to_tree(dmx, species_dir):
    """
    Convert dmx to a nested representation of the
    lower diagonal of the distance matrix. Generate
    a _DistanceMatrix object and construct a tree.
    """
    from ete3 import Tree
    nw_file = os.path.join(species_dir, 'tree.nw')
    m = dmx.as_matrix()
    trild = np.tril(m, k=0)

    nested_dmx = []
    mx_len = len(trild)
    for i in np.arange(0, mx_len):
        tmp = trild[i, :i+1]
        nested_dmx.append(tmp.tolist())

    names = dmx.index.tolist()
    bio_dmx = _DistanceMatrix(names, nested_dmx)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(bio_dmx)
    Phylo.write(tree, nw_file, 'newick')
    tree = Tree(nw_file, 1)
    tree = style_tree(tree)
    return tree


def style_tree(tree):
    from ete3 import NodeStyle
    nstyle = NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 2
    nstyle["fgcolor"] = "black"
    for n in tree.traverse():
        n.set_style(nstyle)
    return tree


def color_clade(tree, criteria, to_color):
    """Color nodes using ete3
    """
    from ete3 import NodeStyle

    colors = {"N_Count": "red",
              "Contigs": "green",
              "MASH": "blue",
              "Assembly_Size": "yellow"}

    for genome in to_color:
        genome = re.match('.*(GCA_\d+\.\d)', genome).group(1)
        n = tree.get_leaves_by_name(genome).pop()
        nstyle = NodeStyle()
        nstyle["fgcolor"] = colors[criteria]
        nstyle["size"] = 5
        n.set_style(nstyle)


def stats_and_filter(species_dir, dmx, filter_ranges):
    stats = generate_stats(species_dir, dmx)
    stats.to_csv(os.path.join(species_dir, 'stats.csv'))
    tree = dmx_to_tree(dmx, species_dir)
    results = filter_all(species_dir, stats, tree, filter_ranges)
    failed, passed_final = results
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
    dmx = os.path.join(species_dir, "dmx.csv")
    filter_log = os.path.join(species_dir, "filter_log.txt")

    info = os.path.join(species_dir, "info")
    if not os.path.isdir(info):
        os.mkdir(info)

    files = [sketch_file, dmx, filter_log]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)
