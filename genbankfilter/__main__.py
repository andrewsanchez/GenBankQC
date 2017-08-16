import os
import argparse
import logging
import pandas as pd
import click
import genbankfilter.filter as gbf
import genbankfilter.mash as mash


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
        print(
            "stats.csv is current. Filtering will be based on existing stats.csv"
        )
        return True
    else:
        print(
            "stats.csv is not current.  MASH will be run and stats.csv will be updated."
        )
        return False


@click.command()
@click.option(
    '--mash-exe',
    help='Path to MASH executable if not in your PATH',
    default='~/usr/bin/mash')
@click.option(
    '-l',
    '--filter-level',
    help='Value to be used for all filters',
    type=float,)
@click.option(
    '-n',
    '--max_n_count',
    help='Maximum number of acceptable unknown bases',
    type=int,
    default=100)
@click.option(
    '-c',
    '--c-range',
    help='Filtering level for number of contigs',
    type=float,
    default=3.0)
@click.option(
    '-s',
    '--s-range',
    help='Filtering level for the assembly size',
    type=float,
    default=3.0)
@click.option(
    '-m', '--m-range', help='Filtering level for MASH distances', type=float,
    default=3.0)
@click.option(
    '-f',
    '--filter-only',
    help="Run filtering without running MASH or stats.  Use this option "
    "if your species directory already contains the distance matrix, "
    "stats.csv, and you want to run the filters with different parameters",
    default=False,
    is_flag=True)
@click.argument('species-dir', type=click.Path(exists=True, file_okay=False))
def cli(mash_exe, filter_level, max_n_count, c_range, s_range, m_range,
        species_dir, filter_only):
    """
    Assess the integrity of your FASTA collection.
    """

    if filter_level:
        filter_ranges = max_n_count, filter_level, filter_level, filter_level
    else:
        filter_ranges = max_n_count, c_range, s_range, m_range
    click.echo('Filtering levels:')
    click.echo('Max Unknowns:  {}'.format(max_n_count))
    click.echo('Contigs:  {}'.format(c_range))
    click.echo('Assembly size:  {}'.format(s_range))
    click.echo('MASH distances:  {}'.format(m_range))

    if not gbf.min_fastas_check(species_dir):
        click.echo("{} contains less than 5 genomes.".format(species_dir))
        pass
    # elif os.path.isfile(os.path.join(species_dir, "info", "stats.csv")):
    #     stats = os.path.join(species_dir, "info", "stats.csv")
    #     stats = pd.read_csv(stats, index_col=0)
    #     if stats_are_current():
    #         filter_med_ad(species_dir, stats)
    #     else:
    #         mash_stats_and_filter()

    if filter_only:
        gbf.filter_only(species_dir, dmx, filter_ranges)
    else:
        dmx = mash.mash(species_dir)
        gbf.stats_and_filter(species_dir, dmx, filter_ranges)
