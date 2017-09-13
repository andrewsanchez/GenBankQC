import click
import genbankfilter.filter as gbf
import genbankfilter.mash as mash


@click.command()
@click.option(
    '--mash-exe',
    help='Path to MASH executable if not in your PATH',
    default='~/usr/bin/mash')
@click.option(
    '-l',
    '--filter-level',
    help='Value to be used for all filters',
    type=float, )
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
    '-m',
    '--m-range',
    help='Filtering level for MASH distances',
    type=float,
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
    _FilteredSpecies = gbf.FilteredSpecies(species_dir)
    print(_FilteredSpecies)
    if not gbf.min_fastas_check(species_dir):
        click.echo("{} contains less than 5 genomes.".format(species_dir))
        pass
    if filter_only:
        gbf._filter_all(_FilteredSpecies)
    else:
        dmx = mash.mash(species_dir)
        gbf.stats_and_filter(species_dir, dmx, filter_ranges)
