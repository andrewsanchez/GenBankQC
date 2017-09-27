import click

import genbankfilter.filter as gbf


@click.command()
@click.option('-l', '--filter-level', type=float,
              help='Value to be used for all filters')
@click.option('-n', '--max_unknowns', type=int, default=200,
              help='Maximum number of acceptable unknown bases')
@click.option('-c', '--c-range', type=float, default=3.0,
              help='Tolerance level for number of contigs',)
@click.option('-s', '--s-range', type=float, default=3.0,
              help='Tolerance level for the assembly size')
@click.option('-m', '--m-range', type=float, default=3.0,
              help='Tolerance level for MASH distances')
@click.option('-f', '--filter-only',
              help="Run filtering without running MASH or stats if current "
              "dmx.csv and stats.csv already exist.",
              default=False, is_flag=True)
@click.option('-d', '--dry-run', is_flag=True)
@click.argument('species-dir', type=click.Path(exists=True, file_okay=False))
def cli(filter_level, max_unknowns, c_range, s_range, m_range,
        path, dry_run, filter_only):
    """ Assess the integrity of your FASTA collection."""
    species = gbf.FilteredSpecies(path, max_unknowns,
                                  c_range, s_range, m_range)
    if dry_run:
        click.echo(print(species))
    # elif not gbf.min_fastas_check(path):
    #     click.echo("{} contains less than 5 genomes.".format(path))
    #     pass
    elif filter_only:
        gbf.filter_all(species)
    # else:
    #     dmx = mash.mash(path)
    #     gbf.stats_and_filter(path, dmx, filter_ranges)
