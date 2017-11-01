import click

from genbank_qc import Species, Genbank


@click.command()
@click.option('-l', '--filter-level', type=float,
              help='Value to be used for all filters')
@click.option('-n', '--max_unknowns', type=int, default=200,
              help='Maximum number of acceptable unknown bases')
@click.option('-c', '--c-tolerance', type=float, default=3.0,
              help='Tolerance level for number of contigs',)
@click.option('-s', '--s-tolerance', type=float, default=3.0,
              help='Tolerance level for the assembly size')
@click.option('-m', '--m-tolerance', type=float, default=3.0,
              help='Tolerance level for MASH distances')
@click.option('-d', '--dry-run', is_flag=True)
@click.option('--subdir', is_flag=True)
@click.argument('path', type=click.Path(exists=True, file_okay=False))
def cli(filter_level, max_unknowns, c_tolerance, s_tolerance, m_tolerance,
        dry_run, subdir, path):
    """Assess the integrity of your FASTA collection."""
    if subdir:
        species = Species(path, max_unknowns, c_tolerance, s_tolerance,
                          m_tolerance)
        species.qc()
    else:
        Genbank(path).qc()
