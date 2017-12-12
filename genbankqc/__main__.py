import click

from genbankqc import Species, Genbank

help_text = """
Assess the integrity of your FASTA collection.

Run genbankqc on subdirectories in parent directory GENKBANK.
Subdirectories should contain 5 or more FASTAs.
"""


@click.command(help=help_text)
@click.option('-n', '--max_unknowns', type=int, default=200,
              help='Maximum number of unknown bases')
@click.option('-c', '--c-deviations', type=float, default=3.0,
              help='Deviations for number of contigs',)
@click.option('-s', '--s-deviations', type=float, default=3.0,
              help='Deviations for the assembly size')
@click.option('-m', '--m-deviations', type=float, default=3.0,
              help='Deviations for MASH distances')
@click.option('-l', '--filter-level', type=float,
              help='Deviations for all metrics')
@click.option('-d', '--dry-run', is_flag=True)
@click.option('--subdir', is_flag=True,
              help='Run on single species')
@click.argument('genbank', type=click.Path(exists=True, file_okay=False))
def cli(filter_level, max_unknowns, c_deviations, s_deviations, m_deviations,
        dry_run, subdir, genbank):
    if subdir:
        species = Species(genbank, max_unknowns, c_deviations, s_deviations,
                          m_deviations)
        species.qc()
    else:
        Genbank(genbank).qc()
