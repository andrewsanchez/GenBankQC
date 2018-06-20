import click
import traceback
from collections import namedtuple
from genbankqc import Genbank
from genbankqc import Species
from genbankqc import Metadata



help_text = """
Assess the integrity of your FASTA collection.

Run genbankqc on subdirectories in parent directory PATH.

Subdirectories should contain 5 or more FASTAs.

Specify a single species directory with the --species flag.
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
@click.option('--species', is_flag=True,
              help='Run on single species')
@click.argument('path', type=click.Path(exists=True, file_okay=False))
def cli(filter_level,
        max_unknowns,
        c_deviations,
        s_deviations,
        m_deviations,
        dry_run,
        species,
        path):
    if species:
        from genbankqc import Species
        try:
            s = Species(path, max_unknowns, c_deviations, s_deviations,
                        m_deviations)
            s.qc()
            print("Completed", s.species)
            print(s)
        except Exception:
            print('Failed', s.species)
            traceback.print_exc()
    else:
        from genbankqc import Genbank
        genbank = Genbank(path)
        genbank.qc()
