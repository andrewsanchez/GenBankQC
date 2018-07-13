import click
import traceback
from collections import namedtuple
from genbankqc import Genbank
from genbankqc import Genome
from genbankqc import Species


class CLIGroup(click.Group):
    def parse_args(self, ctx, args):
        try:
            if args[0] in self.commands:
                if len(args) == 1 or args[1] not in self.commands:
                    args.insert(0, '')
        except IndexError:
            pass
        super(CLIGroup, self).parse_args(ctx, args)


@click.group(invoke_without_command=True, no_args_is_help=True, cls=CLIGroup)
@click.pass_context
@click.argument('path', type=click.Path(), required=False)
def cli(ctx, path):
    """
    Assess the integrity of your genomes through automated analysis of
    species-based statistics and metadata.
    """
    # TODO: Option for basic info about PATH
    _ctx = namedtuple('ctx', ['genbank', 'assembly_summary'])
    genbank = Genbank(path)
    ctx.obj = _ctx(genbank=genbank, assembly_summary=genbank.assembly_summary)
    if ctx.invoked_subcommand is None:
        genbank.qc()


@cli.command()
@click.pass_obj
@click.option('--max_unknowns', '-n', type=int,
              default=200, help='Maximum number of unknown bases')
@click.option('--c-deviations', '-c', type=float,
              default=3.0, help='Deviations for number of contigs',)
@click.option('--s-deviations', '-s', type=float,
              default=3.0, help='Deviations for the assembly size')
@click.option('--m-deviations', '-m', type=float,
              default=3.0, help='Deviations for MASH distances')
@click.option('--filter-level', '-l', type=float,
              help='Deviations for all metrics')
@click.argument('path', type=click.Path(exists=True, file_okay=False))
def species(ctx, max_unknowns, c_deviations,
            s_deviations, m_deviations, filter_level, path):
    """
    Run qc command on given species
    """
    try:
        species = Species(path, max_unknowns, c_deviations, s_deviations,
                          m_deviations, ctx.assembly_summary)
        species.qc()
    except Exception:
        click.echo('Failed', species.species)
        traceback.print_exc()


@cli.command()
@click.pass_obj
@click.argument('path', type=click.Path(exists=True, dir_okay=False))
@click.option('--metadata', help='Get metadata for genome at PATH',
              is_flag=True)
def genome(ctx, path, metadata):
    """
    Get information about a genome or list of genomes.
    """
    genome = Genome(path, ctx.assembly_summary)
    click.echo(genome.metadata)


@cli.command()
def metadata():
    """
    Generate Metadata
    """
    pass
