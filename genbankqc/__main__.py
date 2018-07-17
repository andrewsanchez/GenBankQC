import click
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
@click.argument('path', type=click.Path(exists=True, file_okay=False))
@click.option('--unknowns', '-n',
              type=int, default=200,
              help='Maximum number of unknown bases (not A, T, C, G)')
@click.option('--contigs', '-c',
              type=float, default=3.0,
              help='Acceptable deviations from median number of contigs')
@click.option('--assembly_size', '-s',
              type=float, default=3.0,
              help='Acceptable deviations from median assembly size')
@click.option('--distance', '-d',
              type=float, default=3.0,
              help='Acceptable deviations from median MASH distances')
@click.option('--all',
              type=float,
              help='Acceptable deviations for all metrics')
@click.option('--metadata', is_flag=True,
              help='Get metadata for genome at PATH',)
def species(ctx, path, unknowns, contigs, assembly_size, distance, all,
            metadata):
    """
    Run commands on a single species.
    """
    kwargs = {"max_unknowns": unknowns,
              "contigs": contigs,
              "assembly_size": assembly_size,
              "mash": distance,
              "assembly_summary": ctx.assembly_summary}
    species = Species(path, **kwargs)
    species.qc()
    if metadata:
        species.metadata()


@cli.command()
@click.pass_obj
@click.argument('path', type=click.Path(exists=True, dir_okay=False))
@click.option('--metadata', is_flag=True,
              help='Get metadata for genome at PATH')
def genome(ctx, path, metadata):
    """
    Get information about a single genome.
    """
    genome = Genome(path, ctx.assembly_summary)
    if metadata:
        click.echo(genome.metadata)
