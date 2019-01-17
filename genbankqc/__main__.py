import os
import re
import click
import logbook

from genbankqc import Genbank
from genbankqc import Genome
from genbankqc import Species


class CLIGroup(click.Group):
    def parse_args(self, ctx, args):
        try:
            if args[0] in self.commands:
                if len(args) == 1 or args[1] not in self.commands:
                    args.insert(0, "")
        except IndexError:
            pass
        super(CLIGroup, self).parse_args(ctx, args)


@click.group(invoke_without_command=True, no_args_is_help=True, cls=CLIGroup)
@click.pass_context
@click.argument("path", type=click.Path(), required=False)
def cli(ctx, path):
    """Assess the integrity of your genomes through automated analysis of
    species-based statistics and metadata.
    """
    if ctx.invoked_subcommand is None:
        logbook.set_datetime_format("local")
        handler = logbook.TimedRotatingFileHandler(
            os.path.join(path, ".logs", "qc.log"), backup_count=10
        )
        handler.push_application()
        genbank = Genbank(path)
        genbank.qc()


@cli.command()
@click.argument("path", type=click.Path())
@click.argument("email")
@click.option("--update/--no-update", " /-U", default=True, help="Update metadata")
def metadata(path, email, update):
    """Download assembly_summary.txt and BioSample metadata."""
    logbook.set_datetime_format("local")
    handler = logbook.TimedRotatingFileHandler(
        os.path.join(path, ".logs", "metadata.log"), backup_count=10
    )
    handler.push_application()
    genbank = Genbank(path)
    metadata = genbank.metadata(email=email, update=update)
    genbank.species_metadata(metadata)


@cli.command()
@click.argument("path", type=click.Path(exists=True, file_okay=False))
@click.option(
    "--unknowns",
    "-n",
    type=int,
    default=200,
    help="Maximum number of unknown bases (not A, T, C, G)",
)
@click.option(
    "--contigs",
    "-c",
    type=float,
    default=3.0,
    help="Acceptable deviations from median number of contigs",
)
@click.option(
    "--assembly_size",
    "-s",
    type=float,
    default=3.0,
    help="Acceptable deviations from median assembly size",
)
@click.option(
    "--distance",
    "-d",
    type=float,
    default=3.0,
    help="Acceptable deviations from median MASH distances",
)
@click.option("--all", type=float, help="Acceptable deviations for all metrics")
@click.option("--metadata", is_flag=True, help="Get metadata for genome at PATH")
def species(path, unknowns, contigs, assembly_size, distance, all, metadata):
    """Run commands on a single species"""
    kwargs = {
        "max_unknowns": unknowns,
        "contigs": contigs,
        "assembly_size": assembly_size,
        "mash": distance,
    }
    logbook.set_datetime_format("local")
    handler = logbook.TimedRotatingFileHandler(
        os.path.join(path, ".logs", "qc.log"), backup_count=10
    )
    handler.push_application()
    species = Species(path, **kwargs)
    species.qc()
    if metadata:
        species.metadata()


@cli.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
@click.option("--metadata", is_flag=True, help="Get metadata for genome at PATH")
def genome(path, metadata):
    """ Get information about a single genome."""

    genome = Genome(path)
    if metadata:
        click.echo(genome.metadata)


@cli.command()
@click.pass_obj
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
def log_stats(ctx, path):
    """
    Summarize basic stats for a given log file.
    """
    log_file = os.path.join(ctx.genbank.path, ".logs", path)
    not_enough_genomes = [0, "Not enough genomes"]
    completed_metadata_command = [0, "Completed metadata command"]
    already_complete = [0, "Already complete"]
    tree_already_complete = [0, "Tree already complete"]
    generated_stats = [0, "Generated stats"]
    qc_completed = [0, "qc command completed"]
    stats = [
        log_file,
        not_enough_genomes,
        completed_metadata_command,
        already_complete,
        tree_already_complete,
        generated_stats,
        qc_completed,
    ]
    with open(log_file) as f:
        for line in f:
            if re.search(not_enough_genomes[1], line):
                not_enough_genomes[0] += 1
            elif re.search(completed_metadata_command[1], line):
                completed_metadata_command[0] += 1
            elif re.search(already_complete[1], line):
                already_complete[0] += 1
            elif re.search(tree_already_complete[1], line):
                tree_already_complete[0] += 1
            elif re.search(generated_stats[1], line):
                generated_stats[0] += 1
            elif re.search(qc_completed[1], line):
                qc_completed[0] += 1
    for i in stats:
        click.echo(i[1])
        click.echo(i[0])
