import argparse
import os
import glob
import gzip
import re
import logging

def clean_up(genbank_mirror, path_vars):

    info_dir, slurm, out, logger = path_vars

    latest_assembly_versions = os.path.join(info_dir, "latest_assembly_versions.csv")
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    slurm_script = os.path.join(slurm, "get_latest_assembly_versions.sbatch")
    sync_array = os.path.join(genbank_mirror, ".info", "slurm", "sync_array.txt")
    sync_array_script = os.path.join(slurm, 'sync_array_script.sbatch')
    grab_genomes_script = os.path.join(slurm, 'grab_genomes_script.sbatch')

    for f in [latest_assembly_versions, latest_assembly_versions_array, slurm_script, sync_array_script, grab_genomes_script, sync_array]:
        if os.path.isfile(f):
            os.remove(f)

def get_species_list(assembly_summary, species_list):

    if species_list == "all":

        species_list = assembly_summary.scientific_name[assembly_summary.scientific_name.notnull()]
        species_list = set(species_list.tolist())

    return species_list

def create_species_dirs(genbank_mirror, assembly_summary, logger, species_list):

    for species in species_list:
        try:
            species_dir = os.path.join(genbank_mirror, species)
        except TypeError:
            continue
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)
            logger.info("Directory created: {}".format(species))

def get_local_genomes(genbank_mirror):

    local_genomes = []

    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if re.match('GCA.*fasta', f):
                genome_id = '_'.join(f.split('_')[:2])
                local_genomes.append(genome_id)

    return local_genomes

def get_new_genome_list(genbank_mirror, assembly_summary, local_genomes, species_list):

    # TODO: Faster way to do this in pandas?

    new_genomes = []

    for species in species_list:
        latest_assembly_versions = assembly_summary.index[assembly_summary.scientific_name == species].tolist()
        for genome in latest_assembly_versions:
            if genome not in local_genomes:
                new_genomes.append(genome)

    return new_genomes

def remove_old_genomes(genbank_mirror, assembly_summary, old_genomes, logger):

    # TODO: there might be a faster way to do this with pandas
    for genome_id in old_genomes:
        # Would have to keep the old assembly summary file in order to avoid globbing the species dir
        associated_files = glob.glob("{}/*/{}*".format(genbank_mirror, genome_id)) # globs sketch files as well
        for f in associated_files:
            os.remove(f)
            logger.info("Removed {}".format(f))

def get_sketch_files(genbank_mirror):

    sketch_files = []
    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if re.match('GCA.*msh', f):
                genome_id = re.sub('\.msh', '', f)
                sketch_files.append(genome_id)

    return sketch_files

def get_missing_sketch_files(local_genomes, sketch_files):

    missing_sketch_files = set(local_genomes) - set(sketch_files)

    return missing_sketch_files

def get_old_genomes(genbank_mirror, assembly_summary, local_genomes):

    old_genomes = []

    latest_assembly_versions = assembly_summary.index.tolist()
    for genome_id in local_genomes:
        if genome_id not in latest_assembly_versions:
            old_genomes.append(genome_id)

    return old_genomes

def assess_genbank_mirror(genbank_mirror, assembly_summary, species_list, logger):

    local_genomes = get_local_genomes(genbank_mirror)
    # new_genomes = get_new_genome_list(genbank_mirror, assembly_summary, local_genomes, species_list)
    sketch_files = get_sketch_files(genbank_mirror)
    missing_sketch_files = get_missing_sketch_files(local_genomes, sketch_files)

    logger.info('{} sketch files present in local collection.'.format(len(sketch_files)))
    logger.info('{} sketch files not in local collection.'.format(len(missing_sketch_files)))

    return local_genomes, sketch_files, missing_sketch_files

def unzip_genome(root, f, genome_id):

    """
    Decompress genome and remove the compressed genome.
    """

    zipped_src = os.path.join(root, f)
    zipped = gzip.open(zipped_src)
    decoded = zipped.read()
    unzipped = "{}.fasta".format(genome_id)
    unzipped = os.path.join(root, unzipped)
    unzipped = open(unzipped, "wb")
    zipped.close()
    os.remove(zipped_src)
    unzipped.write(decoded)
    unzipped.close()

def unzip_genbank_mirror(genbank_mirror):

    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if f.endswith("gz"):
                genome_id = "_".join(f.split("_")[:2])
                try:
                    unzip_genome(root, f, genome_id)
                except OSError:
                    continue

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    path_vars = config.instantiate_path_vars(genbank_mirror)
    clean_up(genbank_mirror, path_vars)

if __name__ == "__main__":
    main()
