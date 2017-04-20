#!/usr/bin/env python

import os
import argparse
from urllib.request import urlretrieve
from urllib.error import URLError
from ftplib import error_temp
from time import strftime, sleep

def grab_zipped_genome(genbank_mirror, species, genome_id, genome_url, ext=".fna.gz"):

    """
    Download compressed genome from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
    """

    zipped_path = "{}_genomic{}".format(genome_id, ext)
    zipped_url = "{}/{}".format(genome_url, zipped_path)
    zipped_dst = os.path.join(genbank_mirror, species, zipped_path)
    urlretrieve(zipped_url, zipped_dst)

def get_genome_id_and_url(assembly_summary, accession):

    genome_id = assembly_summary.ftp_path[accession].split('/')[-1]
    genome_url = assembly_summary.ftp_path[accession]

    return genome_id, genome_url

def sync_latest_genomes(genbank_mirror, assembly_summary, new_genomes, logger):

    for accession in new_genomes:
        genome_id, genome_url = get_genome_id_and_url(assembly_summary, accession)
        species = assembly_summary.scientific_name.loc[accession]
        try:
            grab_zipped_genome(genbank_mirror, species, genome_id, genome_url)
            logger.info("Downloaded {}".format(genome_id))
        except error_temp as e:
            logger.exception('error_temp for {}\n{}'.format(genome_id, e))
            sleep(2)
            grab_zipped_genome(genbank_mirror, species, genome_id, genome_url)
            logger.info("Downloaded {}".format(genome_id))
        except URLError as e:
            logger.exception('URLError[1] for {}\n{}'.format(genome_id, e))
            grab_zipped_genome(genbank_mirror, species, genome_id, genome_url, ext=".fasta.gz")
            logger.info("Downloaded {}".format(genome_id))
        except URLError as e:
            logger.exception('URLError[2] for {}\n{}'.format(genome_id, e))
            continue

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-u", "--unzip", action="store_true")
    args = parser.parse_args()
    genbank_mirror = args.genbank_mirror

    assembly_summary = get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt")
    unzip_genbank_mirror(genbank_mirror)
    rename(genbank_mirror, assembly_summary)

if __name__ == "__main__":
    main()
