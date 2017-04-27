#!/usr/bin/env python

import os
import logging
import subprocess
import pandas as pd
import tarfile
from urllib.request import urlretrieve

def get_assembly_summary(genbank_mirror, update=True, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"):

    """Get current version of assembly_summary.txt and load into DataFrame"""

    assembly_summary_dst = os.path.join(genbank_mirror, ".info", "assembly_summary.txt")

    if update:
        urlretrieve(assembly_summary_url, assembly_summary_dst)
        assembly_summary = pd.read_csv(assembly_summary_dst, sep="\t", index_col=0, skiprows=1)

    else:
        assembly_summary = pd.read_csv(assembly_summary_dst, sep="\t", index_col=0)


    return assembly_summary

def update_assembly_summary(genbank_mirror, assembly_summary, names):

    for taxid in names.index:
        scientific_name = names.scientific_name.loc[taxid]
        # get the list of indices that share the same species_taxid in assembly_summary
        ixs = assembly_summary.index[assembly_summary.species_taxid == taxid].tolist()
        assembly_summary.loc[ixs, 'scientific_name'] = scientific_name

    updated_assembly_summary = os.path.join(genbank_mirror, '.info', 'assembly_summary.txt')
    assembly_summary.to_csv(updated_assembly_summary, sep='\t')

    return assembly_summary

def get_scientific_names(genbank_mirror, assembly_summary, taxdump_url="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", update=True):

    """
    Get names.dmp from the taxonomy dump
    """

    info_dir = os.path.join(genbank_mirror, ".info")
    names_dmp = os.path.join(genbank_mirror, ".info", 'names.dmp')

    if update:
        taxdump = urlretrieve(taxdump_url)
        taxdump_tar = tarfile.open(taxdump[0])
        taxdump_tar.extract('names.dmp', info_dir)

    sed_cmd = "sed -i '/scientific name/!d' {}".format(names_dmp) # we only want rows with the scientific name
    subprocess.Popen(sed_cmd, shell='True').wait()
    names = pd.read_csv(names_dmp, sep='\t', index_col=0, header=None, usecols=[0,2])
    names = names.loc[set(assembly_summary.species_taxid.tolist())]
    names.index.name = 'species_taxid'
    names.columns = ['scientific_name']
    names.scientific_name.replace({' ': '_'}, regex=True, inplace=True)
    names.scientific_name.replace({'/': '_'}, regex=True, inplace=True)
    names.to_csv(names_dmp)

    return names

def get_resources(genbank_mirror, logger, update):

    """
    Get assembly_summary.txt for bacteria and taxonomy dump file.
    Parse and load into Pandas DataFrames.
    """

    if update:
        assembly_summary = get_assembly_summary(genbank_mirror, update)
        logger.info('Got new assembly_summary.txt')

        names = get_scientific_names(genbank_mirror, assembly_summary)
        assembly_summary = update_assembly_summary(genbank_mirror, assembly_summary, names)
    else:
        assembly_summary = get_assembly_summary(genbank_mirror, update)

    return assembly_summary
