#!/usr/bin/env python

import os, re, argparse
import pandas as pd
from urllib.request import urlretrieve
from .get_resources import get_assembly_summary

def rm_duplicates(seq):

    """
    remove duplicate strings during renaming
    """

    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def rename(target_dir, assembly_summary):

    """
    Clean up assembly_summary.txt and renamed FASTA's.
    """

    # If infraspecific_name and isolate columns are empty, fill infraspecific_name with "NA"
    assembly_summary.update(assembly_summary['infraspecific_name'][(assembly_summary['infraspecific_name'].isnull()) &\
            (assembly_summary['isolate'].isnull())].fillna('NA'))

    # If infraspecific_name column is empty and isolate column is not empty, fill infraspecific_name with the value of isolate.
    assembly_summary.update(assembly_summary['infraspecific_name'][(assembly_summary['infraspecific_name'].isnull()) &\
            (assembly_summary['isolate'].notnull())].fillna(assembly_summary['isolate']))

    assembly_summary.assembly_level.replace({' ': '_'}, regex=True, inplace=True)
    assembly_summary.organism_name.replace({' ': '_'}, regex=True, inplace=True)
    assembly_summary.organism_name.replace({'[\W]': '_'}, regex=True, inplace=True)
    assembly_summary.infraspecific_name.replace({'[\W]': '_'}, regex=True, inplace=True)

    for root, dirs, files in os.walk(target_dir):
        for f in files:
            if f.startswith("GCA"):
                # accession_id = "_".join(f.split('_')[0:2])
                accession_id = re.sub('.fasta', '', f)
                if accession_id in assembly_summary.index:
                    org_name = assembly_summary.get_value(accession_id, 'organism_name')
                    strain = assembly_summary.get_value(accession_id, 'infraspecific_name')
                    assembly_level  = assembly_summary.get_value(accession_id, 'assembly_level')
                    new_name = '{}_{}_{}_{}.fasta'.format(accession_id, org_name, strain, assembly_level)
                    rm_words = re.compile( r'((?<=_)(sp|sub|substr|subsp|str|strain)(?=_))' )
                    new_name = rm_words.sub('_', new_name)
                    new_name = re.sub(r'_+', '_', new_name )
                    new_name = new_name.split('_')
                    new_name = rm_duplicates(new_name)
                    new_name = '_'.join(new_name)
                    old = os.path.join(root, f)
                    new = os.path.join(root, new_name)
                    os.rename(old, new)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('target_dir', help = 'The folder whose contents will be renamed', type=str)
    parser.add_argument('-s', '--source', help = 'Specify a directory to rename.', action="store_true")
    args = parser.parse_args()
    genbank_mirror = args.target_dir

    rename(genbank_mirror, assembly_summary)

if __name__ == '__main__':
    main()
