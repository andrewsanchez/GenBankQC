from NCBITK import config
from NCBITK import curate
from NCBITK import sync
from NCBITK import get_resources

import unittest
import os
import glob
import tempfile
import shutil
import pandas as pd


class TestCurate(unittest.TestCase):

    def setUp(self):

        self.genbank_mirror = tempfile.mkdtemp(prefix='Genbank_')
        self.assembly_summary = pd.read_csv('NCBITK/test/resources/assembly_summary.txt', sep="\t", index_col=0)
        self.path_vars = config.instantiate_path_vars(self.genbank_mirror)
        self.info_dir, self.slurm, self.out, self.logger = self.path_vars

        self.test_species = 'Acinetobacter_nosocomialis'
        self.test_genomes = self.assembly_summary.index[self.assembly_summary.scientific_name == self.test_species]
        self.species_list = curate.get_species_list(self.assembly_summary, [self.test_species])
        self.species_dir = os.path.join(self.genbank_mirror, self.test_species)

        self.genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        self.local_genomes, self.new_genomes,\
        self.old_genomes, self.sketch_files,\
        self.missing_sketch_files = self.genbank_assessment

    def test_create_species_dirs_all(self):

        species_list = curate.get_species_list(self.assembly_summary, 'all')
        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, species_list)
        local_species = os.listdir(self.genbank_mirror)
        local_species.remove('.info')

        self.assertEqual(len(local_species), len(species_list))

    def test_create_species_dirs_list(self):

        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, self.species_list)
        local_species = os.listdir(self.genbank_mirror)
        local_species.remove('.info')

        self.assertEqual(len(local_species), len(self.species_list))

    def test_create_species_dirs_str(self):
        None
        

    def test_assess_fresh(self):

        self.assertTrue(len(self.new_genomes) == len(self.test_genomes))
        self.assertTrue(len(self.missing_sketch_files) == len(self.test_genomes))
        self.assertTrue(len(self.sketch_files) == 0)
        self.assertTrue(len(self.local_genomes) == 0)
        self.assertTrue(len(self.old_genomes) == 0)

    def test_sync_latest_genomes(self):

        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, self.species_list)

        genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment
        before_sync_new_genomes = new_genomes

        sync.sync_latest_genomes(self.genbank_mirror, self.assembly_summary, new_genomes, self.logger)

        species_dir = os.path.join(self.genbank_mirror, self.test_species)

        for f in glob.glob('{}/GCA*'.format(species_dir)):
            self.assertTrue(os.path.isfile(f))

        genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment

        self.assertTrue(len(before_sync_new_genomes) == len(local_genomes))
        self.assertTrue(len(local_genomes) == len(self.test_genomes))
        self.assertTrue(len(missing_sketch_files) == len(self.test_genomes))
        self.assertTrue(len(new_genomes) == 0)
        self.assertTrue(len(old_genomes) == 0)
        self.assertTrue(len(sketch_files) == 0)

    def test_assess_partial(self):

        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, self.species_list)

        genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment
        before_sync_new_genomes = new_genomes[:10]
        after_sync_new_genomes = new_genomes[10:]

        # TODO: Already testing this function.  Probably unnecessary
        sync.sync_latest_genomes(self.genbank_mirror, self.assembly_summary, before_sync_new_genomes, self.logger)

        genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment

        self.assertTrue(len(local_genomes) == len(before_sync_new_genomes))
        self.assertTrue(len(missing_sketch_files) == len(before_sync_new_genomes) + len(new_genomes))
        self.assertFalse(len(new_genomes) == 0)
        self.assertTrue(len(new_genomes) == len(after_sync_new_genomes))
        self.assertTrue(len(old_genomes) == 0)
        self.assertTrue(len(sketch_files) == 0)

    def test_get_old_genomes(self):

        local_genomes = self.test_genomes
        not_in_assembly_summary = self.test_genomes[:5].tolist()
        self.assembly_summary.drop(not_in_assembly_summary, inplace=True)

        old_genomes = curate.get_old_genomes(self.genbank_mirror, self.assembly_summary, local_genomes)

        self.assertTrue(sorted(not_in_assembly_summary) == sorted(old_genomes))

    def test_get_local_genomes(self):

        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, self.species_list)

        for genome in self.test_genomes:
            dst = os.path.join(self.species_dir, genome)
            tempfile.mkstemp(prefix=genome, dir=self.species_dir)

        local_genomes = curate.get_local_genomes(self.genbank_mirror)

        self.assertTrue(len(local_genomes) == len(self.test_genomes))

    def tearDown(self):
        shutil.rmtree(self.genbank_mirror)

class TestArrays(unittest.TestCase):
    None

if __name__ == '__main__':
    unittest.main()
