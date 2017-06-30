import os, argparse
import pandas as pd
from subprocess import Popen

def reference_genome(sketch_files, reference_genome, threshold,  mash_exe):
    output = reference_genome.strip(".fasta.msh")
    output = "{}_distances.csv".format(output)
    output = os.path.join(sketch_files, output)
    distance_cmd = "{} dist -t {} *msh > {}".format(mash_exe, reference_genome, output)
    Popen(distance_command, shell="True").wait()
    distances = pd.read_csv(output, delimiter="\t", index_col=0, header=0)
    passed = distances[distances <=  threshold]
    passed_log = os.path.join(sketch_files, "{}_{}".format(reference_genome, threshold))
    passed.to_csv(passed_log)

def all(genbank_mirror, species_list='all'):

    """
    Create a distance matrix for all genomes of a given species or for a list of species.
    Defaults to all species in genbank_mirror
    """

    if species_list == 'all':
        for root, dirs, files in os.walk(genbank_mirror):
            for d in dirs:
                out_file = os.path.join(root, d, 'all.msh')
                all_sketch_files = os.path.join(root, d, '*msh')
                if os.path.isfile(out_file):
                    os.remove(out_file)
                paste_cmd = "mash paste {} {}".format(out_file, all_sketch_files)
                Popen(paste_cmd, shell="True")
                print("Created {}".format(out_file))
    else:
        for species in species_list:
            None


def main():
    parser = argparse.ArgumentParser(description = "Identify genomes within a provided distance from a single reference genome")
    parser.add_argument("sketch_files", help = "directory containing MASH sketch files")
    parser.add_argument("reference_genome", help = "directory containing MASH sketch files")
    parser.add_argument("-t", "--threshold", help = "Find genomes that are within this distance", default=.02, type=float)
    parser.add_argument("-x", "--mash_exe", help = "Path to MASH exectuable", default="/common/contrib/bin/mash-Linux64-v1.1.1/mash")
    args = parser.parse_args()

    reference_genome(args.sketch_files, args.reference_genome, args.threshold, args.mash_exe)

if __name__ == "__main__":
    main()
