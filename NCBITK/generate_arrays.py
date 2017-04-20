import os
import argparse
import subprocess
import pandas as pd
import NCBITK.config as config
from re import sub
from time import strftime, sleep

ymd = strftime("%y.%m.%d")

def gen_sbatch_script(genbank_mirror, array, job_name, time):
    None

def gen_sbatch_array_script(genbank_mirror, array, job_name, mem, time, chunk=False):

    info_dir, slurm, out, log_file = config.instantiate_path_vars(genbank_mirror)
    out = os.path.join(out, "{}_%A_%a.out".format(job_name))
    print(job_name)
    sbatch_script = os.path.join(slurm, "{}.sbatch".format(job_name))
    print('Generating {}'.format(sbatch_script))
    array_len = len(list(open(array)))
    with open(sbatch_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time={}\n".format(time))
        f.write("#SBATCH --mem={}\n".format(mem))
        f.write("#SBATCH --job-name={}\n".format(job_name))
        f.write("#SBATCH --output={}\n".format(out)) # can I remove this line to avoid getting out files?
        if chunk:
            f.write("#SBATCH --array=1-{}%3\n".format(array_len))
        else:
            f.write("#SBATCH --array=1-{}\n".format(array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(array))
        f.write("srun $cmd")

    return sbatch_script

def gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array):

    info_dir, slurm, out, log_file = instantiate_path_vars(genbank_mirror)
    out = os.path.join(out, "get_latest_%a.out")
    latest_assembly_versions_script = os.path.join(slurm, "get_latest_assembly_versions.sbatch")
    print('Generating {}'.format(latest_assembly_versions_script))
    array_len = len(list(open(latest_assembly_versions_array)))
    with open(latest_assembly_versions_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=01:00\n")
        f.write("#SBATCH --job-name=get_latest\n")
        f.write("#SBATCH --output={}\n".format(out)) # can I remove this line to avoid getting out files?
        f.write("#SBATCH --array=1-{}%2\n".format(array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(latest_assembly_versions_array))
        f.write("srun $cmd")

    return latest_assembly_versions_script

def gen_sync_array_script(genbank_mirror, get_latest_job_id):
    info_dir, slurm, out, log_file = instantiate_path_vars(genbank_mirror)
    sync_array_script = os.path.join(slurm, 'sync_array_script.sbatch')
    print('Generating {}'.format(sync_array_script))

    # consolidate into function
    with open(sync_array_script, 'a') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=01:00\n")
        f.write("#SBATCH --job-name=gen_sync_array\n")
        f.write("#SBATCH --output={}/gen_sync_array.out\n".format(out))
        f.write("#SBATCH --dependency={}\n".format(get_latest_job_id))
        f.write('cmd="python /common/contrib/tools/NCBITK/slurm/generate_arrays.py {}"\n'.format(genbank_mirror))
        f.write("srun $cmd")

    return sync_array_script

def gen_grab_genomes_script(genbank_mirror, sync_array_job_id):

    info_dir, slurm, out, log_file = instantiate_path_vars(genbank_mirror)
    sync_array = os.path.join(slurm, "sync_array.txt")
    grab_genomes_script = os.path.join(slurm, 'grab_genomes_script.sbatch')
    out = os.path.join(out, 'grab_genomes%a.out')
    print('Generating {}'.format(grab_genomes_script))

    while not os.path.isfile(sync_array):
        sleep(10)
    else:
        sync_array_len = len(list(open(sync_array)))

    with open(grab_genomes_script, 'a') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=04:30\n")
        f.write("#SBATCH --job-name=grab_genomes\n")
        f.write("#SBATCH --output={}\n".format(out))
        f.write("#SBATCH --dependency={}\n".format(sync_array_job_id))
        f.write("#SBATCH --array=1-{}%2\n".format(sync_array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(sync_array))
        f.write("srun $cmd")

    return grab_genomes_script, sync_array_len

def write_grab_genomes_array(genbank_mirror):

    latest_assembly_versions = curate.read_latest_assembly_versions(genbank_mirror)
    grab_genomes_array = os.path.join(genbank_mirror, ".info", "slurm", "grab_genomes_array.txt")
    print('Generating {}'.format(grab_genomes_array))
    new_genomes = curate.get_new_genomes(genbank_mirror, latest_assembly_versions)
    args = []
    for name in new_genomes:
        species = latest_assembly_versions.loc[name, 'species']
        path = latest_assembly_versions.loc[name, 'dir']
        args.append(','.join([species, name, path]))

    groups = [args[n:n+2000] for n in range(0, len(args), 2000)]
    with open(grab_genomes_array, "a") as f:
        for group in groups:
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/ftp_functions.py -g {} {}\n".format(genbank_mirror, ' '.join(group)))

    return grab_genomes_array

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    parser.add_argument("-g", '--grab_genomes', action='store_true')
    args = parser.parse_args()

    if args.grab_genomes:
        write_grab_genomes_array(args.genbank_mirror)

if __name__ == "__main__":
    main()
