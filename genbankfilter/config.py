#!/usr/bin/env python

import os
import time
import logging

def instantiate_path_vars(genbank_mirror):

    info_dir = os.path.join(genbank_mirror, ".info", 'filtering')
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out")
    ymd = time.strftime("%y.%m.%d_%I:%M:%S_%p")

    for d in [info_dir, slurm, out]:
        if not os.path.isdir(d):
            os.mkdir(d)

    log_file = os.path.join(info_dir, 'log_{}.txt'.format(ymd))
    logger = instantiate_logger(log_file)

    return info_dir, slurm, out, logger

def instantiate_logger(log_file):

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    fh.setFormatter(formatter)

    logger.addHandler(fh)

    return logger
