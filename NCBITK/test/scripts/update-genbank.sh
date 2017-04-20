#!/bin/bash

species='Acinetobacter_nosocomialis'
python ~/projects/NCBITK/run.py --update "$@" --use_local --species $species
