#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --account=s3673
#SBATCH --constraint='sky|hasw'
#SBATCH --output=diskframe-%j.log

set -e
source ~/.bashrc
mod_r

Rscript scripts/disk-frame.R
