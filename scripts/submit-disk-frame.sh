#!/usr/bin/env bash
#SBATCH --

. "$HOME/.basrc"
mod_r

Rscript scripts/disk-frame.R
