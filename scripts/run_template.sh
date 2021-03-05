#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --account=s3673
#SBATCH --constraint='sky|hasw'
#SBATCH --array=1-@NRUNS@
#SBATCH --output=logs/hta-%A_%a.log

set -e
source ~/.bashrc

conda activate isofit
python workflow.py $(printf "./ht-trait-inputs/configs/chadwick-%03d.json" $SLURM_ARRAY_TASK_ID)
