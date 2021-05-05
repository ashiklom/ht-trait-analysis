# Hypertrace paper and code

## Setting up the simulations

Run the following:

#. `Rscript scripts/prepare-chadwick-inputs.R` -- Converts Chadwick NEON data into ENVI files. Stored in `data/derived/ht-trait-inputs`
#. `Rscript scripts/prepare-hypertrace.R` -- Creates the Hypertrace configuration files. Stored in `data/derived/ht-trait-inputs`
#. `rsync -avz --progress data/derived/ht-trait-inputs discover:~/projects/sbg-uncertainty/isofit/examples/py-hypertrace/` -- Copy data over to cluster
#. (On cluster, from `py-hypertrace` directory above) `sbatch ht-trait-inputs/run_chadwick.sh`
