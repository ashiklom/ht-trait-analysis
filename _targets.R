## library() calls go here
library(conflicted)
## library(dotenv)
library(targets)
library(fs)
library(stringr)
library(dplyr)
library(raster)
library(tidyr)
library(purrr)
library(progress)
library(ggplot2)
library(ggthemes)
library(distributions3)
library(data.table)
library(forcats)

stopifnot(
  requireNamespace("reshape2", quietly = TRUE),
  requireNamespace("qs", quietly = TRUE),
  requireNamespace("hyperSpec", quietly = TRUE)
)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("between", "dplyr")
conflict_prefer("pdf", "distributions3")
conflict_prefer("quantile", "stats")

resultdir <- path("data", "derived", "ht-output")

# Source all R files
dir_ls("./R/", regexp = "\\.R$") %>%
  str_subset("plan[.]R", negate = TRUE) %>%
  walk(source)

tar_option_set(format = "qs")

tar_pipeline(
  tar_target(
    plsr_files,
    dir_ls(resultdir, regexp = "plsr-results", recurse = TRUE)
  ),
  tar_target(plsr_info, as.data.table(parse_ht_directory(plsr_files))),
  tar_target(raw_plsr_results, rbindlist(lapply(plsr_files, read_plsr_file)),
             format = "fst_dt"),
  tar_target(
    true_plsr,
    do_true_plsr("data/derived/zz-traitanalysis/neon-vegetation.envi")
  )
)
