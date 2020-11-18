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
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(ggthemes)
library(distributions3)
library(data.table)
library(forcats)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("between", "dplyr")
conflict_prefer("pdf", "distributions3")
conflict_prefer("quantile", "stats")

tar_option_set(format = "qs")
options(clustermq.scheduler = "multicore")
future::plan("multicore")

resultdir <- path("data", "derived", "ht-output")

dir_ls("./R/", regexp = "\\.R$") %>%
  str_subset("plan[.]R", negate = TRUE) %>%
  walk(source)

tar_pipeline(
  tar_target(
    refl_paths,
    dir_ls(resultdir, regexp = "estimated-reflectance(-[[:digit:]]+)?$",
           recurse = TRUE)
  ),
  tar_target(
    unc_paths,
    dir_ls(resultdir, regexp = "posterior-uncertainty(-[[:digit:]]+)?$",
           recurse = TRUE)
  ),
  tar_target(reflfiles, refl_paths, format = "file", pattern = map(refl_paths),
             deployment = "main"),
  tar_target(uncfiles, unc_paths, format = "file", pattern = map(unc_paths),
             deployment = "main"),
  tar_target(
    propagated_refl,
    reflectance_uncertainty_fast(reflfiles, uncfiles),
    pattern = map(reflfiles, uncfiles),
    iteration = "list"
  ),
  tar_target(pls_coef_paths, dir_ls(
    "data/plsr-coefficients/singh-ea-supplement",
    regexp = "PLS_coefficients_Mode_A"
  )),
  tar_target(pls_coef_files, pls_coef_paths,
             pattern = map(pls_coef_paths), format = "file"),
  tar_target(plsr_data, read.csv(pls_coef_files),
             pattern = map(pls_coef_files), iteration = "list"),
  tar_target(
    plsr_result,
    spec_plsr(propagated_refl, plsr_data),
    pattern = cross(propagated_refl, plsr_data),
    iteration = "list"
  )
)
