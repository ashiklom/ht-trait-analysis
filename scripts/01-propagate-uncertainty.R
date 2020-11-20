#!/usr/bin/env Rscript

library(fs)
library(stringr)
library(magrittr)
library(foreach)
library(future)
library(doFuture)

stopifnot(
  requireNamespace("qs", quietly = TRUE),
  requireNamespace("hyperSpec", quietly = TRUE),
  requireNamespace("caTools", quietly = TRUE)
)

argv <- commandArgs(trailingOnly = TRUE)
resultdir <- argv[1]
stopifnot(!is.na(resultdir))

reflfiles <- dir_ls(
  resultdir,
  regexp = "estimated-reflectance(-[[:digit:]]+)?$",
  recurse = TRUE
)
uncfiles <- reflfiles %>%
  str_replace("estimated-uncertainty", "posterior-uncertainty")

reflectance_uncertainty_fast <- function(reflfile, uncfile,
                                         n = 20,
                                         overwrite = FALSE) {

  outfile <- str_replace(
    reflfile,
    "estimated-reflectance",
    "propagated-uncertainty"
  ) %>% str_c(".qs")
  if (overwrite) suppressWarnings(file.remove(outfile))
  if (file.exists(outfile)) {
    ## message("Skipping existing file ", outfile)
    return(outfile)
  }
  refl_r <- hyperSpec::read.ENVI(reflfile, paste0(reflfile, ".hdr"))
  rwl <- hyperSpec::wl(refl_r)
  refl <- refl_r$spc
  unc <- caTools::read.ENVI(uncfile, paste0(uncfile, ".hdr"))[1, , seq_along(rwl)]

  tryCatch(
    stopifnot(
      length(dim(refl)) == 2,
      all(dim(unc) == dim(refl))
    ),
    error = function(e) {
      stop("Encountered an error:",
           "\nRefl file: ", reflfile,
           "\nUnc file: ", uncfile,
           "\nError:\n",
           conditionMessage(e))
    }
  )
  nobs <- nrow(refl)
  nwl <- ncol(refl)
  zvals <- qnorm(seq(0.025, 0.975, length.out = n))
  # NOTE: Same order as input
  rsamp <- array(NA_real_, c(n, nobs, nwl))
  for (i in seq_along(zvals)) {
    rsamp[i,,] <- pmax(unc * zvals[i] + refl, 0)
  }

  # Observations, samples, wavelengths
  result <- aperm(rsamp, c(2, 1, 3))
  attr(result, "wavelength") <- rwl
  qs::qsave(result, outfile)
  outfile
}

plan("multicore")
registerDoFuture()
outfiles <- foreach(reflfile = reflfiles, uncfile = uncfiles,
                   .errorhandling = "remove",
                   .combine = c,
                   .packages = c("stringr", "magrittr")) %dopar%
  reflectance_uncertainty_fast(reflfile, uncfile)

out_norm <- normalizePath(outfiles)
writeLines(out_norm, "./propagated-uncertainty-files.txt")
