#!/usr/bin/env Rscript

library(stringr)
library(foreach)
library(fs)
library(future)
library(doFuture)

plsr_list <- function(plsr_coef_file) {
  plsr_data <- read.csv(plsr_coef_file)
  pls_wl <- names(plsr_data) %>%
    stringr::str_subset("WVL_") %>%
    stringr::str_remove("WVL_") %>%
    as.numeric()

  pls_int <- plsr_data[, "Intercept"]
  pls_slope <- plsr_data %>%
    dplyr::select(dplyr::starts_with("WVL_")) %>%
    as.matrix() %>%
    t()
  stopifnot(
    ncol(pls_slope) == length(pls_int),
    nrow(pls_slope) == length(pls_wl)
  )
  list(intercept = pls_int,
       slope = pls_slope,
       wl = pls_wl)
}

fastquant <- function(x, probs) {
  td <- tdigest::tdigest(x)
  tdigest::tquantile(td, probs)
}

spec_plsr <- function(specdata, plsr_list) {

  # Observations, samples, wavelengths
  bdat <- specdata
  bwl <- attr(bdat, "wavelength")

  # Define dimensions, to make things clearer
  nobs <- nrow(bdat)
  nsamp <- ncol(bdat)
  nwl <- length(bwl)

  pls_wl <- plsr_list[[1]][["wl"]]

  # Resample reflectance samples to AVIRIS-Classic wavelengths
  interp_wl <- function(x) {
    tryCatch(
      approx(bwl, x, pls_wl, yleft = 0, yright = 0)$y,
      error = function(e) {
        message(conditionMessage(e))
        return(rep_len(NA_real_, length(pls_wl)))
      }
    )
  }
  # Result: wavelength, nsamp, nobs
  bdat_r <- apply(bdat, c(1, 2), interp_wl) * 10000

  out_array <- array(NA_real_, c(length(plsr_list), nobs, nsamp, 6))
  dimnames(out_array) <- list(names(plsr_list), NULL, NULL, c(
    "Mean", "SD", "Min", "Max", "q025", "q975"
  ))

  for (j in seq_along(plsr_list)) {
    pls_int <- plsr_list[[j]][["intercept"]]
    pls_slope <- plsr_list[[j]][["slope"]]
    npls <- length(pls_int)

    # npls x nobs x nsamp
    pls_est <- array(NA_real_, c(npls, nobs, nsamp))
    for (i in seq_len(nsamp)) {
      pls_est[,,i] <- pls_int + t(bdat_r[,,i]) %*% pls_slope
    }

    # Each of these should be nobs x nsamp; drop first dimension
    pls_est_means <- apply(pls_est, c(2, 3), mean, na.rm = TRUE)
    pls_est_sd <- apply(pls_est, c(2, 3), sd, na.rm = TRUE)
    pls_est_min <- apply(pls_est, c(2, 3), min, na.rm = TRUE)
    pls_est_max <- apply(pls_est, c(2, 3), max, na.rm = TRUE)
    pls_est_q <- apply(pls_est, c(2, 3), fastquant, c(0.025, 0.975))

    out_array[j,,,"Mean"] <- pls_est_means
    out_array[j,,,"SD"] <- pls_est_sd
    out_array[j,,,"Min"] <- pls_est_min
    out_array[j,,,"Max"] <- pls_est_max
    out_array[j,,,"q025"] <- pls_est_q[1,,]
    out_array[j,,,"q975"] <- pls_est_q[2,,]
  }

  out_array
}

do_plsr <- function(reflfile, overwrite = FALSE) {
  outfile <- reflfile %>%
    stringr::str_replace("propagated-uncertainty", "plsr-results")
  if (overwrite) fs::file_delete(outfile)
  if (fs::file_exists(outfile)) {
    return(outfile)
  }
  pls_coef_files <- fs::dir_ls(
    "data/plsr-coefficients/singh-ea-supplement",
    regexp = "PLS_coefficients_Mode_A"
  )
  pls_tags <- fs::path_file(pls_coef_files) %>%
    fs::path_ext_remove() %>%
    stringr::str_remove("PLS_coefficients_Mode_A_")
  pls_list <- lapply(pls_coef_files, plsr_list)
  specdata <- qs::qread(reflfile)
  output <- spec_plsr(specdata, pls_list)
  names(output) <- pls_tags
  qs::qsave(output, outfile)
  outfile
}

prop_files <- readLines("./propagated-uncertainty-files.txt")

do_plsr(prop_files[1])

plan("multicore")
registerDoFuture()
resultfiles <- foreach(reflfile = prop_files,
                       .combine = c,
                       .packages = c("magrittr"),
                       .errorhandling = "remove") %dopar%
  do_plsr(reflfile)

result_norm <- normalizePath(resultfiles)
writeLines(resultfiles, "./plsr-output-files.txt")
