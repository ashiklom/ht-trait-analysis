reflectance_uncertainty <- function(reflfile, uncfile, outfile, n = 10) {
  refl_r <- brick(reflfile)
  unc_r <- brick(uncfile)

  rwl <- get_raster_wl(refl_r)
  ## ibad <- (rwl > 1300 & rwl < 1450) |
  ##   (rwl > 1780 & rwl < 1950)
  refl <- refl_r[]
  stopifnot(length(dim(refl)) == 2)
  unc <- unc_r[][,seq_along(rwl)]
  stopifnot(all(dim(unc) == dim(refl)))

  rcorr <- cor(refl)

  rsamp <- array(NA_real_, c(length(rwl), n, nrow(refl)))
  pb <- progress_bar$new(total = nrow(refl))
  for (i in seq_len(nrow(refl))) {
    rsamp[,,i] <- tryCatch(
      sample_reflectance(refl[i,], unc[i,], rcorr, n = n),
      error = function(e) {
        warning(conditionMessage(e))
        return(NA_real_)
      }
    )
    pb$tick()
  }

  rsamp_r <- brick(aperm(rsamp, c(3, 2, 1)))
  writeRaster(rsamp_r, outfile, format = "ENVI",
              options = c("INTERLEAVE=BIP", "SUFFIX=ADD"),
              overwrite = TRUE)
  write(c("wavelength={", paste(rwl, collapse = ",\n"), "}"),
        paste0(outfile, ".hdr"),
        append = TRUE)
  write(c("bbl={", paste(rep_len(0, length(rwl)), collapse = ",\n"), "}"),
        paste0(outfile, ".hdr"),
        append = TRUE)
  invisible(outfile)
}

# R version
# 67.06

# Cpp version
# 60.92

## sourceCpp("R/cpp-functions.cpp")

rmvnorm <- function(n, mu, Sigma) {
  L <- t(chol(Sigma))
  Z <- matrix(rnorm(length(mu) * n), ncol = n)
  mu + L %*% Z
}

sample_reflectance <- function(refl, unc, cormat, n = 1) {
  dunc <- diag(unc)
  covmat <- dunc %*% cormat %*% dunc
  rmvnorm(n, refl, covmat)
}

if (FALSE) {

  refldir <- "../isofit-calibration/examples/py-hypertrace/zz-calibration/output/neon/"
  reflfiles <- dir_ls(
    refldir,
    regexp = "estimated-reflectance-[[:digit:]]+$",
    recurse = TRUE
  )
  uncfiles <- dir_ls(
    refldir,
    regexp = "posterior-uncertainty-[[:digit:]]+$",
    recurse = TRUE
  )

  reflfile <- reflfiles[1]
  uncfile <- uncfiles[1]
  outfile <- str_replace(reflfile, "estimated-reflectance", "uncertain-reflectance")

}
