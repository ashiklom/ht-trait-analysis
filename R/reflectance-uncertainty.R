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

  # Wavelength, samples, observations
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

  # Observations, samples, wavelengths
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

reflectance_uncertainty_fast <- function(reflfile, uncfile, n = 20) {
  refl_r <- hyperSpec::read.ENVI(reflfile, paste0(reflfile, ".hdr"))
  rwl <- hyperSpec::wl(refl_r)
  refl <- refl_r$spc
  unc <- caTools::read.ENVI(uncfile, paste0(uncfile, ".hdr"))[1, , seq_along(rwl)]

  stopifnot(length(dim(refl)) == 2)
  stopifnot(all(dim(unc) == dim(refl)))
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
  result
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
