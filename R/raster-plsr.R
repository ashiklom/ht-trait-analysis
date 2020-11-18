raster_plsr <- function(rasterfile, outfile, plsr_coef_file) {
  plsr_data <- read.csv(plsr_coef_file)

  pls_wl <- names(plsr_data) %>%
    str_subset("WVL_") %>%
    str_remove("WVL_") %>%
    as.numeric()

  pls_int <- plsr_data[, "Intercept"]
  pls_slope <- plsr_data %>%
    select(starts_with("WVL_")) %>%
    as.matrix() %>%
    t()

  stopifnot(
    ncol(pls_slope) == length(pls_int),
    nrow(pls_slope) == length(pls_wl)
  )

  b <- brick(rasterfile)
  bwl <- get_raster_wl(b)

  bdat <- b[]
  # Row order is b[1,1], b[1,2], b[1,3]...b[2,1], b[2,2], b[2,3]...
  # I.e., Reflectance samples are kept together

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
  bdat_r <- t(apply(bdat, 1, interp_wl)) * 10000

  # Rows = Model samples; Cols = pixels
  pls_est <- pls_int + t(bdat_r %*% pls_slope)
  pls_est_means <- colMeans(pls_est, na.rm = TRUE)
  pls_est_sd <- apply(pls_est, 2, sd, na.rm = TRUE)
  pls_est_min <- apply(pls_est, 2, min, na.rm = TRUE)
  pls_est_max <- apply(pls_est, 2, max, na.rm = TRUE)
  pls_est_lo <- apply(pls_est, 2, quantile, 0.025, na.rm = TRUE)
  pls_est_hi <- apply(pls_est, 2, quantile, 0.975, na.rm = TRUE)

  out_array <- array(NA_real_, c(ncol(b), nrow(b), 6))
  dimnames(out_array) <- list(NULL, NULL, c(
    "Mean", "SD", "Min", "Max", "q025", "q975"
  ))
  out_array[,,"Mean"] <- pls_est_means
  out_array[,,"SD"] <- pls_est_sd
  out_array[,,"Min"] <- pls_est_min
  out_array[,,"Max"] <- pls_est_max
  out_array[,,"q025"] <- pls_est_lo
  out_array[,,"q975"] <- pls_est_hi

  saveRDS(out_array, outfile)
  invisible(outfile)
}

do_plsr <- function(rbrick, plsr_coef_file) {
  plsr_data <- read.csv(plsr_coef_file)
  pls_wl <- names(plsr_data) %>%
    str_subset("WVL_") %>%
    str_remove("WVL_") %>%
    as.numeric()
  pls_int <- plsr_data[, "Intercept"]
  pls_slope <- plsr_data %>%
    select(starts_with("WVL_")) %>%
    as.matrix() %>%
    t()
  stopifnot(
    ncol(pls_slope) == length(pls_int),
    nrow(pls_slope) == length(pls_wl)
  )
  bwl <- get_raster_wl(rbrick)

  bdat <- rbrick[]
  # Row order is b[1,1], b[1,2], b[1,3]...b[2,1], b[2,2], b[2,3]...
  # I.e., Reflectance samples are kept together

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
  bdat_r <- t(apply(bdat, 1, interp_wl)) * 10000

  # Rows = Model samples; Cols = pixels
  pls_est <- pls_int + t(bdat_r %*% pls_slope)
  pls_est_means <- colMeans(pls_est, na.rm = TRUE)
  pls_est_sd <- apply(pls_est, 2, sd, na.rm = TRUE)
  pls_est_min <- apply(pls_est, 2, min, na.rm = TRUE)
  pls_est_max <- apply(pls_est, 2, max, na.rm = TRUE)
  pls_est_lo <- apply(pls_est, 2, stats::quantile, 0.025, na.rm = TRUE)
  pls_est_hi <- apply(pls_est, 2, stats::quantile, 0.975, na.rm = TRUE)

  out_array <- array(NA_real_, c(ncol(rbrick), nrow(rbrick), 6))
  dimnames(out_array) <- list(NULL, NULL, c(
    "Mean", "SD", "Min", "Max", "q025", "q975"
  ))
  out_array[,,"Mean"] <- pls_est_means
  out_array[,,"SD"] <- pls_est_sd
  out_array[,,"Min"] <- pls_est_min
  out_array[,,"Max"] <- pls_est_max
  out_array[,,"q025"] <- pls_est_lo
  out_array[,,"q975"] <- pls_est_hi
  out_array
}

spec_plsr <- function(specdata, plsr_data) {
  pls_wl <- names(plsr_data) %>%
    str_subset("WVL_") %>%
    str_remove("WVL_") %>%
    as.numeric()

  pls_int <- plsr_data[, "Intercept"]
  pls_slope <- plsr_data %>%
    select(starts_with("WVL_")) %>%
    as.matrix() %>%
    t()

  npls <- length(pls_int)

  stopifnot(
    ncol(pls_slope) == length(pls_int),
    nrow(pls_slope) == length(pls_wl)
  )

  # Observations, samples, wavelengths
  bdat <- specdata
  bwl <- attr(specdata, "wavelength")

  # Define dimensions, to make things clearer
  nobs <- nrow(bdat)
  nsamp <- ncol(bdat)
  nwl <- length(bwl)

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
  pls_est_lo <- apply(pls_est, c(2, 3), quantile, 0.025, na.rm = TRUE)
  pls_est_hi <- apply(pls_est, c(2, 3), quantile, 0.975, na.rm = TRUE)

  out_array <- array(NA_real_, c(nobs, nsamp, 6))
  dimnames(out_array) <- list(NULL, NULL, c(
    "Mean", "SD", "Min", "Max", "q025", "q975"
  ))
  out_array[,,"Mean"] <- pls_est_means
  out_array[,,"SD"] <- pls_est_sd
  out_array[,,"Min"] <- pls_est_min
  out_array[,,"Max"] <- pls_est_max
  out_array[,,"q025"] <- pls_est_lo
  out_array[,,"q975"] <- pls_est_hi

  out_array
}
