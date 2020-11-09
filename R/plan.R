htpkgdir <- dir_create(path("data", "derived", "zz-traitanalysis"))

pls_coef_files <- dir_ls("data/plsr-coefficients/singh-ea-supplement",
                         regexp = "PLS_coefficients_Mode_A")
pls_tags <- str_extract(pls_coef_files, "(?<=_A_).*(?=\\.csv)")

resultdir <- path("data", "derived", "ht-output")
reflfiles <- dir_ls(resultdir, regexp = "estimated-reflectance(-[[:digit:]]+)?$", recurse = TRUE)
uncfiles <- dir_ls(resultdir, regexp = "posterior-uncertainty(-[[:digit:]]+)?$", recurse = TRUE)
radfiles <- dir_ls(resultdir, regexp = "toa-radiance(-[[:digit:]]+)?$", recurse = TRUE)

figdir <- dir_create("figures")

ispec <- 12
bad_wl <- function(wl) {
  (wl < 400) |
    (wl > 1300 & wl < 1450) |
    ## (wl > 1780 & wl < 1950) |
    (wl > 1780 & wl < 2000) |
    (wl > 2450)
}

choosefile <- . %>%
  str_subset("cal_calibration_drift") %>%
  str_subset("scale_32") %>%
  str_subset("aod_0.09")

cal_color <- RColorBrewer::brewer.pal(5, "Set1")

pls_refl_grid <- crossing(
  nesting(pls_tags, pls_coef_files),
  rasterfile = reflfiles
) %>%
  mutate(
    targetname = sprintf("plsr-%04d", row_number()),
    outfile = path(path_dir(rasterfile), paste0("ZZPLSR-", pls_tags, ".envi"))
  )

the_plan <- drake_plan(
  # Propagate uncertainty through outputs
  propagated_refl = target(
    reflectance_uncertainty(reflfiles, uncfiles, str_replace(
      reflfiles, "estimated-reflectance", "uncertain-reflectance"
    ) %>% str_c(".envi")),
    format = "file",
    dynamic = map(reflfiles, uncfiles)
  ),

  # Run PLSR on propagated reflectance
  plsr_result = target(
    raster_plsr(
      propagated_refl,
      str_replace(propagated_refl, "uncertain-reflectance", str_extract(pls_coef_files, "(?<=_A_).*(?=\\.csv)")) %>%
        str_replace("\\.envi$", ".rds"),
      pls_coef_files
    ),
    dynamic = cross(propagated_refl, pls_coef_files),
    format = "file"
  ),

  poster_true_refl = {
    true_refl_file <- file_in("data/derived/zz-traitanalysis/neon-vegetation.envi")
    true_refl_r <- brick(true_refl_file)
    true_wl <- get_raster_wl(true_refl_r)
    ibad <- bad_wl(true_wl)
    true_refl <- drop(true_refl_r[1, ispec])
    true_refl[ibad] <- NA_real_

    png(file_out(!!path(figdir, "poster-true-refl.png")),
        width = 6, height = 4, units = "in", res = 300)
    par(mar = c(3, 3, 0.1, 0.3), mgp = c(2, 0.7, 0))
    plot(true_wl, true_refl, type = "l",
         ylab = "Reflectance [0,1]",
         xlab = "Wavelength (nm)")
    dev.off()
  },

  poster_toa_radiance = {
    radfile <- radfiles %>%
      choosefile %>%
      str_subset("radiance$")
    rad_r <- brick(radfile)

    wl <- get_raster_wl(rad_r)
    rad <- drop(rad_r[1, ispec])

    png(file_out(!!path(figdir, "poster-toa-radiance.png")),
        width = 3.7, height = 2.7, units = "in", res = 300)
    par(mar = c(3, 3, 0.1, 0.3), mgp = c(1.7, 0.7, 0))
    plot(wl, rad, type = "l",
         xlab = "Wavelength (nm)",
         ylab = expression("TOA radiance" ~ (W ~ m^-2)))
    dev.off()
  },

  poster_cal_radiance = {
    radfile_l <- radfiles %>%
      choosefile %>%
      str_subset("radiance-[[:digit:]]+")
    rad_rl <- lapply(radfile_l, brick)
    wl <- get_raster_wl(rad_rl[[1]])
    rad <- sapply(rad_rl, function(x) drop(x[1, ispec]))
    png(file_out(!!path(figdir, "poster-cal-radiance.png")),
        width = 3.7, height = 2.7, units = "in", res = 300)
    par(mar = c(3, 3, 0.1, 0.3), mgp = c(1.7, 0.7, 0))
    matplot(wl, rad, type = "l", lty = "solid",
            col = cal_color,
            xlab = "Wavelength (nm)",
            ylab = expression("TOA radiance" ~ (W ~ m^-2)))
    dev.off()
    ii <- wl < 600
    png(file_out(!!path(figdir, "poster-cal-radiance-zoom.png")),
        width = 3.3, height = 2.7, units = "in", res = 300)
    par(mar = c(1.7, 1.7, 0.1, 0.3), mgp = c(1.7, 0.7, 0))
    matplot(wl[ii], rad[ii,], type = "l", lty = "solid",
            col = cal_color,
            xlab = "", ylab = "")
    dev.off()
  },

  poster_cal_reflectance = {
    reflfile_l <- reflfiles %>%
      choosefile %>%
      str_subset("estimated-reflectance-[[:digit:]]+")
    refl_rl <- lapply(reflfile_l, brick)
    wl <- get_raster_wl(refl_rl[[1]])
    ibad <- bad_wl(wl)
    refl <- sapply(refl_rl, function(x) drop(x[1, ispec]))
    refl[ibad,] <- NA_real_
    png((path(figdir, "poster-cal-reflectance.png")),
        width = 3.7, height = 2.7, units = "in", res = 300)
    par(mar = c(3, 3, 0.1, 0.3), mgp = c(1.7, 0.7, 0))
    matplot(wl, refl, type = "l", lty = "solid",
            col = cal_color,
            xlab = "Wavelength (nm)",
            ylab = "Reflectance [0,1]")
    dev.off()
    ii <- wl < 600
    png(file_out(!!path(figdir, "poster-cal-reflectance-zoom.png")),
        width = 3.3, height = 2.7, units = "in", res = 300)
    par(mar = c(1.7, 1.7, 0.1, 0.3), mgp = c(1.7, 0.7, 0))
    matplot(wl[ii], refl[ii,], type = "l", lty = "solid",
            col = cal_color,
            xlab = "", ylab = "")
    dev.off()
  },

  poster_distributions_figure = {
    png(path(figdir, "poster-distribution-conceptual.png"),
        width = 4.2, height = 1.0, units = "in", res = 300)
    par(mar = c(0, 0, 0, 0))
    curve(dnorm(x, 0, 0.2), -3, 3, n = 1000,
          axes = FALSE, xlab = "", ylab = "")
    curve(dnorm(x, 0, 0.5), -3, 3, n = 1000, add = TRUE, col = 2)
    curve(dnorm(x, 0, 1.0), -3, 3, n = 1000, add = TRUE, col = 3)
    dev.off()
  },

  poster_sensitivity_atm = {
    radfile_l <- radfiles %>%
      str_subset("cal_NONE") %>%
      str_subset("szen_0.00__ozen_0.00__saz_0.00__oaz_0.00") %>%
      str_subset("snr_300") %>%
      str_subset("toa-radiance$")
    rad_rl <- lapply(radfile_l, brick)
    wl <- get_raster_wl(rad_rl[[1]])
    rad <- sapply(rad_rl, function(x) drop(x[1, ispec]))
    atm_colors <- RColorBrewer::brewer.pal(ncol(rad), "Paired")
    png(file_out(!!path(figdir, "poster-radiance-atmospheres.png")),
        width = 3.7, height = 2.7, units = "in", res = 300)
    par(mar = c(3, 3, 0.1, 0.3), mgp = c(1.7, 0.7, 0))
    matplot(wl, rad, type = "l",
            col = atm_colors, lty = "solid",
            xlab = "Wavelength (nm)",
            ylab = expression("TOA Radiance" ~ (W ~ m^-2)))
    dev.off()
  },

  poster_sensitivity_geom = {
    radfile_l <- radfiles %>%
      str_subset("cal_NONE") %>%
      str_subset("atm_midlatitude_summer") %>%
      str_subset("aod_0.090") %>%
      str_subset("snr_300") %>%
      str_subset("toa-radiance$")
    rad_rl <- lapply(radfile_l, brick)
    wl <- get_raster_wl(rad_rl[[1]])
    rad <- sapply(rad_rl, function(x) drop(x[1, ispec]))
    atm_colors <- RColorBrewer::brewer.pal(ncol(rad), "Paired")
    png(file_out(!!path(figdir, "poster-radiance-geometries.png")),
        width = 3.7, height = 2.7, units = "in", res = 300)
    par(mar = c(3, 3, 0.1, 0.3), mgp = c(1.7, 0.7, 0))
    matplot(wl, rad, type = "l",
            col = atm_colors, lty = "solid",
            xlab = "Wavelength (nm)",
            ylab = expression("TOA Radiance" ~ (W ~ m^-2)))
    dev.off()
  },

  poster_calibration_noise = {
    noise <- R.matlab::readMat("data/derived/zz-traitanalysis/calibration_drift.mat")
    wl <- drop(noise$wavelengths)
    covar <- noise$Covariance
    draws <- rmvnorm(5, rep(0, length(wl)), covar)

    png(file_out(!!path(figdir, "poster-calibration-noise.png")),
        width = 3.4, height = 1.8, units = "in", res = 300)
    par(mar = c(2.8, 2.8, 0.1, 0.3), mgp = c(1.4, 0.4, 0),
        cex = 0.8, cex.axis = 0.8)
    matplot(wl, draws, type = "l", lty = "solid",
            col = RColorBrewer::brewer.pal(5, "Set1"),
            xlab = "Wavelength (nm)",
            ylab = "Calibration noise")
    abline(h = 0, lty = "dashed", col = "black")
    dev.off()
  }

)
