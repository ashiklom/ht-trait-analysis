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
  },

  poster_reflectance_uncertain = {
    reflfile <- propagated_refl %>%
      str_subset("cal_NONE") %>%
      head(1)
    refl_r <- brick(reflfile)
    wl <- get_raster_wl(refl_r)
    ibad <- bad_wl(wl)
    refl <- t(refl_r[ispec, ])
    refl[refl < 0] <- 0
    refl[ibad,] <- NA
    refl_lo <- apply(refl, 1, min)
    refl_hi <- apply(refl, 1, max)
    refl_mu <- rowMeans(refl)
    plt <- ggplot() +
      aes(x = wl, y = refl_mu, ymin = refl_lo, ymax = refl_hi) +
      geom_ribbon(fill = "gray70") +
      geom_line() +
      theme_base() +
      labs(x = "Wavelength (nm)", y = "Reflectance [0, 1]")
    ggsave(file_out(!!path(figdir, "poster-reflectance-uncertainty.png")), plt,
           width = 3.7, height = 2.7, units = "in", dpi = 300)
  },

  poster_plsr_coefficients = {
    pls_coefs_raw <- map_dfr(pls_coef_files, read.csv, .id = "file") %>%
      mutate(file = path_file(file) %>%
               str_remove("PLS_coefficients_Mode_A_") %>%
               str_remove("\\.csv")) %>%
      as_tibble()
    pls_coefs_tidy <- pls_coefs_raw %>%
      select(-Intercept) %>%
      pivot_longer(
        starts_with("WVL_"),
        names_to = "wavelength",
        values_to = "value",
        names_transform = list(wavelength = ~as.numeric(str_remove(.x, "WVL_")))
      )

    pls_coefs_summary <- pls_coefs_tidy %>%
      group_by(file, wavelength) %>%
      summarize(
        Mean = mean(value),
        lo = quantile(value, 0.025),
        hi = quantile(value, 0.975)
      ) %>%
      ungroup()

    plt <- pls_coefs_summary %>%
      mutate(variable = factor(
        file,
        c("pct_N", "pct_C", "M_area"),
        c("Leaf %N", "Leaf %C", "Leaf mass per area")
      )) %>%
      filter(!is.na(variable)) %>%
      ggplot() +
      aes(x = wavelength, y = Mean, ymin = lo, ymax = hi) +
      geom_ribbon(fill = "gray80") +
      geom_line() +
      facet_grid(vars(variable), scales = "free_y",
                 labeller = label_wrap_gen(10)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_base() +
      theme(axis.text.y = element_blank(),
            axis.ticks = element_blank()) +
      labs(x = "Wavelength (nm)", y = "PLS coefficient")

    ggsave(file_out(!!path(figdir, "poster-PLS-coefficients.png")), plt,
           width = 3.3, height = 4.0, units = "in", dpi = 600)
  },

  poster_plsr_one_result = {
    plsrfiles <- plsr_result %>%
      str_subset("(M_area|pct_C|pct_N)\\.rds") %>%
      head(3)
    dq <- distributions3::quantile
    plsr_l <- lapply(plsrfiles, readRDS)
    dists <- map(plsr_l, ~Normal(.x[1, ispec, "Mean"], .x[1, ispec, "SD"]))
    names(dists) <- path_file(path_ext_remove(plsrfiles))
    png(path(figdir, "poster-plsr-one-result.png"), width = 2, height = 4,
        units = "in", res = 600)
    par(mfrow = c(3, 1), mar = c(3, 0.1, 0.1, 0.1))
    dn <- dists[["pct_N"]]
    curve(pdf(dn, x), dq(dn, 0.001), dq(dn, 0.999),
          n = 1000, xlab = "", ylab = "", yaxt = "n")
    dc <- dists[["pct_C"]]
    curve(pdf(dc, x), dq(dc, 0.001), dq(dc, 0.999),
          n = 1000, xlab = "", ylab = "", yaxt = "n")
    dl <- dists[["M_area"]]
    curve(pdf(dl, x), dq(dl, 0.001), dq(dl, 0.999),
          n = 1000, xlab = "", ylab = "", yaxt = "n")
    dev.off()
  },

  poster_plsr_results = target({
    all_plsrfiles <- dir_ls("data/derived/ht-output",
                            regexp = "\\.rds", recurse = TRUE)
    plsr_info <- parse_ht_directory(all_plsrfiles) %>%
      mutate(plsr_var = path_ext_remove(fname) %>%
               path_file() %>%
               str_remove("-[[:digit:]]+"))

    tidy_plsr <- function(fname) {
      dat <- readRDS(fname)
      rownames(dat) <- sprintf("r%02d", seq_len(nrow(dat)))
      colnames(dat) <- sprintf("i%02d", seq_len(ncol(dat)))
      datm <- apply(dat, 3, c)
      meta <- expand.grid(dimnames(dat)[1:2])
      result <- data.frame(meta, datm) %>%
        as_tibble() %>%
        mutate(irefl = as.numeric(str_remove(Var1, "r")),
               iobs = as.numeric(str_remove(Var2, "i"))) %>%
        select(iobs, irefl, everything()) %>%
        select(-Var1, -Var2)
      result
    }

    all_plsr <- plsr_info %>%
      mutate(data = map(fname, tidy_plsr)) %>%
      unnest(data) %>%
      as.data.table()

    all_plsr
  }, formats = "fst_dt"),

  poster_plsr_summary_dat = {

    setDT(poster_plsr_results)

    var_calibration <- poster_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
      plsr_var, atm, aod, h2o, szen, ozen, saz, oaz, snr,
      calibration, calscale,
      iobs, irefl
    )][
    , .(SD = mean(SD, na.rm = TRUE)),
      .(plsr_var)
    ][, source := "calibration"]

    # Variance across atmospheric correction -- not irefl
    var_isofit <- poster_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
      plsr_var, atm, aod, h2o, szen, ozen, saz, oaz, snr,
      calibration, calscale, calid,
      iobs
    )][
    , .(SD = mean(SD, na.rm = TRUE)),
      .(plsr_var)
    ][, source := "isofit"]

    # Variance by atmosphere -- not atm, aod, h2o
    var_atm <- poster_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
      plsr_var, szen, ozen, saz, oaz, snr,
      calibration, calscale, calid,
      iobs, irefl
    )][
    , .(SD = mean(SD, na.rm = TRUE)),
      .(plsr_var)
    ][, source := "atmosphere"]

    # Variance by SNR
    var_snr <- poster_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
      plsr_var, atm, aod, h2o, szen, ozen, saz, oaz,
      calibration, calscale, calid,
      iobs, irefl
    )][
    , .(SD = mean(SD, na.rm = TRUE)),
      .(plsr_var)
    ][, source := "SNR"]

    # Variance by geometry  -- not szen, ozen, saz, oaz
    var_geom <- poster_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
      plsr_var, atm, aod, h2o, snr,
      calibration, calscale, calid,
      iobs, irefl
    )][
    , .(SD = mean(SD, na.rm = TRUE)),
      .(plsr_var)
    ][, source := "geometry"]

    # Variance by algorithm
    var_algorithm <- poster_plsr_results[, .(SD = mean(SD, na.rm = TRUE)), .(
      plsr_var
    )][, source := "algorithm"]

    var_summary <- rbindlist(list(
      var_calibration, var_isofit, var_atm, var_snr,
      var_geom, var_algorithm, mutate(true_plsr_var, source = "input")
    ))

    var_summary
  },
  poster_plsr_summary_fig = {
    summary_dat <- poster_plsr_summary_dat %>%
      mutate(
        source = factor(
          source,
          c("input", "calibration", "isofit", "algorithm",
            "atmosphere", "SNR", "geometry"),
          c("Observations", "Calibration", "Atm. correction", "Algorithm",
            "Atmosphere", "Signal-noise ratio", "Sun-sensor geometry")
        ) %>%
          fct_reorder(SD, .desc = TRUE, .fun = mean) %>%
          fct_relevel("Observations"),
        plsr_var = factor_pls(plsr_var)) %>%
      filter(!is.na(plsr_var))

    plt <- ggplot(summary_dat) +
      aes(x = source, y = SD) +
      geom_bar(stat = "identity") +
      facet_grid(vars(plsr_var), scales = "free") +
      labs(x = "Variance source", y = "Average standard deviation") +
      scale_x_discrete(labels = label_wrap_gen(12)) +
      theme_base() +
      theme(axis.text.x = element_text(size = rel(0.9)))

    ggsave(file_out(!!path(figdir, "poster-variance-barplot.png")), plt,
           width = 8, height = 6, units = "in", dpi = 450)
  },

  true_plsr = {
    # load input reflectance
    input_refl_file <- file_in("data/derived/zz-traitanalysis/neon-vegetation.envi")
    input_refl_r <- brick(input_refl_file)
    wl <- get_raster_wl(input_refl_r)
    input_refl <- input_refl_r[]
    # run PLSR models to get "true" values
    plsr_result <- map(pls_coef_files, do_plsr, rbrick = input_refl_r)
    true_plsr <- plsr_result %>%
      map(drop) %>%
      map(as_tibble) %>%
      bind_rows(.id = "fname") %>%
      group_by(fname) %>%
      mutate(iobs = row_number()) %>%
      ungroup() %>%
      mutate(plsr_var = path_file(fname) %>%
               path_ext_remove() %>%
               str_remove("PLS_coefficients_Mode_A_"))
    true_plsr
  },

  # get true observed variance; replace this in above figure
  true_plsr_var = {
    true_plsr %>%
      group_by(plsr_var) %>%
      summarize(SD = sd(Mean, na.rm = TRUE)) %>%
      ungroup()
  },

  # calculate biases
  plsr_bias = target({
    true_plsr_dt <- as.data.table(true_plsr)
    setkey(true_plsr_dt, plsr_var, iobs)
    setkey(poster_plsr_results, plsr_var, iobs)
    biases <- true_plsr_dt[poster_plsr_results]
    biases[, bias := i.Mean - Mean]
    biases
  }, format = "fst_dt"),

  # plot biases
  poster_plsr_bias_byatm = {
    bias_byatm <- plsr_bias[, .(
      mean_bias = mean(bias, na.rm = TRUE),
      lo_bias = quantile(bias, 0.025, na.rm = TRUE),
      hi_bias = quantile(bias, 0.975, na.rm = TRUE)
    ), .(plsr_var, atm, aod, h2o)]

    plt <- bias_byatm %>%
      arrange(atm, aod, h2o) %>%
      mutate(plsr_var = factor_pls(plsr_var),
             atmosphere = label_atmosphere(atm, aod, h2o) %>%
               fct_inorder()) %>%
      filter(!is.na(plsr_var)) %>%
      ggplot() +
      aes(x = aod, y = mean_bias,
          ymin = lo_bias, ymax = hi_bias, color = atm) +
      geom_pointrange() +
      geom_hline(yintercept = 0, lty = "dashed") +
      facet_grid(vars(plsr_var), scales = "free_y") +
      scale_x_log10() +
      labs(x = "Aerosol optical depth", y = "Bias (mean and 95% CI)",
           color = "Atm. profile type") +
      theme_base() +
      theme(legend.position = c(0, 1),
            legend.justification = c(0, 1),
            legend.background = element_blank())
    ggsave(file_out(!!path(figdir, "poster-plsr-bias-byatm.png")), plt,
           width = 7, height = 6, units = "in", dpi = 450)
  }

)

if (FALSE) {


    ggplot(var_summary) +
      aes(x = source, y = SD) +
      geom_bar(stat = "identity") +
      facet_wrap(vars(plsr_var), scales = "free_y") +
      theme_bw()

    var_summary %>%
      group_by(plsr_var) %>%
      mutate(x = -rank(SD)) %>%
      ggplot() +
      aes(x = x, y = ) +
      scale_x_discrete()

    input_var <- all_plsr %>%
      group_by(plsr_var, iobs) %>%
      summarize(Mean = mean(Mean, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(variance_source = "input")

    # Variance from calibration
    cal_var <- all_plsr %>%
      filter(calibration != "NONE") %>%
      group_by(plsr_var, calscale, calid) %>%
      summarize(Mean = mean(Mean, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(variance_source = "calibration")

    # Variance from atmospheric correction
    isofit_var <- all_plsr %>%
      group_by(plsr_var, irefl) %>%
      summarize(Mean = mean(Mean, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(variance_source = "atmospheric correction")

    cal_var %>%
      filter(calscale == 32) %>%
      bind_rows(input_var) %>%
      ggplot() +
      aes(x = variance_source, y = Mean) +
      geom_boxplot() +
      facet_wrap(vars(plsr_var), scales = "free_y")


    bias_

}
