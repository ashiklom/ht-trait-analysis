aod_uncertainty_example <- function(example_spec_dt, input_refl) {

  specdat <- example_spec_dt
  specdat[(wavelength > 1300 & wavelength < 1450), value := NA_real_]
  specdat[(wavelength > 1780 & wavelength < 2000), value := NA_real_]

  plotdat <- specdat[irefl %in% 10:11, .(value = mean(value)), .(
    fname, atm, aod, h2o, szen, ozen, saz, oaz, snr, calibration,
    calscale, calid, wavelength
  )]

  input_sub <- input_refl[iobs == 1]
  input_sub[wl > 1300 & wl < 1450, value := NA_real_]
  input_sub[wl > 1780 & wl < 2000, value := NA_real_]

  plt <- ggplot(plotdat) +
    aes(x = wavelength, y = value, group = fname, color = aod) +
    geom_line() +
    geom_line(aes(x = wl, y = value, linetype = "Input"), data = input_sub, inherit.aes = FALSE,
              color = "red", size = 1) +
    scale_color_viridis_c() +
    theme_bw() +
    labs(x = "Wavelength (nm)", y = "Reflectance [0, 1]",
         color = "AOD", linetype = "")

  figfile <- path(figdir, "aod-reflectance.png")
  ggsave(figfile, plt, width = 5, height = 3, units = "in", dpi = 300)
  invisible(figfile)
}

aod_uncertainty_plsr_example_data <- function(raw_plsr_results, true_plsr) {
  plsr_sub <- raw_plsr_results[
    iobs == 1 & irefl %in% 10:11,
    lapply(.SD, mean),
    .(fname, variable)
  ]

  plsr_meta <- unique(plsr_sub[["fname"]]) %>%
    parse_ht_directory() %>%
    as.data.table()
  setkey(plsr_meta, fname)
  setkey(plsr_sub, fname)
  plsr_sub2 <- plsr_meta[plsr_sub]

  plsr_n <- plsr_sub2[variable == "pct_N"]

  true_plsr_sub <- true_plsr[iobs == 1 & variable == "pct_N"]

  list(plsr_n = plsr_n, true_plsr_sub = true_plsr_sub)
}

aod_uncertainty_plsr_example_figure <- function(datalist) {
  plsr_n <- datalist[["plsr_n"]]
  true_plsr_sub <- datalist[["true_plsr_sub"]]

  plt <- ggplot(plsr_n) +
    aes(x = aod, y = Mean, color = aod) +
    geom_jitter() +
    geom_hline(yintercept = true_plsr_sub$Mean,
               linetype = "dashed", color = "red") +
    scale_color_viridis_c() +
    scale_x_sqrt() +
    labs(x = "AOD", y = "Leaf %N estimate", color = "AOD") +
    theme_bw() +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.background = element_blank())

  figfile <- path(figdir, "aod-plsr-example.png")
  ggsave(figfile, plt, width = 4, height = 4, units = "in", dpi = 300)

  plt2 <- ggplot(plsr_n) +
    aes(x = aod, y = Mean, color = aod, ymin = q025, ymax = q975) +
    geom_pointrange(position = position_jitter()) +
    geom_hline(yintercept = true_plsr_sub$Mean,
               linetype = "dashed", color = "red") +
    scale_color_viridis_c() +
    scale_x_sqrt() +
    labs(x = "AOD", y = "Leaf %N estimate", color = "AOD") +
    theme_bw() +
    guides(color = FALSE)
    ## theme(legend.position = c(0, 1),
    ##       legend.justification = c(0, 1),
    ##       legend.background = element_blank())

  figfile2 <- path(figdir, "aod-plsr-example-errors.png")
  ggsave(figfile2, plt2, width = 4, height = 4, units = "in", dpi = 300)

  invisible(figfile)
}



if (FALSE) {
  tar_load(raw_plsr_results)
  tar_load(true_plsr)
  tar_load("input_refl")
  tar_load(c("true_plsr", "example_spec_1"))
  tar_load(example_spec_dt)
}
