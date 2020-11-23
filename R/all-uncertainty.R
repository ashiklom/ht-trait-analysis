plsr_results_wide_withmeta <- function(raw_plsr_results) {
  raw_plsr_results[, `:=`(
    Min = NULL,
    Max = NULL,
    q025 = NULL,
    q975 = NULL)]
  raw_plsr_results <- raw_plsr_results[variable == "M_area"]
  gc()
  plsr_meta <- unique(raw_plsr_results[["fname"]]) %>%
    parse_ht_directory() %>%
    as.data.table()
  setkey(plsr_meta, fname)
  setkey(raw_plsr_results, fname)
  raw_plsr_results <- plsr_meta[raw_plsr_results]
}

plsr_results_withmeta <- function(raw_plsr_results) {
  raw_plsr_results[, `:=`(
    Min = NULL,
    Max = NULL,
    q025 = NULL,
    q975 = NULL)]
  raw_plsr_results[, atm := str_extract(fname, "(?<=atm_).*?(?=__)")]
  raw_plsr_results[, aod = as.numeric(str_extract(fname, "(?<=aod_)[[:digit:].]+(?=__)"))]
  raw_plsr_results[, h2o = as.numeric(str_extract(fname, "(?<=h2o_)[[:digit:].]+(?=__)"))]

  ## raw_plsr_results[, `:=`(
  ##   ## Min = NULL,
  ##   ## Max = NULL,
  ##   ## q025 = NULL,
  ##   ## q975 = NULL,
  ##   ## atm = str_extract(fname, "(?<=atm_).*?(?=__)"),
  ##   ## aod = as.numeric(str_extract(fname, "(?<=aod_)[[:digit:].]+(?=__)")),
  ##   ## h2o = as.numeric(str_extract(fname, "(?<=h2o_)[[:digit:].]+(?=/)")),
  ##   szen = as.numeric(str_extract(fname, "(?<=szen_)[[:digit:].]+(?=__)")),
  ##   ozen = as.numeric(str_extract(fname, "(?<=ozen_)[[:digit:].]+(?=__)")),
  ##   saz = as.numeric(str_extract(fname, "(?<=saz_)[[:digit:].]+(?=__)")),
  ##   oaz = as.numeric(str_extract(fname, "(?<=oaz_)[[:digit:].]+(?=/)")),
  ##   snr = as.numeric(str_extract(fname, "(?<=snr_)[[:digit:].]+(?=/)")),
  ##   calibration = str_extract(fname, "(?<=cal_).*?(?=__)"),
  ##   calscale = as.numeric(str_extract(fname, "(?<=scale_)[[:digit:].]+(?=/)")),
  ##   calid = as.numeric(str_extract(path_ext_remove(fname), "(?<=-)[[:digit:]]{2}$")),
  ##   fname = NULL
  ## )]
  raw_plsr_results[, fid := .GRP, fname]
  plsr_files <- raw_plsr_results[, .N, .(fid, fname)]
  ## raw_plsr_results[, `:=`(Min = NULL, Max = NULL, q025 = NULL, q975 = NULL,
  ##                         fname = NULL)]
  ## gc()
  plsr_meta <- parse_ht_directory(raw_plsr_results[["fname"]]) %>%
    as.data.table()

  ## setkey(plsr_meta, fid)
  ## setkey(raw_plsr_results, fid)

  ## raw_plsr_results <- plsr_meta[raw_plsr_results]
  raw_plsr_results
}

plsr_uncertainty_bytype <- function(raw_plsr_results, plsr_var) {

  raw_plsr_results[, `:=`(
    Min = NULL,
    Max = NULL,
    q025 = NULL,
    q975 = NULL)]
  raw_plsr_results <- raw_plsr_results[variable == plsr_var]
  gc()
  plsr_meta <- unique(raw_plsr_results[["fname"]]) %>%
    parse_ht_directory() %>%
    as.data.table()
  setkey(plsr_meta, fname)
  setkey(raw_plsr_results, fname)
  raw_plsr_results <- plsr_meta[raw_plsr_results]

  var_calibration <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    variable, atm, aod, h2o, szen, ozen, saz, oaz, snr,
    calibration, calscale,
    iobs, irefl
  )][
  , .(SD = mean(SD, na.rm = TRUE)),
    .(plsr_var)
  ][, source := "calibration"]

  # Variance across atmospheric correction -- not irefl
  var_isofit <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    variable, atm, aod, h2o, szen, ozen, saz, oaz, snr,
    calibration, calscale, calid,
    iobs
  )][
  , .(SD = mean(SD, na.rm = TRUE)),
    .(plsr_var)
  ][, source := "isofit"]

  # Variance by atmosphere -- not atm, aod, h2o
  var_atm <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    variable, szen, ozen, saz, oaz, snr,
    calibration, calscale, calid,
    iobs, irefl
  )][
  , .(SD = mean(SD, na.rm = TRUE)),
    .(plsr_var)
  ][, source := "atmosphere"]

  # Variance by SNR
  var_snr <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    variable, atm, aod, h2o, szen, ozen, saz, oaz,
    calibration, calscale, calid,
    iobs, irefl
  )][
  , .(SD = mean(SD, na.rm = TRUE)),
    .(plsr_var)
  ][, source := "SNR"]

  # Variance by geometry  -- not szen, ozen, saz, oaz
  var_geom <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    variable, atm, aod, h2o, snr,
    calibration, calscale, calid,
    iobs, irefl
  )][
  , .(SD = mean(SD, na.rm = TRUE)),
    .(plsr_var)
  ][, source := "geometry"]

  # Variance by algorithm
  var_algorithm <- raw_plsr_results[, .(SD = mean(SD, na.rm = TRUE)), .(
    variable
  )][, source := "algorithm"]

  var_all <- rbindlist(list(
    var_calibration, var_isofit, var_atm,
    var_snr, var_geom, var_algorithm
  ))

  var_all

}

plsr_uncertainty_figure <- function(plsr_uncertainties, true_plsr) {
  input_var <- true_plsr[, .(SD = sd(Mean, na.rm = TRUE)), variable][
    , source := "input"
  ]

  plotdat <- rbind(plsr_uncertainties, input_var)
  plotdat[, relsd := SD / sum(SD), variable]
  plotdat[, source := fct_reorder(source, relsd, .desc = TRUE)]
  plotdat[, source := fct_relevel(source, "input")]
  # TODO: Better x-axis labels (source)
  plotdat[, variable := factor(variable, c(
    "M_area", "pct_N", "pct_C", "pct_Cellulose",
    "pct_ADL", "pct_ADF", "del15N"
  ), c(
    "LMA", "%N", "%C", "%Cellulose",
    "%ADL", "%ADF", "d15N"
  ))]
  plt <- ggplot(plotdat) +
    aes(x = source, y = SD) +
  # TODO: Input variance as horizontal line
    geom_col() +
    facet_wrap(vars(variable), scales = "free_y") +
    theme_bw() +
    labs(x = "Uncertainty source", y = "Standard deviation")

  figfile <- path(figdir, "all-plsr-uncertainty-barplot.png")
  ggsave(figfile, plt, width = 7, height = 7, units = "in", dpi = 300)
  invisible(figfile)
}
