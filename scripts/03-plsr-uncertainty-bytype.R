library(data.table)
library(targets)
stopifnot(requireNamespace("fst", quietly = TRUE))

tar_load("plsr_vars")

resultlist <- list()

# Targets can't intelligently handle the memory load of this, so we chunk the
# execution by hand.

fstfile <- "_targets/objects/raw_plsr_results"
for (plsr_var in plsr_vars) {
  message(plsr_var)
  rawdat <- fst::fst(fstfile)
  raw_plsr_results <- as.data.table(rawdat[rawdat$variable == plsr_var, c(
    "fname", "iobs", "irefl", "Mean", "SD", "q025", "q975"
  )])
  plsr_meta <- unique(raw_plsr_results[["fname"]]) %>%
    parse_ht_directory() %>%
    as.data.table()
  setkey(plsr_meta, fname)
  setkey(raw_plsr_results, fname)
  raw_plsr_results <- plsr_meta[raw_plsr_results]
  fst::write_fst(raw_plsr_results, file.path("data", "derived", paste0(
    "plsr-", plsr_var, ".fst"
  )))

  var_calibration <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    atm, aod, h2o, szen, ozen, saz, oaz, snr,
    calibration, calscale,
    iobs, irefl
  )][, .(SD = mean(SD, na.rm = TRUE))][, source := "calibration"]

  # Variance across atmospheric correction -- not irefl
  var_isofit <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    atm, aod, h2o, szen, ozen, saz, oaz, snr,
    calibration, calscale, calid,
    iobs
  )][, .(SD = mean(SD, na.rm = TRUE))][, source := "isofit"]

  # Variance by atmosphere -- not atm, aod, h2o
  var_atm <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    szen, ozen, saz, oaz, snr,
    calibration, calscale, calid,
    iobs, irefl
  )][, .(SD = mean(SD, na.rm = TRUE))][, source := "atmosphere"]

  # Variance by SNR
  var_snr <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    atm, aod, h2o, szen, ozen, saz, oaz,
    calibration, calscale, calid,
    iobs, irefl
  )][, .(SD = mean(SD, na.rm = TRUE))][, source := "SNR"]

  # Variance by geometry  -- not szen, ozen, saz, oaz
  var_geom <- raw_plsr_results[, .(SD = sd(Mean, na.rm = TRUE)), .(
    atm, aod, h2o, snr,
    calibration, calscale, calid,
    iobs, irefl
  )][, .(SD = mean(SD, na.rm = TRUE))][, source := "geometry"]

  # Variance by algorithm
  var_algorithm <- raw_plsr_results[, .(SD = mean(SD, na.rm = TRUE))][, source := "algorithm"]

  var_all <- rbindlist(list(
    var_calibration, var_isofit, var_atm,
    var_snr, var_geom, var_algorithm
  ))

  resultlist[[plsr_var]] <- var_all
}

results <- rbindlist(resultlist, idcol = "variable")
qs::qsave(results, "data/derived/plsr-variance.qs")
