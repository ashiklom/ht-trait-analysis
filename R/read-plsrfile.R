read_plsr_file <- function(plsrfile) {
  # Variable, Obs, Refl, Stat
  plsr_raw <- qs::qread(plsrfile)
  rownames(plsr_raw) <- rownames(plsr_raw) %>%
    path_file() %>%
    path_ext_remove() %>%
    str_remove("PLS_coefficients_Mode_A_")

  plsr_dt_long <- reshape2::melt(plsr_raw, varnames = c(
    "variable", "iobs", "irefl", "stat"
  )) %>% as.data.table()

  plsr_dt <- dcast(plsr_dt_long, iobs + irefl + variable ~ stat)
  plsr_dt[, fname := plsrfile]
  plsr_dt
}
