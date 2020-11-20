do_true_plsr <- function(input_refl_file) {
  input_refl_r <- hyperSpec::read.ENVI(
    input_refl_file,
    paste0(input_refl_file, ".hdr")
  )
  input_wl <- hyperSpec::wl(input_refl_r)
  input_refl <- input_refl_r$spc
  # run PLSR models to get "true" values
  pls_coef_files <- dir_ls("data/plsr-coefficients/", recurse = TRUE,
                           regexp = "PLS.*Mode_A")
  pls_data_l <- map(pls_coef_files, read.csv)
  plsr_result <- map(pls_data_l, spec_plsr,
                     bdat = input_refl, bwl = input_wl)

  true_plsr_raw <- vctrs::vec_c(!!!plsr_result)
  dimnames(true_plsr_raw) <- dimnames(plsr_result[[1]])
  rownames(true_plsr_raw) <- pls_coef_files %>%
    path_file() %>%
    path_ext_remove() %>%
    str_remove("PLS_coefficients_Mode_A_")

  plsr_dt_long <- reshape2::melt(true_plsr_raw, varnames = c(
    "variable", "iobs", "irefl", "stat"
  )) %>% as.data.table()
  plsr_dt <- dcast(plsr_dt_long, iobs + irefl + variable ~ stat)
  plsr_dt[, `:=`(fname = "input")]
  plsr_dt
}
