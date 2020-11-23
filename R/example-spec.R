example_spec <- function(iobs) {
  reflfiles <- dir_ls(resultdir, recurse = TRUE,
                      regexp = "propagated-uncertainty(-[[:digit:]]+)?\\.qs")
  refl_data_all <- parse_ht_directory(reflfiles)

  in1 <- qs::qread(reflfiles[1])
  wl <- attr(in1, "wavelength")
  refl_data <- refl_data_all %>%
    mutate(raw_spectra = map(fname, safely(read_runc), i = iobs))

  refl_data_proc <- refl_data %>%
    mutate(all_spectra = map(raw_spectra, "result")) %>%
    select(-raw_spectra) %>%
    filter(map_lgl(all_spectra, negate(is.null)))
  refl_data_proc
}

read_spec <- function(f, i) {
  dat <- hyperSpec::read.ENVI(f, paste0(f, ".hdr"))
  drop(dat$spc[i,])
}

read_runc <- function(f, i) {
  drop(qs::qread(f)[i,,])
}
