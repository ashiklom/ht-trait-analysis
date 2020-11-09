tidy_neon_csv <- function(dat, wavelength_file = "data/aviris-ng.txt") {
  waves <- read_table(
    wavelength_file,
    col_names = FALSE,
    col_types = "ddd"
  ) %>%
    transmute(wavelength = X2 * 1000, fwhm = X3 * 1000) %>%
    arrange(wavelength) %>%
    mutate(band = paste0("B", row_number() - 1))

  dat %>%
    rename(id = "system:index") %>%
    select(-.geo) %>%
    pivot_longer(-id, "band", "value") %>%
    left_join(waves, "band") %>%
    select(id, band, wavelength, value) %>%
    arrange(id, wavelength)
}
