tidy_example_spec <- function(example_spec_1) {
  wl <- read.table("data/derived/zz-traitanalysis/aviris-ng.txt") %>%
    pull(V2) %>%
    `*`(1000)
  dat2 <- example_spec_1 %>%
    mutate(
      irefl = list(1:20),
      all_spectra = map(all_spectra, as_tibble, .name_repair = ~as.character(wl))
    ) %>%
    unnest(c(irefl, all_spectra))
  dat3 <- dat2 %>%
    pivot_longer(
      matches("^[[:digit:].]+$"),
      names_to = "wavelength",
      values_to = "value",
      names_transform = list(wavelength = as.numeric)
    ) %>%
    as.data.table()
  dat3
}
