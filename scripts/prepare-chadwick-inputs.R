#!/usr/bin/env Rscript

# Create all inputs for trait analysis workflow
# (1) Convert NEON image CSV samples to rasters
# (2) Create Hypertrace configuration files

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

library(vroom)
library(raster)
library(fs)
library(dplyr)
library(tidyr)

neon_wl <- readLines(path(
  "data", "chadwick", "kdchadwick-east_river_trait_modeling-5f6b260",
  "data", "neon_wavelengths.txt"
)) %>% as.numeric()

bad_wl <- function(wl) {
  (wl < 400) | (wl > 1350 & wl < 1400) | (wl > 1815 & wl < 1960) | (wl > 2460)
}

waves <- tibble(
  wavelength = neon_wl,
  band = sprintf("refl_B_%d", seq_along(wavelength))
)

crown_data <- vroom(path(
  "data", "chadwick", "kdchadwick-east_river_trait_modeling-5f6b260",
  "data", "20200214_CRBU2018_AOP_Crowns_extraction.csv"
))

crown_data_long <- crown_data %>%
  select(
    -starts_with("_obs_"),
    -starts_with("ered"),
    -starts_with("_tch_"),
    -starts_with("wtr"), -starts_with("h_me")
  ) %>%
  mutate(spec_id = row_number(), .before = 1) %>%
  pivot_longer(starts_with("refl_"), names_to = "band", values_to = "reflectance") %>%
  left_join(waves, "band") %>%
  mutate(reflectance = reflectance / 10000)

if (FALSE) {
  library(ggplot)
  ids <- crown_data_long %>% distinct(ID) %>% pull() %>% sort()
  crown_data_long %>%
    filter(ID %in% head(ids)) %>%
    mutate(reflectance = if_else(bad_wl(wavelength), NA_real_, reflectance)) %>%
    ggplot() +
    aes(x = wavelength, y = reflectance, group = spec_id, color = factor(ID)) +
    geom_line()
}

wl <- crown_data_long %>%
  distinct(wavelength) %>%
  pull()
bbl <- as.integer(bad_wl(wl))
spec_wide <- crown_data_long %>%
  select(-band) %>%
  pivot_wider(names_from = "wavelength", values_from = "reflectance")
spec_mat <- spec_wide %>%
  select(-(spec_id:Y_UTM)) %>%
  as.matrix() %>%
  `rownames<-`(spec_wide$spec_id)
spec_mat[spec_mat < 0] <- 0
spec_array <- array(NA_real_, c(1, nrow(spec_mat), ncol(spec_mat)))
spec_array[1,,] <- spec_mat
spec_raster <- brick(spec_array)

outdir <- dir_create(path("data", "derived", "ht-trait-inputs"))
outfile <- path(outdir, "chadwick.envi")
if (file_exists(outfile)) file_delete(outfile)
writeRaster(spec_raster, outfile, format = "ENVI",
            options = c("INTERLEAVE=BIL", "SUFFIX=ADD"),
            overwrite = TRUE)
write(c("wavelength={", paste(wl, collapse = ",\n"), "}"),
      paste0(outfile, ".hdr"),
      append = TRUE)
write(c("bbl={", paste(bbl, collapse = ","), "}"),
      paste0(outfile, ".hdr"),
      append = TRUE)
