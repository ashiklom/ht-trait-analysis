library(conflicted)
conflict_prefer("Map", "rgee")

library(fs)
library(sf)
library(rgee)
ee_Initialize()

neonrefl <- ee$ImageCollection("users/pgbrodrick/SFA/collections/ciacorn_priority")
mosaic <- neonrefl$mosaic()$divide(10000)

rgb <- mosaic$select(paste0("B", c(60, 26, 12)))

roi <- ee$Geometry$Polygon(list(c(-107.165, 39.059),
                                c(-107.165, 38.798),
                                c(-106.817, 38.798),
                                c(-106.817, 39.059)))

samps <- read.csv("point-samples.csv") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

samps_ee <- sf_as_ee(samps)

bbox <- samps_ee$geometry()$bounds()

Map$centerObject(bbox, zoom = 14)
m1 <- Map$addLayer(rgb, list(), "mosaic")
m2 <- Map$addLayer(samps_ee, list(color = "red"), "samples")
m1 + m2

# Extract pixels from image
message("Extracting pixels...")
rawdat <- ee_extract(mosaic, samps_ee, scale = 1)

library(dplyr)
library(tidyr)
avng_bands <- read.table(path("data", "aviris-ng.txt")) %>%
  select(wavelength = V2, fwhm = V3) %>%
  arrange(wavelength) %>%
  mutate(wavelength = wavelength * 1000,
         fwhm = fwhm * 1000,
         band = paste0("B", row_number() - 1)) %>%
  as_tibble()

procdat <- rawdat %>%
  mutate(id = row_number()) %>%
  pivot_longer(-id, names_to = "band", values_to = "value") %>%
  left_join(avng_bands, "band") %>%
  select(id, band, wavelength, value) %>%
  arrange(id, wavelength)

outdir <- dir_create(path("data", "derived"))
write.table(procdat, path(outdir, "neon-pixels.csv"),
            row.names = FALSE, sep = ",", quote = FALSE)
