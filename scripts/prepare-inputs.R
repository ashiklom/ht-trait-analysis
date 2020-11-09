#!/usr/bin/env Rscript

# Create all inputs for trait analysis workflow
# (1) Convert NEON image CSV samples to rasters
# (2) Create Hypertrace configuration files

library(conflicted)
conflict_prefer("select", "dplyr")

library(vroom)
library(raster)
library(fs)
library(magrittr, include.only = "%>%")
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)

neon_wl <- vroom(path("data", "neon-wl.txt"), col_types = "cdd")

tidy_spec <- function(fname) {
  neon_spec <- cols("system:index" = "c", .default = "d", ".geo" = "_")
  dat <- vroom(fname, col_types = neon_spec) %>%
    rename(id = "system:index") %>%
    pivot_longer(-id, names_to = "band", values_to = "value") %>%
    left_join(neon_wl, "band") %>%
    select(-fwhm) %>%
    arrange(id, wavelength)
}

dat2raster <- function(dat, outfile) {
  wl <- dat %>%
    distinct(wavelength) %>%
    pull()
  spec_wide <- dat %>%
    select(-band) %>%
    pivot_wider(names_from = "wavelength", values_from = "value")
  spec_mat <- spec_wide %>%
    select(-id) %>%
    as.matrix() %>%
    `rownames<-`(spec_wide$id)
  spec_mat[spec_mat < 0] <- 0
  spec_array <- array(NA_real_, c(1, nrow(spec_mat), ncol(spec_mat)))
  spec_array[1,,] <- spec_mat
  spec_raster <- raster::brick(spec_array)
  raster::writeRaster(spec_raster, outfile, format = "ENVI",
                      options = c("INTERLEAVE=BIP", "SUFFIX=ADD"),
                      overwrite = TRUE)
  write(c("wavelength={", paste(wl, collapse = ",\n"), "}"),
        paste0(outfile, ".hdr"),
        append = TRUE)
  write(c("bbl={", paste(rep_len(0, length(wl)), collapse = ",\n"), "}"),
        paste0(outfile, ".hdr"),
        append = TRUE)
  invisible(outfile)
}

# Vegetation samples
veg <- tidy_spec(path("data", "neon-samples", "neon-vegetation-samples.csv"))
snow <- tidy_spec(path("data", "neon-samples", "neon-snow-samples.csv"))
ground <- tidy_spec(path("data", "neon-samples", "neon-bareground-samples.csv"))

outdir <- dir_create(path("data", "derived", "zz-traitanalysis"))
dat2raster(veg, path(outdir, "neon-vegetation.envi"))
dat2raster(snow, path(outdir, "neon-snow.envi"))
dat2raster(ground, path(outdir, "neon-ground.envi"))

# Add input files if they don't exist
file_copy(path("..", "hypertrace-data", "priors", "aviris-ng", "surface_EMIT.mat"),
          path(outdir, "surface_EMIT.mat"), overwrite = TRUE)
file_copy(path("..", "hypertrace-data", "wavelengths", "aviris-ng.txt"),
          path(outdir, "aviris-ng.txt"), overwrite = TRUE)
file_copy(path("..", "hypertrace-data", "other", "20201006_calibration_drift.mat"),
          path(outdir, "calibration_drift.mat"), overwrite = TRUE)

# Set up Hypertrace config files
atmospheres <- vroom("data/atmosphere-clusters.csv", col_types = "ddddd") %>%
  mutate(atm = if_else(tropical > 0, "tropics", "midlatitude_summer"))

atmosphere_list <- atmospheres %>%
  select(atm, aod = AOD, h2o = H2O) %>%
  pmap(list) %>%
  map(unname)

DISCOVER <- TRUE
ht_common <- list(
  wavelength_file = "./zz-traitanalysis/aviris-ng.txt",
  rtm_template_file = "./lrt-template.inp",
  lutdir = "./luts",
  outdir = "./output/traitanalysis",
  isofit = list(
    forward_model = list(
      instrument = list(SNR = 300, integrations = 1),
      surface = list(surface_category = "multicomponent_surface"),
      radiative_transfer = list(
        lut_grid = list(AOT550 = c(0.001, 0.123, 0.6),
                        H2OSTR = c(1.0, 2.5, 3.25)),
        statevector = list(
          AOT550 = list(bounds = c(0.01, 0.6),
                        scale = 0.01,
                        prior_mean = 0.05, prior_sigma = 0.2,
                        init = 0.05),
          H2OSTR = list(bounds = c(1.0, 3.25),
                        scale = 0.01,
                        prior_mean = 1.75, prior_sigma = 1.0,
                        init = 1.75)
        ),
        unknowns = list(H2O_ABSCO = 0.01),
        radiative_transfer_engines = list(vswir = list(
          engine_name = "libradtran",
          engine_base_dir = "~/projects/models/libRadtran-2.0.3",
          environment = if (DISCOVER) paste(
            "module purge",
            "module load comp/intel/19.1.0.166 mpi/impi/20.0.0.166 netcdf4/4.7.4",
            "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ashiklom/.local/lib",
            sep = "; "
          ) else paste(
            "source ~/.bashrc",
            "conda activate libradtran",
            "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib",
            sep = "; "
          ),
          lut_names = c("H2OSTR", "AOT550"),
          statevector_names = c("H2OSTR", "AOT550")
        ))
      )
    ),
    implementation = list(inversion = list(windows = rbind(
      c(400, 1300), c(1450, 1780), c(1950, 2450)
    )))
  ),
  hypertrace = list(
    calibration_uncertainty_file = c(NA, "./zz-traitanalysis/calibration_drift.mat"),
    solar_zenith = c(0, 5, 10, 20, 30, 40),
    observer_zenith = c(0, 5, 10, 20, 30, 40),
    solar_azimuth = c(0, 45, 90, 180),
    n_calibration_draws = I(5),
    calibration_scale = c(2, 4, 8, 16, 32),
    atm_aod_h2o = atmosphere_list,
    inversion_mode = I("inversion"),
    snr = c(200, 300, 400, 500),
    surface_file = I("./zz-traitanalysis/surface_EMIT.mat")
  )
)

quickrun <- list(hypertrace = list(
  solar_zenith = I(0),
  observer_zenith = I(0),
  solar_azimuth = I(0),
  atm_aod_h2o = atmosphere_list[1],
  snr = I(300),
  calibration_scale = c(8, 32)
))

outdir_remote <- "./zz-traitanalysis"
veg_config <- ht_common %>%
  modifyList(list(
    reflectance_file = path(outdir_remote, "neon-vegetation.envi")
  ))
write(toJSON(veg_config, auto_unbox = TRUE, pretty = TRUE),
      path(outdir, "neon-full-vegetation.json"))
veg_config_quick <- modifyList(veg_config, quickrun)
write(toJSON(veg_config_quick, auto_unbox = TRUE, pretty = TRUE),
      path(outdir, "neon-quick-vegetation.json"))

snow_config <- ht_common %>%
  modifyList(list(
    reflectance_file = path(outdir_remote, "neon-snow.envi")
  ))
write(toJSON(veg_config, auto_unbox = TRUE, pretty = TRUE),
      path(outdir, "neon-full-snow.json"))
snow_config_quick <- modifyList(snow_config, quickrun)
write(toJSON(snow_config_quick, auto_unbox = TRUE, pretty = TRUE),
      path(outdir, "neon-quick-snow.json"))

ground_config <- ht_common %>%
  modifyList(list(
    reflectance_file = path(outdir_remote, "neon-ground.envi")
  ))
write(toJSON(veg_config, auto_unbox = TRUE, pretty = TRUE),
      path(outdir, "neon-full-ground.json"))
ground_config_quick <- modifyList(ground_config, quickrun)
write(toJSON(ground_config_quick, auto_unbox = TRUE, pretty = TRUE),
      path(outdir, "neon-quick-ground.json"))
