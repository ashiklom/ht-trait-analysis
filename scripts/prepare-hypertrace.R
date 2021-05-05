library(conflicted)
conflict_prefer("select", "dplyr")

library(vroom)
library(dplyr)
library(tibble)
library(fs)
library(purrr)
library(jsonlite)

outdir <- path("data", "derived", "ht-trait-inputs")

reflectance_files <- path(".", "ht-trait-inputs", c(
  "chadwick.envi"
))

noise_files <- path(".", "ht-trait-inputs", sprintf("noise-%d.txt", c(200, 300, 400)))

atmospheres <- vroom("data/atmosphere-clusters.csv", col_types = "ddddd") %>%
  mutate(atm = if_else(tropical > 0, "ATM_TROPICAL", "ATM_MIDLAT_SUMMER"))

atmosphere_list <- atmospheres %>%
  select(atm, aod = AOD, h2o = H2O) %>%
  pmap(list) %>%
  map(unname)

ht_common <- list(
  wavelength_file = "./hypertrace-data/wavelengths/aviris-ng.txt",
  rtm_template_file = "/dev/null",
  lutdir = "./luts/sRTMnet-ht-trait",
  outdir = "./output/ht-trait-analysis",
  isofit = list(
    forward_model = list(
      instrument = list(SNR = 300, integrations = 1),
      surface = list(surface_category = "multicomponent_surface"),
      radiative_transfer = list(
        lut_grid = list(AOT550 = c(0.01, 0.1, 0.258, 0.505, 0.752, 1.0),
                        H2OSTR = c(0.1, 0.6125, 1.325, 2.55, 3.775, 5.0)),
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
          engine_name = "simulated_modtran",
          engine_base_dir = "./6sv-2.1",
          aerosol_template_path = "../../data/aerosol_template.json",
          earth_sun_distance_file = "../../data/earth_sun_distance.txt",
          irradiance_file = "../../examples/20151026_SantaMonica/data/prism_optimized_irr.dat",
          emulator_aux_file = "./sRTMnet_v100/sRTMnet_v100_aux.npz",
          emulator_file = "./sRTMnet_v100/sRTMnet_v100",
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
    longitude = I(0),
    latitude = I(0),
    elevation_km = c(0.01, 0.5, 1, 2),
    atm_aod_h2o = atmosphere_list,
    inversion_mode = I("inversion"),
    surface_file = I("./hypertrace-data/priors/aviris-ng/surface_EMIT.mat"),
    outdir_scheme = I("hash")
  )
)
observer_zenith <- c(0, 10, 30)
observer_azimuth <- c(0, 90, 180)
localtime <- c(8.5, 9.5, 10, 10.5, 12, 15)

ht_iterations <- prod(lengths(ht_common$hypertrace))
ht_loops <- prod(length(noise_files), length(observer_zenith), length(observer_azimuth), length(localtime))
message("Total iterations per job: ", ht_iterations)
message("Total jobs: ", ht_loops)
message("Input reflectances: ", length(reflectance_files))
message("Total output size: ", ht_iterations * ht_loops * length(reflectance_files))

for (rfile in reflectance_files) {
  i <- 0
  for (noisefile in noise_files) {
    for (ozen in observer_zenith) {
      for (oaz in observer_azimuth) {
        for (lt in localtime) {
          i <- i + 1
          outfile <- path(outdir, "configs", sprintf("%s-%03d.json", path_ext_remove(basename(rfile)), i))
          out_nc <- path(ht_common[["outdir"]], sprintf("%s-%03d.nc", path_ext_remove(basename(rfile)), i))
          ht_settings <- modifyList(ht_common, list(
            reflectance_file = rfile,
            outfile = out_nc,
            hypertrace = list(observer_zenith = I(ozen),
                              observer_azimuth = I(oaz),
                              noisefile = I(noisefile),
                              localtime = I(lt))
          ))
          write(toJSON(ht_settings, auto_unbox = TRUE, pretty = TRUE), outfile)
        }
      }
    }
  }
}

runfile <- readLines("scripts/run_template.sh")
newfile <- gsub("@NRUNS@", ht_loops, runfile)
writeLines(newfile, "data/derived/ht-trait-inputs/run_chadwick.sh")
