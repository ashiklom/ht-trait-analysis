library(ncdf4)
library(fs)
library(data.table)
library(jsonlite)
library(magrittr)
library(purrr)
library(disk.frame)

setup_disk.frame()
options(future.globals.maxSize = Inf)

df_outdir <- dir_create("./data/diskframe")

raw_output_dir <- "~/projects/sbg-uncertainty/isofit/examples/py-hypertrace/output/ht-trait-analysis"
ncfiles <- dir_ls(raw_output_dir, glob = "*.nc")

parse_atmosphere <- function(l) {
  atm <- l$atm_aod_h2o
  lout <- l
  lout$atm_aod_h2o <- NULL
  modifyList(lout, list(
    atm_profile = atm[1],
    aod = as.numeric(atm[2]),
    h2o = as.numeric(atm[3])
  ))
}

tidy_ht_info <- function(ht_vec, ht_config) {
  ht_config_l <- fromJSON(ht_config$value)
  ht_vec_l <- map(ht_vec, fromJSON) %>%
    map(setNames, names(ht_config_l$hypertrace)) %>%
    map(parse_atmosphere)
  out <- rbindlist(ht_vec_l)
  out[, iht := .I]
  return(out)
}

out_rad_path <- path(df_outdir, "toa_radiance.df")
out_rad <- disk.frame(out_rad_path)
nrad <- nchunks(out_rad)

out_refl_path <- path(df_outdir, "estimated_reflectance.df")
out_refl <- disk.frame(out_refl_path)
nrefl <- nchunks(out_refl)

## for (ncf in ncfiles) {
for (i in seq_along(ncfiles)) {
  ncf <- ncfiles[i]
  message(ncf)
  do_rad <- i > nrad
  do_refl <- i > nrefl
  if (!(do_rad || do_refl)) {
    message("Skipping iteration ", i)
  }
  nc <- nc_open(ncf)

  ht_vec <- ncvar_get(nc, "hypertrace")
  ht_config <- ncatt_get(nc, 0, "hypertrace_config")
  ht_info <- tidy_ht_info(ht_vec, ht_config)

  # hypertrace, wavelength, sample
  if (do_rad) {
    rad_raw <- ncvar_get(nc, "toa_radiance")
    rad_tidy <- as.data.table(rad_raw)
    setnames(rad_tidy, c("iht", "wavelength", "isample", "value"))
    rad_tidy2 <- ht_info[rad_tidy, on = "iht"]
    rad_df <- add_chunk(out_rad, rad_tidy2[,!"iht"])
    rm(rad_tidy, rad_tidy2, rad_raw)
  } else {
    message("Skipping existing rad at iteration ", i)
  }

  if (do_refl) {
  # same for reflectance estimate
  refl_raw <- ncvar_get(nc, "estimated_reflectance")
  refl_tidy <- as.data.table(refl_raw)
  setnames(refl_tidy, c("iht", "wavelength", "isample", "value"))
  refl_tidy2 <- ht_info[refl_tidy, on = "iht"]
  refl_df <- add_chunk(out_refl, refl_tidy2[,!"iht"])
  rm(refl_tidy, refl_tidy2, refl_raw)
  } else {
    message("Skipping existing refl at iteration ", i)
  }
  nc_close(nc)
}

message("Processing complete!")
