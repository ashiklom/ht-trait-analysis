## library() calls go here
library(conflicted)
## library(dotenv)
library(targets)
library(fs)
library(stringr)
library(dplyr)
library(raster)
library(tidyr)
library(purrr)
library(progress)
library(ggplot2)
library(ggthemes)
library(distributions3)
library(data.table)
library(forcats)

stopifnot(
  requireNamespace("reshape2", quietly = TRUE),
  requireNamespace("qs", quietly = TRUE),
  requireNamespace("hyperSpec", quietly = TRUE)
)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("between", "dplyr")
conflict_prefer("pdf", "distributions3")
conflict_prefer("quantile", "stats")

resultdir <- path("data", "derived", "ht-output")

# Source all R files
dir_ls("./R/", regexp = "\\.R$") %>%
  str_subset("plan[.]R", negate = TRUE) %>%
  walk(source)
figdir <- dir_create(path("figures", "agu"))

tar_option_set(format = "qs")

tar_pipeline(
  tar_target(
    plsr_files,
    dir_ls(resultdir, regexp = "plsr-results", recurse = TRUE)
  ),
  tar_target(plsr_info, as.data.table(parse_ht_directory(plsr_files))),
  tar_target(raw_plsr_results, rbindlist(lapply(plsr_files, read_plsr_file)),
             format = "fst_dt"),
  tar_target(input_refl_file, "data/derived/zz-traitanalysis/neon-vegetation.envi",
             format = "file"),
  tar_target(input_refl, {
    indat <- hyperSpec::read.ENVI(input_refl_file,
                                  paste0(input_refl_file, ".hdr"))
    wl <- hyperSpec::wl(indat)
    spec <- as.data.table(t(indat$spc))
    spec[, wl := wl]
    spec_long <- melt(spec, id.vars = "wl")
    spec_long[, iobs := as.integer(str_remove(variable, "V"))]
    spec_long[, variable := NULL]
    spec_long[, .(iobs, wl, value)]
  }, format = "fst_dt"),
  tar_target(true_plsr, do_true_plsr(input_refl_file), format = "fst_dt"),
  tar_target(example_spec_1, example_spec(1)),
  tar_target(example_spec_dt, tidy_example_spec(example_spec_1),
             format = "fst_dt"),
  # AOD uncertainty from single spectrum (reflectance)
  tar_target(
    aod_reflectance_fig,
    aod_uncertainty_example(example_spec_dt, input_refl),
    format = "file"
  ),
  # TODO: Retrieval uncertainty (range) from single spectrum (reflectance)
  # PLSR uncertainty from example spectrum
  tar_target(
    aod_plsr_fig_data,
    aod_uncertainty_plsr_example_data(raw_plsr_results, true_plsr)
  ),
  tar_target(
    aod_plsr_fig,
    aod_uncertainty_plsr_example_figure(aod_plsr_fig_data),
    format = "file"
  ),
  # PLSR uncertainties across all spectra
  tar_target(plsr_uncertainties, qs::qread("data/derived/plsr-variance.qs")),
  tar_target(
    plsr_uncertainty_fig,
    plsr_uncertainty_figure(plsr_uncertainties, true_plsr),
    format = "file"
  ),
  tar_target(input_refl_fig, {
    input_summary <- input_refl[, .(
      lo = min(value),
      hi = max(value)
    ), wl]
    iselected <- input_refl[between(wl, 998, 1003)][order(value)][seq.int(50, 450, 50)]
    input_sub <- input_refl[iobs %in% iselected$iobs]
    input_sub[screen_wl(wl), value := NA_real_]
    input_summary[screen_wl(wl), `:=`(lo = NA_real_, hi = NA_real_)]
    plt <- ggplot(input_summary) +
      aes(x = wl) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray70") +
      geom_line(aes(y = value, group = iobs), data = input_sub) +
      labs(x = "Wavelength (nm)", y = "Reflectance [0,1]") +
      theme_bw()
    figfile <- path(figdir, "input-reflectance.png")
    ggsave(figfile, plt, width = 5, height = 3.5, units = "in", dpi = 300)
    figfile
  }, format = "file")
  # TODO: PLSR trends with AOD
)
