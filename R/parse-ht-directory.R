parse_ht_directory <- function(fname) {
  tibble::tibble(
    fname = fname,
    atm = str_extract(fname, "(?<=atm_).*?(?=__)"),
    aod = as.numeric(str_extract(fname, "(?<=aod_)[[:digit:].]+(?=__)")),
    h2o = as.numeric(str_extract(fname, "(?<=h2o_)[[:digit:].]+(?=/)")),
    szen = as.numeric(str_extract(fname, "(?<=szen_)[[:digit:].]+(?=__)")),
    ozen = as.numeric(str_extract(fname, "(?<=ozen_)[[:digit:].]+(?=__)")),
    saz = as.numeric(str_extract(fname, "(?<=saz_)[[:digit:].]+(?=__)")),
    oaz = as.numeric(str_extract(fname, "(?<=oaz_)[[:digit:].]+(?=/)")),
    snr = as.numeric(str_extract(fname, "(?<=snr_)[[:digit:].]+(?=/)")),
    calibration = str_extract(fname, "(?<=cal_).*?(?=__)"),
    calscale = as.numeric(str_extract(fname, "(?<=scale_)[[:digit:].]+(?=/)")),
    calid = as.numeric(str_extract(path_ext_remove(fname), "(?<=-)[[:digit:]]{2}$"))
  )
}
