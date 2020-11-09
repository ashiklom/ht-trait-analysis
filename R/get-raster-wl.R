get_raster_wl <- function(b) {
  names(b) %>%
    str_extract("(?<=[.][.])[[:digit:].]+[[:digit:]]+") %>%
    as.numeric()
}
