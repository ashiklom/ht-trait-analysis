label_atmosphere <- function(atm, aod, h2o) {
  sprintf(
    "%sA%dW%d",
    str_to_upper(str_sub(atm, 1, 1)),
    round(aod * 100),
    round(h2o * 10)
  )
}

factor_pls <- function(plsr_var) {
  factor(
    plsr_var,
    c("pct_N", "pct_C", "M_area"),
      ## "pct_Cellulose", "pct_ADL", "pct_ADF",
      ## "del15N"),
    c("Leaf %N", "Leaf %C", "LMA")
      ## "Leaf %Cellulose")
  )
}
