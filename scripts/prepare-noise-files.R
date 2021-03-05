library(conflicted)
library(vroom)
library(dplyr)

noise_400 <- vroom("data/sbg_noise_coeffs.txt", skip = 1,
                   col_names = c("#wvl", "a", "b", "c", "rmse"))

noise_200 <- noise_400 %>%
  mutate(across(c(a, b, c, rmse), ~.x * 2))

noise_300 <- noise_400 %>%
  mutate(across(c(a, b, c, rmse), ~.x * 4/3))

vroom_write(noise_400, "data/derived/ht-trait-inputs/noise-400.txt", delim = " ")
vroom_write(noise_300, "data/derived/ht-trait-inputs/noise-300.txt", delim = " ")
vroom_write(noise_200, "data/derived/ht-trait-inputs/noise-200.txt", delim = " ")
