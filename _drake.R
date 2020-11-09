## Load your packages, e.g. library(drake).
source("./packages.R")

## Load your R files
invisible(
  lapply(list.files("./R", ".R$", full.names = TRUE), source)
)

## _drake.R must end with a call to drake_config().
## The arguments to drake_config() are basically the same as those to make().
## lock_envir allows functions that alter the random seed to be used. The biggest
## culprits of this seem to be interactive graphics e.g. plotly and mapdeck.
## options("clustermq.scheduler" = "multicore")
## future::plan(future::multicore)
drake_config(the_plan,
             ## parallelism = "clustermq",
             ## parallelism = "future",
             ## jobs = 4,
             ## max_expand = 50,
             lock_envir = FALSE)

## drake::drake_cache()$unlock()
