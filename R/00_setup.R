# R/00_setup.R
rm(list = ls(all = TRUE))
set.seed(1)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(rstan)
  library(doParallel)
  library(bayesplot)
  library(HDInterval)
  library(modeest)
  library(here)
})

# Parallel setup
total_cores <- parallel::detectCores()
cl <- makeCluster(max(1, total_cores - 1))
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

options(mc.cores = total_cores)
rstan_options(auto_write = TRUE)
