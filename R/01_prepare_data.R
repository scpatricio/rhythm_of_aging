# R/01_prepare_data.R
source(here::here("R", "00_setup.R"))

# ---- Option A: load shared derived data (default)
csv_path <- here::here("data", "cohort_data_HMD.csv")
if (file.exists(csv_path)) {
  data_full <- read.csv2(csv_path)
  message("Loaded derived dataset: data/cohort_data_HMD.csv")
  return(invisible(data_full))
}

# ---- Option B: rebuild from HMD (requires user-provided files)
# NOTE: HMD raw files cannot be redistributed; user must download separately.
hmd_root <- here::here("data-raw", "HMD")
if (!dir.exists(hmd_root)) {
  stop("No data found. Either add data/cohort_data_HMD.csv or place HMD files in data-raw/HMD/.")
}

# You used country_names.RData; keep it project-local
load(here::here("data-raw", "country_names.RData"))

countries <- list.files(path = hmd_root)

read_country_data <- function(country) {
  mx_file <- file.path(hmd_root, country, "STATS", "cMx_1x1.txt")
  ex_file <- file.path(hmd_root, country, "STATS", "cExposures_1x1.txt")

  if (!file.exists(mx_file) || !file.exists(ex_file)) {
    message(sprintf("Skipping %s: file not found.", country))
    return(NULL)
  }

  mx_data <- read.csv(mx_file, sep = "", skip = 2)
  ex_data <- read.csv(ex_file, sep = "", skip = 2)

  Ex <- melt(setDT(ex_data), id.vars = c("Year", "Age"),
             variable.name = "Sex", value.name = "Ex")
  mx <- melt(setDT(mx_data), id.vars = c("Year", "Age"),
             variable.name = "Sex", value.name = "mx")

  df <- merge(as.data.frame(mx), as.data.frame(Ex), by = c("Year", "Age", "Sex"))
  df$Age[df$Age == "110+"] <- 110
  df$Age <- as.numeric(df$Age)

  df <- df %>%
    filter(mx != "." & Ex != ".") %>%
    mutate(
      Dx = ceiling(as.numeric(Ex) * as.numeric(mx)),
      Ex = as.numeric(Ex),
      mx = as.numeric(mx)
    ) %>%
    filter(Sex != "Total", Ex > 0) %>%
    rename(Cohort = Year)

  df$Code <- country
  df$Country <- country_names[[country]]
  df
}

data_list <- lapply(countries, read_country_data)
data_full <- do.call(rbind, data_list)

dir.create(here::here("data"), showWarnings = FALSE, recursive = TRUE)
write.csv2(data_full, here::here("data", "cohort_data_HMD.csv"), row.names = FALSE)
message("Rebuilt and saved: data/cohort_data_HMD.csv")

invisible(data_full)
