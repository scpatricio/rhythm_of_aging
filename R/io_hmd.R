# R/io_hmd.R

read_country_data <- function(country_code, hmd_root, country_names_map) {
  mx_file <- file.path(hmd_root, country_code, "STATS", "cMx_1x1.txt")
  ex_file <- file.path(hmd_root, country_code, "STATS", "cExposures_1x1.txt")

  if (!file.exists(mx_file) || !file.exists(ex_file)) {
    message(sprintf("Skipping %s: file not found.", country_code))
    return(NULL)
  }

  mx_data <- read.csv(mx_file, sep = "", skip = 2)
  ex_data <- read.csv(ex_file, sep = "", skip = 2)

  Ex <- data.table::melt(data.table::setDT(ex_data),
                         id.vars = c("Year", "Age"),
                         variable.name = "Sex",
                         value.name = "Ex")
  mx <- data.table::melt(data.table::setDT(mx_data),
                         id.vars = c("Year", "Age"),
                         variable.name = "Sex",
                         value.name = "mx")

  df <- merge(as.data.frame(mx), as.data.frame(Ex), by = c("Year", "Age", "Sex"))

  df$Age[df$Age == "110+"] <- 110
  df$Age <- as.numeric(df$Age)

  df <- df |>
    dplyr::filter(mx != "." & Ex != ".") |>
    dplyr::mutate(
      Dx = ceiling(as.numeric(Ex) * as.numeric(mx)),
      Ex = as.numeric(Ex),
      mx = as.numeric(mx)
    ) |>
    dplyr::filter(Sex != "Total", Ex > 0) |>
    dplyr::rename(Cohort = Year)

  df$Code <- country_code
  df$Country <- country_names_map[[country_code]]

  df
}

build_hmd_cohort_dataset <- function(hmd_root, country_names_rdata) {
  load(country_names_rdata)  # loads `country_names`
  country_codes <- list.files(path = hmd_root)

  data_list <- lapply(country_codes, read_country_data,
                      hmd_root = hmd_root,
                      country_names_map = country_names)

  dplyr::bind_rows(data_list)
}
