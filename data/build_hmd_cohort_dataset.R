###############################################################################
# Build cohort-level dataset from HMD cohort files (cMx_1x1, cExposures_1x1)
# Output: PNAS_cohort_data_HMD.csv (semicolon-friendly via write.csv2)
###############################################################################

# ---- Packages ----
library(data.table)
library(dplyr)

# ---- User inputs ----
# Path to the folder that contains one subfolder per HMD population
# Example structure: HMD_ROOT / "Denmark" / "STATS" / "cMx_1x1.txt"
HMD_ROOT <- "PATH/TO/HMD"  # <-- CHANGE THIS

# RData file that contains an object called `country_names`
# mapping HMD folder codes -> readable population names
COUNTRY_NAMES_FILE <- "country_names.RData"  # <-- CHANGE THIS (or keep if in working dir)

# Populations to include in the paper
countries_include <- c(
  "Denmark", "Finland", "France (Total)", "Italy", "Netherlands",
  "England & Wales (Total)", "Norway", "Sweden", "Australia",
  "Canada", "United States", "Japan"
)

# Output file
OUT_FILE <- "PNAS_cohort_data_HMD.csv"

# ---- Load country name mapping ----
load(COUNTRY_NAMES_FILE)  # expects `country_names` in the workspace

# ---- Helper: read one population ----
read_country_data <- function(country_code, hmd_root, country_names) {

  mx_file <- file.path(hmd_root, country_code, "STATS", "cMx_1x1.txt")
  ex_file <- file.path(hmd_root, country_code, "STATS", "cExposures_1x1.txt")

  if (!file.exists(mx_file) || !file.exists(ex_file)) {
    message(sprintf("Skipping %s: required file(s) not found.", country_code))
    return(NULL)
  }

  # HMD text files have 2 header rows; whitespace-delimited
  mx_data <- read.csv(mx_file, sep = "", skip = 2)
  ex_data <- read.csv(ex_file, sep = "", skip = 2)

  # Reshape to long format: columns become Sex (Male/Female/Total)
  mx_long <- melt(
    setDT(mx_data),
    id.vars = c("Year", "Age"),
    variable.name = "Sex",
    value.name = "mx"
  )

  ex_long <- melt(
    setDT(ex_data),
    id.vars = c("Year", "Age"),
    variable.name = "Sex",
    value.name = "Ex"
  )

  # Merge mortality rates with exposures
  df <- merge(as.data.frame(mx_long), as.data.frame(ex_long), by = c("Year", "Age", "Sex"))

  # Clean age and convert types
  df$Age[df$Age == "110+"] <- 110
  df$Age <- as.numeric(df$Age)

  df <- df %>%
    filter(mx != ".", Ex != ".") %>%
    mutate(
      mx = as.numeric(mx),
      Ex = as.numeric(Ex),
      Dx = ceiling(Ex * mx)  # reconstructed deaths
    ) %>%
    filter(Sex != "Total", Ex > 0) %>%
    rename(Cohort = Year) %>%
    mutate(
      Code = country_code,
      Country = country_names[[country_code]]
    )

  df
}

# ---- Read all available populations in the local HMD folder ----
country_codes <- list.files(path = HMD_ROOT)

data_list <- lapply(country_codes, read_country_data, hmd_root = HMD_ROOT, country_names = country_names)
data_full <- do.call(rbind, data_list)

# ---- Subset to paper populations and export ----
data_paper <- data_full %>%
  filter(Country %in% countries_include)

write.csv2(data_paper, OUT_FILE, row.names = FALSE)

message(sprintf("Saved: %s", OUT_FILE))
