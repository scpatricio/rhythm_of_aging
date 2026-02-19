# Data note (HMD)

This repository includes `data/cohort_data_HMD.csv`, which is a derived dataset used for the analyses.

The underlying Human Mortality Database (HMD) input files (e.g., `cMx_1x1.txt`, `cExposures_1x1.txt`) cannot be redistributed here due to HMD terms of use.  
To rebuild `cohort_data_HMD.csv` from raw HMD files, download the HMD data using your own access and place country folders under:

`data-raw/HMD/COUNTRY/STATS/`

Then run `R/01_prepare_data.R`.
