# tsatools
Automates unit root testing and procedures for multivariate time series analysis, including model diagnostics and calculations of inferential statistics (e.g., long-run multipliers, impulse response functions, standard errors)

## Installation
To install `tsatools`, you can run the following R code:
```
install.packages("devtools")
devtools::install_github("nathanamorse/tsatools")
```

## Usage

See the [vignette](https://nmorse.com/code/tsatools/) for details on using this package. It includes the following functions:

- `adf_tests`: runs a series of augmented Dickey-Fuller tests with various lags (adapted from code from Suzanna Linn)
- `bg_tests`: runs tests for serial correlation
- `irfs`: calculates impulse response functions from error correction models
- `lrms`: calculates long-run multipliers from error correction models
- `phi_tests`: conducts tests for unit roots when the model deterministics are unknown
- `ur_tests`: conducts tests for unit roots when the model deterministics are known
