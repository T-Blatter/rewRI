# rewRI: Reference Interval Estimation from the Re-Weighted BC or YJ-Transformation

This package supports the fully automated estimation of direct reference intervals (RIs) from clinical cohort data. An accompanying graphical output further helps to assess the accuracy of the inferred estimates.

## Overview

**rewRI** provides a robust statistical framework for estimating Reference Intervals (RIs) by combining power transformations with adaptive reweighting. By utilizing **Box-Cox (BC)** or **Yeo-Johnson (YJ)** transformations alongside the weighting schemes from the `cellWise` package, it effectively handles non-normal distributions and mitigates the impact of pathological outliers in clinical datasets.

This package implements the methodology detailed in:
> Blatter, T.U.; Nakas, C.T.; Leichtle, A.B. (2024). Direct, age- and gender-specific reference intervals: applying a modified M-estimator of the Yeo-Johnson transformation to clinical real-world data. _JLM_. [doi:10.1515/labmed-2024-0076](https://doi.org/10.1515/labmed-2024-0076)

## Installation

You can install the **rewRI** R package directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("T-Blatter/rewRI")

```

## Quick Start

The core function `rewRI()` allows you to choose between transformation types and automatically computes the re-weighted interval.

```r
library(rewRI)

# Example sample data with potential outlier(s) and sample skewness
obs <- c(rlnorm(200, meanlog = 2, sdlog = 0.5), rnorm(20, mean = 25, sd = 1))

# Estimate reweighted RI using Yeo-Johnson transformation directly from the data
ri_results <- rewRI(obs, method = "YJ", conf_level = 0.95)

# Summary of the estimated limits
summary(ri_results)

```

## Requirements

* **R** (>= 4.0.0)
* **cellWise** (>= 2.5.2)

## Citation


```bibtex
@article{BlatterNakasLeichtle+2024+239+250,
	url = {https://doi.org/10.1515/labmed-2024-0076},
	title = {Direct, age- and gender-specific reference intervals: applying a modified M-estimator of the Yeo-Johnson transformation to clinical real-world data},
	author = {Tobias Ueli Blatter and Christos Theodoros Nakas and Alexander Benedikt Leichtle},
	pages = {239--250},
	volume = {48},
	number = {5},
	journal = {Journal of Laboratory Medicine},
	doi = {doi:10.1515/labmed-2024-0076},
	year = {2024}
}
```

## License

This package is licensed under the **GPL-3** License.