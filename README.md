# pseudoabs_transferability

Code accompanying the manuscript **Temporal and spatial transferability in telemetry-based dynamic species distribution models: The effects of algorithms and pseudo-absence techniques**

## Overview

Species distribution models often rely on pseudo-absence data when true absence observations are unavailable. Numerous pseudo-absence generation methods and algorithms exist, yet relatively little is known about how these choices affect a model's ability to transfer across novel spatial and temporal domains.

This repository provides an end-to-end workflow to:

* Generate multiple pseudo-absence datasets;
* Fit a range of statistical and machine-learning SDMs;
* Quantify spatial and temporal transferability;
* Assess environmental extrapolation;
* Compare model performance across pseudo-absence strategies and modelling frameworks.

---

## Repository Structure

```text
pseudoabs_transferability/
│
├── R/                                 # Supporting R functions
├── figures/                           # Manuscript figures and outputs
│
├── 00_thinning.R                      # Data thinning
├── 01_create_background_samples.R     # Background sampling
├── 02_create_CRWs.R                   # Correlated random walks
├── 03a_calculate_step_lengths.R       # Step length calculation to inform buffer sampling
├── 03b_create_buffer_samples.R        # Buffer sampling
├── 04_extract_variables.R             # Extract environmental variables
│
├── 05_temporal_*.R                    # Temporal transferability analyses
├── 06_spatial_*.R                     # Spatial transferability analyses
│
├── 99_combine_spatial_outputs.R       # Combine all spatial transfer test results
├── 99_combine_temporal_outputs.R      # Combine all temporal transfer test results
├── 99_increase_temporal_distance_*.R  # Supplementary analyses of results when temporal interval is increased
└── 99_sample_sizes.R                  # Compute sample sizes of tests
```


## Requirements

Analyses are conducted entirely in **R**.

---

## Citation

If you use this repository, please cite the associated manuscript:
Wilson et al. (2026) Temporal and spatial transferability in telemetry-based dynamic species distribution models: The effects of algorithms and pseudo-absence techniques. _Ecological Modelling_ [ADD DOI WHEN CREATED]

---

## Author

**Joshua C. Wilson**
PhD Researcher
University of Southampton and British Antarctic Survey
