[![CRAN status](https://www.r-pkg.org/badges/version/hespdiv)](https://cran.r-project.org/package=hespdiv)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)

# HespDiv

hespdiv is an R package implementing the HespDiv framework for hierarchical spatial data subdivision.
The package provides tools for identifying topologically contiguous spatial clusters and analysing hierarchical regionalization structures in palaeobiogeography and related spatial datasets.

The methodological framework is described in:

Daumantas, L. & Spiridonov, A. (2024).
hespdiv: an R package for spatially constrained, hierarchical and contiguous regionalization in palaeobiogeography.
Palaeontology, 67(3), e12702.
https://doi.org/10.1111/pala.12702
See the 'HespDiv Walkthrough: Case of US Miocene Mammals' vignette for the application example.

# Installation
You can install the development version from GitHub:
```{r eval = FALSE}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Liudas-Dau/hespdiv")
```
# Example and walkthrough
A full demonstration of the workflow is provided in the vignette:

[**“HespDiv Walkthrough: Case of US Miocene Mammals”**](https://htmlpreview.github.io/?https://github.com/Liudas-Dau/hespdiv/blob/master/inst/Walkthrough.html)

The vignette illustrates how to:

- prepare spatial occurrence data
- perform hierarchical spatial subdivision
- analyse subdivision results
- visualise hierarchical regionalization patterns

# Example datasets

The walkthrough uses example datasets distributed in a separate package:

```{r eval = FALSE}
# install devtools is you don't have it!
devtools::install_github("Liudas-Dau/hespdiv_data")
```

# Citation
If you use *hespdiv* in your research, please cite:

Daumantas, L. & Spiridonov, A. (2024).
hespdiv: an R package for spatially constrained, hierarchical and contiguous regionalization in palaeobiogeography.
Palaeontology, 67(3), e12702.
https://doi.org/10.1111/pala.12702
