# HespDiv

`hespdiv` is an R package implementing the **HespDiv framework for hierarchical spatial data subdivision**.
The package provides tools for identifying topologically contiguous spatial clusters and analysing hierarchical regionalization structures in palaeobiogeography and related spatial datasets.

HespDiv represents a suitable framework for identifying hierarchically nested,
spatially contiguous regions (geobiomes *sensu* Spiridonov & Eldredge(2024)) from taxonomic occurrence data.

The theoretical foundations of HespDiv are laid out in my PhD dissertation (link to be added) and in:

Daumantas, L. & Spiridonov, A. (2024).
hespdiv: an R package for spatially constrained, hierarchical and contiguous regionalization in palaeobiogeography.
Palaeontology, 67(3), e12702.
https://doi.org/10.1111/pala.12702

Theoretical foundations of the HespDiv are lied out in my Phd dissertation (link to be added!) and in

Spiridonov, A., & Eldredge, N. (2024). 
The Bretskyan hierarchy, multiscale allopatry, and geobiomes—on the nature of evolutionary things. 
Paleobiology, 1-20.
https://doi.org/10.1017/pab.2023.37 

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
