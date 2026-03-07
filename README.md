# Welcome

This is a "hespdiv" package, built for hierarchical spatial data subdivision (or HespDiv in in short) and analysis. It currently contains one HespDiv method - *hespdiv*, allowing to subdivide data in 2D space. 
The methodological framework implemented is presented in:
Daumantas, L., & Spiridonov, A. (2024). hespdiv: an R package for spatially constrained, hierarchical and contiguous regionalization in palaeobiogeography. Palaeontology, 67(3), e12702. https://doi.org/https://doi.org/10.1111/pala.12702 

See the 'HespDiv Walkthrough: Case of US Miocene Mammals' vignette for the application example.

# Installing "hespdiv"
```{r eval = FALSE}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Liudas-Dau/hespdiv")
```

# Installing "hespdiv" data package

Examples in the walkthrough are based on datasets from a separate data package - "HDData".

```{r eval = FALSE}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Liudas-Dau/hespdiv_data")
```
