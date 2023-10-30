# Welcome

This is a "hespdiv" package, built for hierarchical spatial data subdivision (or HespDiv in in short) and analysis. It currently contains one HespDiv method - *hespdiv*, allowing to subdivide data in 2D space. See the 'HespDiv: Concepts and Methodology' vignette for the explanation of the method and the 'HespDiv Walkthrough: Case of US Miocene Mammals' vignette for the application example.

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
