# Implementing Convex Optimization in R: Two Econometric Examples

This repository hosts code for estimation methods in 

* Zhan Gao and Zhentao Shi: [*Implementing Convex Optimization in R: Two Econometric Examples*](https://arxiv.org/abs/1806.10423).

After local installation of the convex solver [`Mosek`](https://docs.mosek.com/9.1/install/index.html) and the R package `Rmosek`, one can  run the R script `gao-shi-rmosek.R` to replicate the simulation studies.

## Install Rmosek

First, `MOSEK` must have been installed.

Next, in `R` (verified with success in R 3.6.3):

```{r}
install.packages("Rmosek")
library(Rmosek)
mosek_attachbuilder("path_to_the_bin_folder_of_MOSEK")
install.rmosek()
```
It is done.
