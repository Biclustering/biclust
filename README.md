
<!-- README.md is generated from README.Rmd. Please edit that file -->
Biclust R package
=================

Installation
------------

#### Release (1.2.0)

``` r
install.packages("biclust")
```

#### Development (1.2.1)

``` r
install.packages("devtools") # If not yet installed on your R Version
devtools::install_github("Biclustering/biclust")
```

Description
-----------

The main function biclust provides several algorithms to find biclusters in two-dimensional data: Cheng and Church, Spectral, Plaid Model, Xmotifs and Bimax. In addition, the package provides methods for data preprocessing (normalization and discretisation), visualisation, and validation of bicluster solutions.

Example
-------

``` r
data(BicatYeast)
res<-biclust(BicatYeast, method=BCPlaid(), verbose=FALSE)
res
```
