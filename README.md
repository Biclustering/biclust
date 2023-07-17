
<!-- README.md is generated from README.Rmd. Please edit that file -->
Biclust R package
=================

General
------------
This repository is archived since the development will be done in CRAN directly (https://cran.r-project.org/web/packages/biclust/index.html) . This is need due to a change in the lead maintainer which will be Friedrich Leisch from now on.

Installation
------------

#### Release (2.0.3)

``` r
install.packages("biclust")
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
