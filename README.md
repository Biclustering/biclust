
<!-- README.md is generated from README.Rmd. Please edit that file -->
Biclust R package
=================

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
