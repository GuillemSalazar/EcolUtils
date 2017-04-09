#**EcolUtils**: Utilities for community ecology analysis

The R package **EcolUtils** provides tools and extensions for community ecology analysis not present in up-to-date available packages. Designed with molecular 16S-derived data in mind.

**EcolUtils** depends on **vegan** and **spaa** package which can be installed from CRAN.


##Content
The current version contains the following functions:

+ **rrarefy.perm()** 
Rarefaction of a community matrix with permutations
+  **adonis.pair()** 
Pairwise comparisons for Permutational Multivariate Analysis of Variance Using Distance Matrices
+ **spec.gen()**
Specialist/Generalist classification of OTUs based on niche width and permutation algorithms
+ **seasonality.test()**
Seasonality test based on autocorrelation and null communities
+ **niche.val()^^
Niche value computation for OTUs in a community through abundance-weighted mean and matrix randomization

##Authors

Guillem Salazar

##Citation

To see the preferable citation of the package, type:
```r
citation("EcolUtils")
```
##Software Versions

*R*: version 3.2.1  
*RStudio*: Version 0.99.467

##Contact Info

Guillem Salazar (salazar@icm.csic.es)


## Installation of the package

* Install the latest version of **vegan** and **spaa** from Bioconductor:
```r
install.packages("vegan")
install.packages("spaa")
```
* Install EcolUtils' current development version from Github:
```r
devtools::install_github("GuillemSalazar/EcolUtils")
```
* After attaching the package, you are ready to start:
```r
library(EcolUtils)
```


