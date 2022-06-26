
# CovCoExpNets

<!-- badges: start -->
<!-- badges: end -->

We integrate glmnet within the CovCoExpNets R package to create covariate-specific GCNs based on relevant predictors (i.e. predictors with non-zero coefficients) as detected by the LASSO regularization mechanism. Each hub gene leads to exactly a new module in the covariate-specific GCN. To complete the composition of those gene modules, for each hub gene s, we create a basic linear model age ~ s + b0 and keep adding genes to the model, choosing amongst the genes with highest correlation with the hub s. We do this while the R2 of the model improves. Each new module´s gene set is annotated with gProfiler2 with the sources GO, KEGG and REAC.

We applied this new tool on a snRNA-seq dataset human post-mortem substantia nigra pars compacta tissue of 13 controls and 14 Parkinson’s disease (PD) cases (18 males and 9 females) with 30-99 years. For more details, please see the link to the manuscript:

## Installation


You can install the development version of CovCoExpNets like so:

``` r
devtools::install_github('aliciagp/CovCoExpNets')
```


## Credits

The development of this package has been possible under the supervision of Juan A. Botía from the University of Murcia and Sebastian Guelfi and many collaborators from Verge Genomics start-up.




