## RobBayMA

A repository to store and manage the code of the robbayma package.

#### Setup

The installation currently can only be run manually by cloning the repo, running 

```
Sys.setenv("PKG_CXXFLAGS"=paste0(" -DARMA_64BIT_WORD=1 -I ",Sys.getenv("R_LIBS_USER"),"/Rcpp/include -I ", Sys.getenv("R_LIBS_USER"),"/RcppArmadillo/include"))
```

and then building it via command line or RStudio.
