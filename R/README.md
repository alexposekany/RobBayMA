# RobBayMA

## Notes

1. approximatelyEqualRel.cpp function needs copy right
taken from [1]

2. parameter su only stored in sampling iterations i.e. when diagnostics functions
need nrun as input use ndraw

## TODOs
[ ] add roxygen2-complient comments to ALL functions [4]

[ ] please read and write data (csv, .RData, PDF) to the `data`-directory (best practice [2])  

[x] PKG_CXXFLAGS should be set via makefile or something similar [3] **Could be done by adding "LinkingTo" in DESCRIPTION [5]**


## Sources
[1] https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/  

[2] https://r-pkgs.org/data.html  

[3] https://cran.r-project.org/doc/manuals/R-exts.html#Configure-and-cleanup  

[4] https://r-pkgs.org/man.html

[5] https://r-pkgs.org/src.html#make
