# ----------------------------- Load Algorithm ------------------------------- #
library(Rcpp)
library(RcppArmadillo)

# set compiler options
Sys.setenv("PKG_CXXFLAGS"=paste0(" -DARMA_64BIT_WORD=1 -I ",Sys.getenv("R_LIBS_USER"),"/Rcpp/include -I ", Sys.getenv("R_LIBS_USER"),"/RcppArmadillo/include"))

# main algorithm
source('R/stmcmclm_cpp.R')

# classes
source('R/model_input_class.R')
source('R/model_output_class.R')
source('R/mcmc_output_class.R')

# functions
# updating
sourceCpp('src/update_taue.cpp')
sourceCpp('src/update_betai.cpp')
sourceCpp('src/update_rjnufine.cpp')
sourceCpp('src/update_mhnuphi.cpp')
# helper functions
sourceCpp('src/group_data.cpp')
sourceCpp('src/save_iteration.cpp')

# some functions are sub-functions
# they are sourced by # include "subfunction.h" in updating functions
# update_wdbetai, update_rjbetai, update_lambda.h, matrix_rgamma.h, approximatelyEqualRel.h

# --------------------------- test algorithm --------------------------------- #
### Algorithm settings
nburn_in <- 100
ndraw_in <- 1000
nrun <- nburn_in *2 + ndraw_in + 1

### data
# Test Data set ("data/platinumSpike_balanced.csv") is in directory "test" and a
# subset of the Platinum Spike in data set.
# File contains additional information. Only expression values are passed to algorithm
data_full <- read.table("data/platinumSpike_balanced.csv", header = T, sep = '\t')
input_data <- as.matrix(data_full[,4:21], dimnames = NULL)
group_vec <- as.factor( c(rep(0, 9), rep(1,9) ))
GId <- row.names(data_full)


### run algorithm
set.seed(1)
# robust analysis
output_rob <- stmcmclm_cpp(input_data, group_vec, GId, numax = 45,
                       nburn = nburn_in, ndraw = ndraw_in, filename="test_rob", rob = 1)
# non-robust analysis
output_gau <- stmcmclm_cpp(input_data, group_vec, GId, numax = 45,
                           nburn = nburn_in, ndraw = ndraw_in, filename="test_gau", rob = 0)

# -------------------------- Analytics --------------------------------------- #
### A view basic analysis tools. Most are still unfinished
# parameter su only stored in sampling iterations i.e. when diagnostics functions
# need nrun as input use ndraw

# analysisBasic: 1) plot traceplots for nu and taue 2) plot su/nrun 3) density nu
source("R/analysisBasic.R")
analysisBasic(output_rob@nus, output_rob@taues, output_rob@sus, ndraw_in)

# analyse_compar: comparison between robust and non-robust model
source("R/analyse_compar.R")
analyse_compar(output_rob@sus, output_gau@sus, ndraw_in)

# analyse_GIDout and analyse_GIDout_variable: number of genes classified as DE
# in both analysis modes (shared) and only in one mode (different)
# analyse_GIDout takes file input where analyse_GIDout_variable takes variables
source("R/analyse_GIDout_variable.R")
# source("analytics/analyse_GIDout.R")
analyse_GIDout_variable(output_rob@sus, output_gau@sus, GId, cut_value=0.8, nrun, write_result = "F", outdir = "")

# analysemat: proportion of nu values within +/-1 of true nu value
source("R/analysemat.R")
analysemat(output_rob@nus, ndraw_in, 4)

# compareWrapper: calls analysemat, analyse_compar, analyse_GIDout_variable
source("R/compareWrapper.R")
compareWrapper(output_rob@sus, output_gau@sus, mean(output_rob@nus), ndraw_in, 45, 4, GId, 0.8)
