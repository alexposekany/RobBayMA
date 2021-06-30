# --------------------- Description ---------------------- #
# Wrapper function for all analysis functions
#----------------------------------------------------------#

# ---------------------- DEV STATUS ---------------------- #
# correct sourcing
# -------------------------------------------------------- #

compareWrapper <- function(suROB, suGAU, nuROB, nrun, numax, nutrue, GID, cut_value){
  source('R/analyse_compar.R')
  source('R/analyse_GIDout_variable.R')
  source('R/analysemat.R')

  l_out <- analyse_GIDout_variable(suROB, suGAU, GID, cut_value, nrun)
  l_out["nutrue"] <- analysemat(nuROB, nrun, nutrue)
  analyse_compar(suROB, suGAU, nrun, numax)

  return(l_out)
}
