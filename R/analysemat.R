# --------------------- DESCRIPTION ---------------------- #
# analysemat1: percent of genes with correctly classifed nu
#              calls analysisBasis.R
# analysemat2: percent of genes with correctly classifed nu
# ---------------------------------------------------------#

# ---------------------- DEV STATUS ---------------------- #
# missing: final touch (sourcing, var names, return values)
# -------------------------------------------------------- #


analysemat <- function(nu, nrun, nutrue){
  t <- round(sum((nutrue-1)<nu&nu<(nutrue+1))/nrun, 2)
  # print(paste('within true value +/-1: ', t))
  return(t)
}