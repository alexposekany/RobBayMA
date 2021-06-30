# --------------------- Description ---------------------- #
# important diagnostic plots 
#----------------------------------------------------------#

# ---------------------- DEV STATUS ---------------------- #
# rename variables
# specify return value/print statement
# decide to use histogram or desityplot
# -------------------------------------------------------- #


# Basic Analysis Function
analysisBasic <- function(nu, taue, su, nrun, title = ""){
  library(coda)
  # nu mean & median
  # print(paste("nu mean: ", round(mean(nu), 2)))
  # print(paste("nu median: ", round(median(nu), 2)))
  
  # traceplots nu & taue 
  par(mfrow=c(2,2))
  coda::traceplot(as.mcmc(nu), ylab = "nu", smooth = F)
  coda::traceplot(as.mcmc(taue), ylab = "taue", smooth = F)
  # Fraction of su
  plot(sort(su/nrun, decreasing = T), ylab="su/nrun")
  # histogram of nu
  # hist(nu, freq = F, main = "")
  # abline(v=mean(nu), col='blue')
  # abline(v=median(nu), col='red')
  densplot(as.mcmc(nu))
  mean_txt <- paste("mean: ", round(mean(nu), 2) )
  md_txt <- paste("median: ", round(median(nu), 2) )
  # browser()
  txt <- paste(mean_txt, md_txt, sep = '\n')
  legend("topleft", legend = txt, bty='n')
  mtext(title, side = 3, line = -2, outer = TRUE)
  par(mfrow=c(1,1))
  }


# analyseGIDout <- function(filename, GID, cut){
#   
# }
