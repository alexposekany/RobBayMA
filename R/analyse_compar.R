# --------------------- DESCRIPTION ---------------------- #
# Visual comparison of two su outputs
# 1st graph: course of sus as line diagram  
# 2nd graph: course of su_1 as line diag and su_2 as points
# 3nd graph: course of su_2 as line diag and su_1 as points
# Interpretation: how similar are differnt su values 
# ---------------------------------------------------------#

# ---------------------- DEV STATUS ---------------------- #
# missing: final comparison with ML
# missing: plots in gray-scale
# missing: captions
# missing: final touch
# -------------------------------------------------------- #



analyse_compar <- function(suROB, suGAU, nrun, numax){
  # browser()
  # load('analysisBasic.R')
  
    
    # ML: m wird returniert
  # m <- mean(nu) 
  # suROB <- su
  # nrun<-dim(nu)[1]
  # in ML: write nu & taue to file & call analysebasic
  # dlmwrite(strcat(strcat(char(filename),'_rob'),'nu.txt'),nu)
  # dlmwrite(strcat(strcat(char(filename),'_rob'),'taue.txt'),taue)
  # analysebasic(nrun,su,taue,nu,numax)
  # suGAU<-su
  # in ML: write nu & taue to file & call analysebasic
  # dlmwrite(strcat(strcat(char(filename),'_norm'),'nu.txt'),nu)
  # dlmwrite(strcat(strcat(char(filename),'_norm'),'taue.txt'),taue)
  # analysebasic(nrun,su,taue,nu,numax)
  
  # plot T-4 & normal su 
  par(mfrow=c(1,3))
  # 1st plot
  plot(sort(suROB/nrun, decreasing = T), lty = "solid", ylab="su/nrun", type='l')
  lines(sort(suGAU/nrun, decreasing = T), lty = "dashed")
  legend_names <- c("rob", "normal")
  legend("bottomleft",legend=legend_names, lty = c('solid', 'dashed'), bty='n', lwd = 1) 
  # 2nd plot
  Gnorm_sort<-sort(suGAU/nrun, decreasing = T, index.return=T)
  Grob<-suROB[Gnorm_sort$ix]/nrun
  plot(sort(suGAU/nrun, decreasing = T), col = 'red', ylab="su/nrun", type = 'l')
  lines(Grob, col = 'blue', type = "p")
  legend_names <- c( "normal", "rob")
  legend("bottomleft",legend=legend_names, col = c('red', 'blue'),lty=1, bty='n', lwd = 1) 
  # 3rd plot
  Grob_sort<-sort(suROB/nrun, decreasing = T, index.return=T)
  Gnorm<-suGAU[Grob_sort$ix]/nrun
  plot(sort(suROB/nrun, decreasing = T), col = 'blue', ylab="su/nrun", type = 'l')
  lines(Gnorm, col = 'red', type = "p")
  legend_names <- c( "rob", "normal")
  legend("bottomleft",legend=legend_names, col = c('blue', 'red'),lty=1, bty='n', lwd = 1) 
  par(mfrow=c(1,1))
  # return(m)
  
}