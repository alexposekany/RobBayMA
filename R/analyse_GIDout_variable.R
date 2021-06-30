# --------------------- Description ---------------------- #
# returns number shared and different genes
#   shared genes: genes passing threshold in both data sets
#   different genes: genes passing threshold in one data set
# Settings: write_results = T to write reults to file
#   all_f1: all genes and p-values identified in data set1
#   all_f2: all genes and p-values identified in data set2
#   above_f1: all genes & p-values passing threshold in data set1  
#   above_f2: all genes & p-values passing threshold in data set2
#----------------------------------------------------------#

# ---------------------- DEV STATUS ---------------------- #
# rename variables
# specify return value/print statement
# -------------------------------------------------------- #

analyse_GIDout_variable <- function(su1_in, su2_in, GID, cut_value, nrun, write_result = "F", outdir = ""){
  
    # analyse file1
  su1 <- su1_in/nrun
  sort_index1 <- sort(su1, decreasing = TRUE, index.return=TRUE) # 1st element-> sorted values , 2nd element -> index
  GIsort1 <- GID[sort_index1$ix]
  logical_cut1 <- sort_index1$x>=cut_value
  s1_passingcut <- sort_index1$x[logical_cut1] # wird nicht weiter verwendet
  GI_1 <- GIsort1[logical_cut1]
  
  # analyse file2 
  su2 <- su2_in/nrun
  sort_index2 <- sort(su2, decreasing = TRUE, index.return=TRUE) # 1st element-> sorted values , 2nd element -> index 
  GIsort2 <- GID[sort_index2$ix]
  logical_cut2 <- sort_index2$x>=cut_value
  s2_passingcut <- sort_index2$x[logical_cut2]
  GI_2 <- GIsort2[logical_cut2]
  
  # Analyse Differences between robust and normal analysis
  G_identical <- intersect(GI_1, GI_2)
  G_different <- c(setdiff(GI_1, GI_2), setdiff(GI_2, GI_1)) # 1st expr.: items in 1st set not present in second set, 2md expr.: items in 2nd set not present in 1st set  
  # print(paste("shared genes: ", length(G_identical)))
  # print(paste("different genes: ", length(G_different)))
  
  # output list
  out_l <- list(shared = length(G_identical), different = length(G_different))
  # browser()
  
  if (write_result == TRUE){
    # All genes with p-values
    # file 1
    ch <- as.character(t(GIsort1))
    m <- matrix(0, nrow = length(GIsort1), ncol = 2)
    for (i in 1:dim(m)[1]){
      m[i,] <- c(as.integer(GIsort1[i]), round(sort_index1$x[i], 5))
    }
    # browser()
    filename_f1 = paste(outdir, "_all", "_f1.txt", sep = "")
    write.table(m, file = filename_f1, quote = T, sep = "\t", row.names = F, col.names = c("GID", "p-val"))

    # file2
    ch <- as.character(t(GIsort2))
    m <- matrix(0, nrow = length(GIsort2), ncol = 2)
    for (i in 1:dim(m)[1]){
      m[i,] <- c(as.integer(GIsort2[i]), round(sort_index2$x[i], 5))
    }
    filename_f2 = paste(outdir, "_all", "_f2.txt", sep = "")
    write.table(m, file = filename_f2, quote = T, sep = "\t", row.names = F, col.names = c("GID", "p-val"))

    # genes passing cutoff with p-values
    # file 1
    ch <- as.character(t(GI_1))
    m <- matrix(0, nrow = length(GI_1), ncol = 2)
    for (i in 1:dim(m)[1]){
      m[i,] <- c(as.integer(GI_1[i]), round(sort_index1$x[i], 5))
    }
    filename_f1 = paste(outdir, "_above_", cut_value, "_f1.txt", sep = "")
    write.table(m, file = filename_f1, quote = T, sep = "\t", row.names = F, col.names = c("GID", "p-val"))

    # file 2
    ch <- as.character(t(GI_2))
    m <- matrix(0, nrow = length(GI_2), ncol = 2)
    for (i in 1:dim(m)[1]){
      m[i,] <- c(as.integer(GI_2[i]), round(sort_index2$x[i], 5))
    }
    filename_f2 = paste(outdir, "_above_", cut_value, "_f2.txt", sep = "")
    write.table(m, file = filename_f2, quote = T, sep = "\t", row.names = F, col.names = c("GID", "p-val"))
  }
    # 
    # To Do: 
    # filename
    # dir 
    # vergleich mit ml output 
    
  
  return(out_l)
}


