# --------------------- Description ---------------------- #
# File version of analyse_GID_variable.R
#----------------------------------------------------------#

# ---------------------- DEV STATUS ---------------------- #
# rename variables
# specify return value/print statement
# -------------------------------------------------------- #

robust_in <- '../mcmc_profiles/analyse_GIDout_rob.RData' # test 2020_02_05.RData
norm_in <- '../mcmc_profiles/analyse_GIDout_norm.RData' # test 2020_02_06.RData
numax_in <- 45
cut_value_in <- 0.85
nrun_in <- 500 # Achtung: Fix für nicht test daten
GID_in <- as.character(1:500) # Achtung: Fix für nicht test daten


analyse_GIDout <- function(file1, file2, GID, cut_value, nrun){
  # load and anylse robust Data
  load(file1)
  su_1 <- su/nrun
  # list element 1: sorted su values , list element 2: indexes
  su1_sort_index <- sort(su_1, decreasing = TRUE, index.return=TRUE)
  su1_GI <- GID[su1_sort_index$ix]
  su1_cut_logical <- su1_sort_index$x>=cut_value
  su1_passing <- su1_sort_index$x[su1_cut_logical]
  GI_1_rob <- su1_GI[su1_cut_logical]

  # Load Analysis that is not robust
  load(file2)
  su2 <- su/nrun
  # list elemt 1 is sorted , elemt 2 are the idexes
  sort_index2 <- sort(su2, decreasing = TRUE, index.return=TRUE)
  GID <- as.character(1:500)
  GIsort2 <- GID[sort_index2$ix]
  su1_cut_logical2 <- sort_index2$x>=cut_value
  s2_passingcut <- sort_index2$x[su1_cut_logical2]
  GI_2 <- GIsort2[su1_cut_logical2]

  # Analyse Differences between robust and normal analysis
  G_identical <- intersect(GI_1_rob, GI_2)
  G_different <- c(setdiff(GI_1_rob, GI_2), setdiff(GI_2, GI_1_rob)) # Achtung: check this quick fix
  print(paste("shared genes: ", length(G_identical)))
  print(paste("differnt genes: ", length(G_different)))

  # save genes above cutoff for robust
  Gt <- t(GI_1_rob)
  ch <- as.character(Gt)
  # Achtung: output contains x as header?
  # Achtung: Hier muss output file dynamisch benannt werden / an unterschiedlichem ort gespeichert werden können
  write.table(ch, file = paste("../mcmc_profiles/analyse_GIDout", cut_value, "_rob.txt", sep = ""), sep = "\n")

  # All genes with p-values
  ch <- as.character(t(su1_GI))
  # Achtung wieso ist su aufeinamal p-value?
  # unfertig
  # To Do:
  # Genes für norm
  # p-werte für rob & norm

}

# analyse_GIDout(robust_in, norm_in, numax_in, cut_value_in, nrun_in)
