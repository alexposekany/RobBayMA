# --------------------------- DESCRIPTION ------------------------------------ #
# stmcmclm_cpp: Main Algorithm calling differnt updating functions 

## Variables passed to stmcmclm_cpp:
# input_data    GENESxEXPERIMENTS-matrix of measurements
# group_vec     N-dim. vector containing the experimental group indicator
# GId           string vector of gene IDs 
# filename      string of start of filename '.RData' wll be added to save
#               results during run-time for easier recovery
# rob           indicator variable, if a robust (1) or Gaussian analysis (0) 

## prior parameters passed to stmcmclm_cpp:
# hyper_a, hyper_b       hyperparameters of beta prior for p
# hyper_c, hyper_d       hyperparameters of gamma prior for lambda
# hyper_g, hyper_h       hyperparameters of gamma prior for taue


## MCMC related parameters passed to stmcmclm_cpp: 
# numax     max. value of nu before the approximation to normal distribution is considered to be sufficient
# nburn     number of draws till burn-in
# ndraw     total number of draws (after burn-in)


# ---------------------------------------------------------------------------- #

stmcmclm_cpp <- function(input_data, group_vec, GId, hyper_a = 0.5, hyper_b = 0.5, hyper_c = 0, hyper_d  = 0, hyper_g = 0, hyper_h = 0, numax, nburn, ndraw, filename, rob){

    # build input structure
    n_groups <- nlevels(group_vec) 
    model_input <- new("ModelInput", m_data = input_data, m_group_vec = group_vec, m_geneID = GId, 
                 m_hyper_a = hyper_a, m_hyper_b = hyper_b, m_hyper_g = hyper_g, m_hyper_h = hyper_h, 
                 m_hyper_c = hyper_c, m_hyper_d = hyper_d, m_numax = numax,
                 m_grouped_data = group_data1(input_data, group_vec), m_total_mu = mean(colMeans(input_data)), # mean is super slow -> use differnt function!
                 m_design_mat = get_design(group_vec, n_groups), m_n_groups = n_groups
    )
  # initialize algorithmic input / output structure   
  nu <- sample(1:6, size = 6, replace = F) # initialize nu with a random number 
  n_genes <- dim(model_input@m_data)[1]
  n_exp <- dim(model_input@m_data)[2]
  
  model_output <- new("ModelOutput")
  if (rob == 1){
      model_output@m_nu <- 10 + nu[1]*5
      # phi|nu ~ Ga(nu/2, nu/2)
      model_output@m_phi_scaling <- matrix(rgamma(n = (n_genes*n_exp), shape = matrix(model_output@m_nu/2, n_genes, n_exp),
                                                  scale = matrix(2/model_output@m_nu, n_genes, n_exp)), n_genes, n_exp)
     
  } else {
      model_output@m_nu <- numax
      model_output@m_phi_scaling <- matrix(1, n_genes, n_exp)
  }
  
  # taue: precision parameter tau for error e 
  model_output@m_taue=rgamma(1, shape = 0.5, scale = 1)  
  
  # n_groups x n_genes matrix with normally distributed data 
  B=matrix(rnorm(n_groups*n_genes), n_groups, n_genes)
  # Mean of experimental data sorted by experimental group used to initialize model_output@m_gene_mu_grouped
  ms = matrix(0, n_genes, n_groups)
  for(i in 1:n_groups){
      # calculate means of each row of data split by group i.e. mean of all replicates 
      ms[,i] <- rowMeans(model_input@m_grouped_data[[i]])
  }
  
  # Mean of experimental data sorted by experimental group
  model_output@m_gene_mu_grouped = t(ms)
  
  # Vektor der alle Zustände der Indikator Variable speichert. Bei initialisierung wird vector mit 0 & 1 aufgefüllt 
  model_output@m_Ig_vec = sample(x = n_genes, size = n_genes, replace = F) %% 2 
  
  # p
  model_output@m_p = rbeta(1, shape1 = model_input@m_hyper_a, shape2 = model_input@m_hyper_b)
  
  model_output@m_tau = matrix(1, n_groups, n_genes)
  
  # lambda:prior precision for beta_g, lambda ~ Ga(c,d) 
  model_output@m_lambda = rgamma(1, shape = model_input@m_hyper_c, scale = model_input@m_hyper_d)  
  # Betas/beta_g: ANOVA parameter vector for beta_g, Distribution dependent on indicator variable Ig 
  # n_groups x n_genes matrix with normally distributed data 
  model_output@m_Betas = B / model_output@m_tau + model_output@m_gene_mu_grouped  # element wise matrix division
  
  # Count of all genes that are differntially expressed in one run 
  model_output@m_total_Ig1 = round(n_genes/2)
  
  # Initialisierung von McmcOutput
  # su: genewise count of Ig=1 for every run  
  su = matrix(0, 1, n_genes) 
  # i1: Sum of all genes that are differntially expressed in a run 
  i1 = matrix(0, 2*nburn+ndraw, 1)
  nus = matrix(0, 2*nburn+ndraw, 1) 
  taues = matrix(0, 2*nburn+ndraw, 1)
  lambdas = matrix(0, 2*nburn+ndraw, 1)
  ps = matrix(0, 2*nburn+ndraw, 1) # matrix aller p-Werte
  b0 = matrix(0, 1, n_genes) # b0: beta aller gene mit Indikator Ig==0
  bg = matrix(0, 2, n_genes) # bg: Betas aller Gene die Ig == 1
  
  # Initialisierung von Speicherstruktur 
  mcmc_output <- new("McmcOutput")
  mcmc_output@nus <- model_output@m_nu
  mcmc_output@taues <- model_output@m_taue
  mcmc_output@ps <- model_output@m_p
  mcmc_output@Igs <- list(model_output@m_Ig_vec)
  mcmc_output@sus <- su
  
  
  ### Burn-in phase 1
  for( current_run in 1:nburn ) {
      if(rob==1) {
          model_output <- update_mhnuphi_cpp(model_input, model_output)
      } else {
          model_output@m_total_Ig1 <- as.numeric(model_output@m_Ig_vec %*% matrix(1, n_genes, 1))
      }
      
      
      model_output@m_p <- rbeta(1, shape1 = model_input@m_hyper_a + model_output@m_total_Ig1, shape2 = model_input@m_hyper_b + (n_genes - model_output@m_total_Ig1))
      model_output <- update_betai_cpp(model_input, model_output)
      model_output@m_taue <- update_taue_cpp(model_input, model_output)
      
      # save current run
      i1[current_run] <- model_output@m_total_Ig1
      nus[current_run] <- model_output@m_nu
      taues[current_run] <- model_output@m_taue
      lambdas[current_run] <- model_output@m_lambda
      # su <- su + t(model_output@m_Ig_vec) # save su in burn-in or not ACHTUNG variable
      b0[which(model_output@m_Ig_vec == 0)] <- b0[which(model_output@m_Ig_vec == 0)] + model_output@m_Betas[which(model_output@m_Ig_vec == 0)]
      bg[,which(model_output@m_Ig_vec == 1)] <- bg[which(model_output@m_Ig_vec == 1)] + model_output@m_Betas[which(model_output@m_Ig_vec == 1)]
      mcmc_output <- save_iteration_cpp( mcmc_output, model_output, su) # taue, p, I, lambda, ,  nrun = current_run
  }
  
  for( current_run in 1:nburn ) {
      if(rob==1) {
          model_output <- update_rjnufine_cpp(model_input, model_output)
      } else {
          model_output@m_total_Ig1 <- as.numeric(model_output@m_Ig_vec %*% matrix(1, n_genes, 1))
      }
      
      model_output@m_p <- rbeta(1, shape1 = model_input@m_hyper_a + model_output@m_total_Ig1, shape2 = model_input@m_hyper_b + (n_genes - model_output@m_total_Ig1))
      model_output <- update_betai_cpp(model_input, model_output)
      model_output@m_taue <- update_taue_cpp(model_input, model_output)
      
      # save current run
      i1[current_run+nburn] <- model_output@m_total_Ig1
      nus[current_run+nburn] <- model_output@m_nu
      taues[current_run+nburn] <- model_output@m_taue
      lambdas[current_run+nburn] <- model_output@m_lambda
      # su <- su + t(model_output@m_Ig_vec) # save su in burn-in or not ACHTUNG variable
      ps[current_run+nburn] <- model_output@m_p
      b0[which(model_output@m_Ig_vec == 0)] <- b0[which(model_output@m_Ig_vec == 0)] + model_output@m_Betas[which(model_output@m_Ig_vec == 0)]
      bg[,which(model_output@m_Ig_vec == 1)] <- bg[which(model_output@m_Ig_vec == 1)] + model_output@m_Betas[which(model_output@m_Ig_vec == 1)]
      mcmc_output <- save_iteration_cpp( mcmc_output, model_output, su) # taue, p, I, lambda, ,  nrun = current_run
  }
  
  #####SAMPLING START
  for(current_run in 1:ndraw){
      if(rob==1){
          model_output<-update_rjnufine_cpp(model_input,model_output)
      } else {
          model_output@m_total_Ig1 = as.numeric(model_output@m_Ig_vec %*% matrix(1, n_genes, 1))
      }
      
      model_output@m_p <- rbeta(1, shape1 = model_input@m_hyper_a + model_output@m_total_Ig1, shape2 = model_input@m_hyper_b + (n_genes - model_output@m_total_Ig1))
      model_output <- update_betai_cpp(model_input, model_output)
      model_output@m_taue <- update_taue_cpp(model_input, model_output)
      
      i1[current_run+2*nburn] <- model_output@m_total_Ig1
      nus[current_run+2*nburn] <- model_output@m_nu
      taues[current_run+2*nburn] <- model_output@m_taue
      lambdas[current_run+2*nburn] <- model_output@m_lambda
      su <- su + t(model_output@m_Ig_vec)
      ps[current_run+2*nburn] <- model_output@m_p
      b0[which(model_output@m_Ig_vec == 0)] <- b0[which(model_output@m_Ig_vec == 0)] + model_output@m_Betas[which(model_output@m_Ig_vec == 0)]
      bg[,which(model_output@m_Ig_vec == 1)] <- bg[which(model_output@m_Ig_vec == 1)] + model_output@m_Betas[which(model_output@m_Ig_vec == 1)]
      
      mcmc_output <- save_iteration_cpp( mcmc_output, model_output, su) # taue, p, I, lambda, ,  nrun = current_run
      
      # save to file after every 250ten sampling run
      if( (current_run+1)%%250 == 0) {
          save(su, nus, taues, i1, lambdas, ps, b0, file = paste(filename, '.RData', sep = ''))
          save(bg, file = paste(filename, '_beta', '.RData', sep = ''))
          save(model_output, file = paste(filename, '_str', '.RData', sep = ''))
      }
  }
    #####SAMPLING STOP
    return(mcmc_output)
    }

