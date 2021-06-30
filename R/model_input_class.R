################################################################################
# Structurizes data input as S4 class which is used to pass data to all subfunctions.

# Variables used:
# m_data:           Data matrix with GxN dimensions containing measurements for y_g,n
# m_group_vec:      N-dim. vector with biological class indicators s_n
# m_grouped_data:   data Y grouped according to grouping vector s
# m_total_mu:       mean of betas, the mean of data means 
# m_design_mat:     matrix containing design vectors
# m_geneId:         string vector of gene IDs 
# m_numax:          max. nu value before it can be approximated with normal distribution
# m_n_groups:       number of groups 

# Hyperparameters:
# m_hyper_a, m_hyper_b: hyperparameters for p-DE, following beta distribution
# m_hyper_g, m_hyper_h: hyperparameters for error taue, following gamma distribution 
# m_hyper_c, m_hyper_d: hyperparameters for precision parameter lambda, following gamma distribution 

################################################################################

setClass(Class = "ModelInput",
         slots = c(
             # Variables
             m_data = "matrix", m_group_vec = "factor", 
             m_grouped_data = "list", # lists contains matrices 
             m_design_mat = "matrix", m_total_mu = "numeric", m_geneID = "character",
             m_numax = "numeric", m_n_groups = "numeric",
             
             # Hyperparameters
             m_hyper_a = "numeric", m_hyper_b = "numeric", m_hyper_g = "numeric",
             m_hyper_h = "numeric", m_hyper_c = "numeric", m_hyper_d = "numeric"
      
      
  )
  
)

  