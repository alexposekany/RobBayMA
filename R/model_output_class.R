################################################################################
# Structurizes and stores output current iteration as S4 class 

# updated variables
# m_nu: degree of freedom (nu) (updated in: update_mhnuphi.cpp / update_rjnufine.cpp)
# m_phi_scaling: scaling factor phi connecting t-distribution and NV (updated in:
#   update_mhnuphi.cpp / update_rjnufine.cpp)
# m_taue: precision parameter tau for error e (taue|g,h~Ga(g,h), taueup.cpp)
# m_Ig_vec: vector of indicator variable Ig decoding if gene is 
#   differntially expressed (1) or not (0), (update_betai.cpp)
# m_p: probability of gene being differntially expressed, (p ~ Be(a,b), stmcmclm.R)
# m_lambda: prior precision for beta_g, (lambda ~ Ga(c,d), update_lambda.cpp)
# m_tau: TODO: unclear what this is, difference btw tau and taue
# m_Betas: ANOVA parameter vector for gene g (update_betai.cpp)

# other variables 
# m_gene_mu_grouped: group-wise gene mean
# m_total_Ig1: sum of all genes where indicator variable Ig is equal to 1 (i.e. DE),
#   calculated based on m_Ig_vec

################################################################################

# model_outputing Class definition 
setClass(Class = "ModelOutput",
         slots = c(
             # updated variables 
             m_nu = "numeric", m_phi_scaling = "matrix", m_taue = "numeric", 
             m_Ig_vec = "numeric", m_p = "numeric", m_tau = "matrix", 
             m_lambda = "numeric", m_Betas = "matrix", 
             
             # other variables 
             m_gene_mu_grouped = "matrix", m_total_Ig1 = "numeric" 
             )
         )

