################################################################################
# Structurizes and stores algorithm output of all iteration as S4 class 

# nus, taues, ps, Igs, lambdas explanation see model_output_class.R
# sus: counts how often one gene had Ig=1
################################################################################

# define class 
setClass(Class = "McmcOutput",
         slots = c(nus = "numeric", taues = "numeric", ps = "numeric",
                   Igs = "list", sus = "matrix", lambdas = "numeric")
         )

# setClass(Class = "McmcOutput2",
#          slots = c(nus = "numeric", taues = "numeric", ps = "numeric",
#                    Igs = "list", sus = "matrix", lambdas = "numeric")
# )

# reserve name of method with setGeneric
setGeneric(name = "saving_current_run", 
           def = function(mcmcstr_object, model_output_object, su, nrun) { #taue, p, I, lambda, 
             standardGeneric("saving_current_run")
             }
           )

setMethod(f = "saving_current_run", 
          signature = c("McmcOutput", "ModelOutput"), 
          definition = function(mcmcstr_object, model_output_object, su, nrun){ 
              # browser()
            mcmcstr_object@nus <- append(mcmcstr_object@nus, model_output_object@m_nu)
            mcmcstr_object@taues <- append(mcmcstr_object@taues, model_output_object@m_taue)
            mcmcstr_object@ps <- append(mcmcstr_object@ps, model_output_object@m_p)
            # mcmcstr_object@Igs[nrun,] <- model_output_object@m_Ig_vec
            mcmcstr_object@Igs[nrun] <- list(model_output_object@m_Ig_vec)
            mcmcstr_object@lambdas <- append(mcmcstr_object@lambdas, model_output_object@m_lambda)
            mcmcstr_object@sus <- su
            
            return(mcmcstr_object)
            }
          )


