# Sensitivity Analysis

# Set working directory
setwd("C:/Users/nrios4/OneDrive - George Mason University - O365 Production/Papers in Progress/OnlineOofA")

# Source code
source("OnlineOofASource.R")

# libraries
library(dplyr)

###### Simulation Parameters ################################################
set.seed(123) # random seed 
true_model = "cp" # sensitivity analysis fixed at cp
sigma = 0.5 # sensitivity analysis fixed at sigma = 0.5
m = 8 # sensitivity analysis fixed at m = 8
tol_vec = c(0.25,0.5,0.75) # tolerances for proposed method
n_check_vec = c(3,5,7) # n_checks for proposed method
learning_rate_vec = c(0.05,0.10,0.15) # learning rates for proposed method
if(m == 6){
  n1 = 24 
}
if(m == 7){
  n1 = 36
}
if(m == 8){
  n1 = 48
}
n_sim = 100 # number of times simulation is repeated

#############################################################################

namestr = paste("Sec3SimRes_",true_model,"_m_",m,"_n1_",n1,"_sigma_",sigma, ".csv",sep="")

# this is a quick helper function that just gives the true optimal orders each m
# to avoid having to re-compute them each time
lookup_optimal_order <- function(true_model, m){
  if(true_model == "pwo" || true_model == "lcpf" || true_model == "cp"){
    true_opt_order <- m:1
  }
  if(true_model == "lcpf+qcpf"){
    if(m == 10){
      true_opt_order = c(9,3,5,7,8,6,4,2,1,10) # found using complete enumeration in C
    }
    if(m == 9){
      true_opt_order = c(8,3,5,7,6,4,2,1,9)
    }
    if(m == 8){
      true_opt_order = c(7,3,5,6,4,2,1,8)
    }
    if(m == 7){
      true_opt_order = c(6,2,4,5,3,1,7)
    }
    if(m == 6){
      true_opt_order = c(5,2,4,3,1,6)
    }
  }
  return(true_opt_order)
}

true_opt_order = lookup_optimal_order(true_model,m)
true_opt_y = generate_response(true_opt_order, true_model = true_model,
                               sigma = 0)




sensitivity_param_grid = expand.grid(n_check_vec, tol_vec, learning_rate_vec)
colnames(sensitivity_param_grid) = c("n_c","delta","eta")
sensitivity_param_grid = sensitivity_param_grid[,3:1]

tableB2means = matrix(NA, nrow = nrow(sensitivity_param_grid), ncol = 3)
tableB2sds = matrix(NA, nrow = nrow(sensitivity_param_grid), ncol = 3)

for(row in 1:nrow(sensitivity_param_grid)){
  cat("Starting row : ", row, "\n")
  learning_rate = sensitivity_param_grid$eta[row]
  tol = sensitivity_param_grid$delta[row]
  n_check = sensitivity_param_grid$n_c[row]
  
  
  proposed_sample_sizes = numeric(n_sim)
  proposed_hamming_dist = numeric(n_sim)
  proposed_rel_optimality_gap = numeric(n_sim)
  for(i in 1:n_sim){
    
    cat("Beginning simulation: ", i, "\n")
    
    proposed_out = online_oofa_sampling(n1 = n1, n_max = 1000, m = m, 
                                        true_model = true_model, sigma = sigma,
                                        tol = tol, learning_rate = learning_rate,
                                        models = c("pwo","lcpf","lcpf+qcpf"),
                                        n_check = n_check)
    
    proposed_sample_sizes[i] = proposed_out$sample_size
    proposed_hamming_dist[i] = hamming_distance(proposed_out$best_order, 
                                                true_opt_order)
    proposed_rel_optimality_gap[i] = abs((proposed_out$best_y - true_opt_y)/true_opt_y)
    
    
    
  }
  
  mean_sample_size = round(mean(proposed_sample_sizes),1)
  sd_sample_size = round(sd(proposed_sample_sizes),2)
  mean_hamming = round(mean(proposed_hamming_dist),2)
  sd_hamming = round(sd(proposed_hamming_dist),2)
  mean_rel_opt = round(mean(proposed_rel_optimality_gap),2)
  sd_rel_opt = round(sd(proposed_rel_optimality_gap),2)
  
  tableB2means[row,] = c(mean_sample_size, mean_hamming, mean_rel_opt)
  tableB2sds[row,] = c(sd_sample_size, sd_hamming, sd_rel_opt)
  
}


combined_table = matrix(paste(tableB2means, " (", tableB2sds, ")", sep = ""), nrow = nrow(tableB2sds), ncol = ncol(tableB2sds))
tableB2 = cbind(sensitivity_param_grid, combined_table)
colnames(tableB2) = c(colnames(sensitivity_param_grid), c("Sample Size", "h dist", "rel opt gap"))

# Uncomment to save to file 
# write.csv(tableB2, "sensitivity_analysis.csv")

kable(tableB2, format = "latex", vline = "", linesep = "")