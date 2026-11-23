# Section 3.2 Simulation

# Source code
source("OnlineOofASource.R")

# libraries
library(dplyr)

###### Simulation Parameters ################################################
set.seed(123) # random seed 
true_model = "lcpf+qcpf" # one of "pwo", "cp", "lcpf" (FP1), "lcpf+qcpf" (FP2)
sigma = 0.5 # 0.1, 0.5
m = 10 # m = 9,10
tol = 0.5 # tolerance for proposed method
n1 = choose(m,2) + 2 # minimum sample size 
n_sim = 100 # number of times simulation is repeated
learning_rate = 0.1 # learning rate for proposed method
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



proposed_sample_sizes = numeric(n_sim)
proposed_hamming_dist = numeric(n_sim)
proposed_rel_optimality_gap = numeric(n_sim)
quicksort_sample_sizes = numeric(n_sim)
quicksort_hamming_dist = numeric(n_sim)
quicksort_rel_optimality_gap = numeric(n_sim)

for(i in 1:n_sim){
  
  cat("Beginning simulation: ", i, "\n")
  
  proposed_out = online_oofa_sampling(n1 = n1, n_max = 1000, m = m, 
                                      true_model = true_model, sigma = sigma,
                                      tol = tol, learning_rate = learning_rate,
                                      models = c("pwo","lcpf","lcpf+qcpf"))
  
  proposed_sample_sizes[i] = proposed_out$sample_size
  proposed_hamming_dist[i] = hamming_distance(proposed_out$best_order, 
                                              true_opt_order)
  proposed_rel_optimality_gap[i] = abs((proposed_out$best_y - true_opt_y)/true_opt_y)
  
  
  quicksort_out = quicksort_algorithm(n1 = 20, m = m, 
                                      true_model = true_model, sigma = sigma)
  
  quicksort_sample_sizes[i] = quicksort_out$sample_size
  quicksort_hamming_dist[i] = hamming_distance(quicksort_out$best_order, 
                                               true_opt_order)
  quicksort_rel_optimality_gap[i] = abs((quicksort_out$best_y - true_opt_y)/true_opt_y)
  
}

results_df = data.frame(method = c(rep("Proposed", n_sim), rep("Quicksort", n_sim) ),
                        sample_size = c(proposed_sample_sizes, quicksort_sample_sizes),
                        hamming_distance = c(proposed_hamming_dist, quicksort_hamming_dist),
                        rel_optimality_gap = c(proposed_rel_optimality_gap, quicksort_rel_optimality_gap)
)

# uncomment to save the results
# write.csv(results_df,namestr)

results_df %>% group_by(method) %>% summarize(sample_size_mean = mean(sample_size), sample_size_sd = sd(sample_size),
                                              h_dist_mean = mean(hamming_distance), h_dist_sd = sd(hamming_distance),
                                              rel_opt_gap_mean = mean(rel_optimality_gap), rel_opt_gap_sd = sd(rel_optimality_gap))

# uncomment to see boxplots
# boxplot(sample_size ~ method, data = results_df)
# boxplot(hamming_distance ~ method, data = results_df)
# boxplot(rel_optimality_gap ~ method, data = results_df)


