# Section 4

rm(list = ls()) # clear envrionment

set.seed(1234)

source("OnlineOofASource.R")

m = 15
# generate processing times
p_vec = round(rchisq(n = m, df = 1, ncp = 0),4)
print(p_vec)
# generate weights
omega_vec = round(rchisq(n = m, df = 1, ncp = 0),4)
print(omega_vec)


# Example parameters
m <- 15
n1 <- 37
n_sim = 200
learning_rate = 0.1; tol = 0.5; n_check = 7
sigma = 0.5

namestr = paste("initial_design_m15_job_scheduling_n1_",n1, ".csv", sep = "")

if(file.exists(namestr)){
  initial_design = read.csv(namestr)[,-1]
  colnames(initial_design) = c()
} else{
  # Find initial design for Phase I
  init_runs <- t(replicate( n1, sample(seq_len(m))))
  # choose(m,2) = 105
  # Partition PWO effects into 3 non-overlapping groups of 35
  pwo1 = seq(from = 1, to = 105, by = 3)
  pwo2 = seq(from = 2, to = 105, by = 3)
  pwo3 = seq(from = 3, to = 105, by = 3)
  
  # Find initial design under PWO1
  opt_design <- find_opt_design_SA_reduced(
    init_runs = init_runs,
    pwo_cols = pwo1,
    n_iter = 2000
  )
  
  # Check that it is full rank under all 3 reduced PWO models
  Z_opt_design = gen.pf(opt_design, "pwof")
  X1 = Z_opt_design[,pwo1]
  X2 = Z_opt_design[,pwo2]
  X3 = Z_opt_design[,pwo3]
  qr(X1)$rank == ncol(X1)
  qr(X2)$rank == ncol(X2)
  qr(X3)$rank == ncol(X3)
  # All of the above should be TRUE. 
  # If not, re-run under different seed/initial design.
  
  # Save the initial design to save time in future runs
  write.csv(opt_design, namestr)
  initial_design = opt_design
}

initial_design = as.matrix(initial_design)

quicksort_costs = numeric(n_sim)
quicksort_sample_sizes = numeric(n_sim)
proposed_costs = numeric(n_sim)
proposed_sample_sizes = numeric(n_sim)

for(i in 1:n_sim){
  cat("Beginning Simulation:", i, "\n")
  
  
  quicksort_out = quicksort_algorithm_Section4(n1 = n1, m = m, 
                                               p_vec = p_vec, omega_vec = omega_vec, sigma = sigma)
  quicksort_costs[i] = quicksort_out$best_y
  quicksort_sample_sizes[i] = quicksort_out$sample_size
  
  proposed_out = online_oofa_sampling_Section4(initial_design = initial_design,
                                               n_max = 1000, m = m, p_vec = p_vec,
                                               omega_vec = omega_vec, sigma = sigma,
                                               tol = tol, learning_rate = learning_rate,
                                               n_check = n_check,
                                               models = c("pwo1","pwo2","pwo3","lcpf","lcpf+qcpf"))
  
  proposed_costs[i] = proposed_out$best_y
  proposed_sample_sizes[i] = proposed_out$sample_size
  
}


### get baseline for comparison, use 5,000 randomly selected orders
random_orders = t(replicate(5000, sample(seq_len(m))))
random_order_costs = apply(random_orders, 1, get_order_cost, p_vec = p_vec,
                           omega_vec = omega_vec, sigma = 0)
summary(random_order_costs)

res_df = data.frame( method = c(rep("quicksort",n_sim), rep("proposed",n_sim),
                                rep("random", length(random_order_costs))),
                                sample_size = c(quicksort_sample_sizes, proposed_sample_sizes,
                                                rep(NA, length(random_order_costs))),
                                cost = c(quicksort_costs, proposed_costs,
                                                         random_order_costs))
  
# uncomment to save the results                              
# write.csv(res_df, "Section4Results_nc7.csv")

par(mfrow = c(1,2))
boxplot(sample_size ~ method, data = res_df, ylab = "Sample Size")
boxplot(cost ~ method, data = res_df, ylab = "True Cost")

library(dplyr)
res_df %>% group_by(method) %>% summarize(med_size = median(sample_size, na.rm = T),
                                          IQR_size = IQR(sample_size, na.rm = T),
                                          med_cost = median(cost, na.rm = T),
                                          IQR_cost = IQR(cost, na.rm = T))
