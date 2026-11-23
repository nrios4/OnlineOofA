# MakeTable2.R
# This code makes Table 2.

rm(list = ls())

library(dplyr)
library(kableExtra)

sigma_vec = c(0.1,0.5) 
m_vec = c(9,10)
true_models = c("pwo","cp","lcpf","lcpf+qcpf")

mean_table <- matrix(nrow = 16, ncol = 6)
sd_table <- matrix(nrow = 16, ncol = 6)

row_counter <- 1

for(sigma in sigma_vec){
  
  
  for(i in 1:length(true_models)){
    for(m in m_vec){
      if(m == 9){
        n1 = 38 
      }
      if(m == 10){
        n1 = 47
      }
      
      
      true_model = true_models[i]
      namestr = paste("Sec3SimRes_",true_model,"_m_",m,"_n1_",n1,"_sigma_",sigma, ".csv",sep="")
      results_df_tmp = read.csv(file.path("SimRes", namestr))[,-1]
      tmp_table <- results_df_tmp %>% group_by(method) %>% summarize(sample_size_mean = mean(sample_size), sample_size_sd = sd(sample_size),
                                                                     h_dist_mean = mean(hamming_distance), h_dist_sd = sd(hamming_distance),
                                                                     rel_opt_gap_mean = mean(rel_optimality_gap), rel_opt_gap_sd = sd(rel_optimality_gap))
      mean_table[row_counter, 1:2] <- round(tmp_table$sample_size_mean,1)
      mean_table[row_counter, 3:4] <- round(tmp_table$h_dist_mean, 2)
      mean_table[row_counter, 5:6] <- round(tmp_table$rel_opt_gap_mean, 2)
      sd_table[row_counter, 1:2] <- round(tmp_table$sample_size_sd, 2)
      sd_table[row_counter, 3:4] <- round(tmp_table$h_dist_sd, 2)
      sd_table[row_counter, 5:6] <- round(tmp_table$rel_opt_gap_sd, 2)
      
      
      row_counter <- row_counter + 1
      
      
    }
  }
  
}





table2 = data.frame(sigma = c(rep(0.1,8), rep(0.5,8)),
                    true_model = c(rep("pwo",2), rep("cp",2),rep("FP1",2),rep("FP2",2),
                                   rep("pwo",2), rep("cp",2),rep("FP1",2),rep("FP2",2)),
                    m = c(m_vec, m_vec, m_vec,m_vec,m_vec,m_vec,m_vec,m_vec),
                    n1 = c(38,47,38,47,38,47,38,47,38,47,38,47,38,47,38,47)
)

combined_table = matrix(paste(mean_table, " (", sd_table, ")", sep = ""), nrow = nrow(mean_table), ncol = ncol(mean_table))

table2 = cbind(table2, combined_table)
table2
colnames(table2) = c("sigma", "true_model", "m", "n1", "Alg3 Sample Size",
                      "QS Sample Size", "Alg3 h dist", "QS h dist",
                      "Alg 3 rel opt gap", "QS rel opt gap")
table2 
kable(table2, format = "latex", vline = "", linesep = "")
