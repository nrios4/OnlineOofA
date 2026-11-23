# MakeTables1andB1.R
# This code makes Table 1 and Table B1.

rm(list = ls())

library(dplyr)
library(kableExtra)


sigma = 0.1 # 0.5 for Table 1, 0.1 for Table B1
m_vec = c(6,7,8)
true_models = c("pwo","cp","lcpf","lcpf+qcpf")

mean_table <- matrix(nrow = 12, ncol = 6)
sd_table <- matrix(nrow = 12, ncol = 6)

row_counter <- 1
for(i in 1:length(true_models)){
  for(m in m_vec){
    if(m == 6){
      n1 = 24 
    }
    if(m == 7){
      n1 = 36
    }
    if(m == 8){
      n1 = 48
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



table1 = data.frame(true_model = c(rep("pwo",3), rep("cp",3),rep("FP1",3),rep("FP2",3)),
                        m = c(m_vec, m_vec, m_vec,m_vec),
                        n1 = c(24,36,48,24,36,48,24,36,48,24,36,48)
                        )

combined_table = matrix(paste(mean_table, " (", sd_table, ")", sep = ""), nrow = nrow(mean_table), ncol = ncol(mean_table))

table1 = cbind(table1, combined_table)
colnames(table1) = c("true_model", "m", "n1", "Alg3 Sample Size",
                     "QS Sample Size", "Alg3 h dist", "QS h dist",
                     "Alg 3 rel opt gap", "QS rel opt gap")
table1 # again, if sigma = 0.1, this is Table B1
kable(table1, format = "latex", vline = "", linesep = "")
