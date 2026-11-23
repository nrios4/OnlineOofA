# MakeTableB2.R
# This code makes Table B2.

rm(list = ls())

library(dplyr)
library(kableExtra)

tableB2 = read.csv(file.path("SimRes","sensitivity_analysis.csv"))[,-1]

kable(tableB2, format = "latex", vline = "", linesep = "")