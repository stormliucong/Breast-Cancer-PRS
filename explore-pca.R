# Last updated: 2020-03-23
# Author: Cong Liu
# explore how many PCs should be included in each population.

rm(list = ls())
library(glm2)
library(dplyr)
library(data.table)
library(survival)
library(epitools)
library(epiDisplay)
library(survminer)

# read merged table.
merged_table = fread('~/2019-PRS-pipeline/workplace/merged_table.txt')
colnames(merged_table)

# define helper function.

# prs normalization to SD = 1 MEAN = 0 in the control group.
prs_norm = function(dt, prs_col, case_control_col = 'phekb_case_control') {
  prs_mean = mean(dt[get(case_control_col) == 'control', get(prs_col)])
  prs_sd = sd(dt[get(case_control_col) == 'control', get(prs_col)])
  dt[, prs_score := (get(prs_col) - prs_mean) / prs_sd]
  dt[, case_control := factor(get(case_control_col), levels = c('control','case'))]
  return(dt)
}


primary_subset_col = c(1, 4, 5, 17:27, 28, 32, 33, 50:52, 65:82)
prs_names = colnames(merged_table)[c(50:52, 65:82)]
primary_subset = merged_table[, primary_subset_col, with = FALSE]
primary_subset = primary_subset[phekb_case_control %in% c('case', 'control')]
primary_subset = primary_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
primary_subset = primary_subset[sex %in% c('male', 'female')]
primary_analysis_results = NULL
primary_analysis_results_quantile = NULL

for (a in c('eu', 'aa', 'la')) {
  cat(paste0('ancestry : ', a,'\n'))
  working_subset = copy(primary_subset)
  working_subset = working_subset[emerge_ancestry == a]
  working_subset = prs_norm(working_subset, prs_col = p, case_control_col = 'phekb_case_control')
  working_subset = working_subset[,.(case_control,ancestry_specific_pc1,ancestry_specific_pc2,
                                     ancestry_specific_pc3,ancestry_specific_pc4,ancestry_specific_pc5,
                                     ancestry_specific_pc6,ancestry_specific_pc7,ancestry_specific_pc8,
                                     ancestry_specific_pc9,ancestry_specific_pc10)]
  working_subset = clean_column(working_subset)
  glm_model = glm(
    case_control ~ .,
    family = binomial,
    data = working_subset
  )
  
  print(summary(glm_model)$coefficients)
}
