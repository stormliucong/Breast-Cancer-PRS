# Last updated: 2020-03-30
# Author: Cong Liu
# cox analysis for PRS project

rm(list = ls())
library(glm2)
library(dplyr)
library(data.table)
library(survival)
library(epitools)
library(epiDisplay)
library(survminer)
library(fmsb)

# read merged table.
merged_table = fread('~/2019-PRS-pipeline/workplace/merged_table.txt')
colnames(merged_table)[4] = 'sex'
merged_table = merged_table[demo_is_consist == 1]

# define helper function.

# prs normalization to SD = 1 MEAN = 0 in the control group.
prs_norm = function(dt, prs_col, case_control_col = 'phekb_case_control') {
  prs_mean = mean(dt[get(case_control_col) == 'control', get(prs_col)])
  prs_sd = sd(dt[get(case_control_col) == 'control', get(prs_col)])
  dt[, prs_score := (get(prs_col) - prs_mean) / prs_sd]
  dt[, case_control := factor(get(case_control_col), levels = c('control','case'))]
  return(dt)
}

clean_column = function(dt){
  dt = copy(dt)
  dt[case_control %in% c('case','control')]
  for(col in colnames(dt)){
    if(!is.na(col)){
      if(class(dt[,get(col)]) %in% c('integer')){
        dt[[col]] = as.numeric(dt[[col]])
      }
      if(class(dt[,get(col)]) %in% c('character')){
        dt[[col]] = factor(dt[[col]])
      }
      if(class(dt[,get(col)]) %in% c('factor')){
        if(length(unique(dt[,get(col)])) < 2){
          dt = dt[,!col, with=FALSE]
        }
      }
    }
  }
  return(dt)
}

# clean glm results.
clean_results = function(glm_model = NULL,
                         working_subset,
                         prs_col,
                         ancestry = 'ALL',
                         sex = 'ALL',
                         family_history = 'ALL',
                         site = 'ALL',
                         subtype = 'ALL',
                         density = 'ALL') {
  dt = data.table(NULL)
  dt[, prs_model := prs_col]
  dt[, ancestry := ancestry]
  dt[, sex := sex]
  dt[, family_history := family_history]
  dt[, site := site]
  dt[, subtype := subtype]
  dt[, density := density]
  dt[, case_number := sum(working_subset$case_control == 'case',na.rm = TRUE)]
  dt[, control_number := sum(working_subset$case_control == 'control',na.rm = TRUE)]
  
  or_display = tryCatch({
    logistic.display(glm_model, simplified = TRUE)
  }, warning = function(w) {
    # warning-handler-code
  }, error = function(e) {
    # error-handler-code
  }, finally = {
    NA
  })

  
  if (!is.null(or_display)) {
    s = summary(glm_model)
    dt = cbind(dt, t(s$conf.int[1, c(1,3,4)]) %>% data.table())
  }
  return(dt)
}

run_glm = function(dt){
  glm_model = tryCatch({
    coxph(Surv(y,as.numeric(case_control))~.,data = dt)
  }, warning = function(w) {
    # warning-handler-code
    # WARNINING 1: glm.fit: fitted probabilities numerically 0 or 1 occurred
    # WARNINING 2: simpleWarning: glm.fit: fitted probabilities numerically 0 or 1 occurred
    print(paste0("warninig in running glm:  ",w))
    print(paste0("return NULL model"))
    return(NULL)
  }, error = function(e) {
    # error-handler-code
    # ERROR 1: case_control not found
    print(paste0("error in running glm:  ",e))
    print(paste0("return NULL model"))
    return(NULL)
  }, finally = {
    NA
  })
}


# primary analysis
primary_subset_col = c(1, 4, 5, 17:27, 28, 32, 33, 50:52, 65:82,40)
prs_names = colnames(merged_table)[c(50:52, 65:82)]
primary_subset = merged_table[, primary_subset_col, with = FALSE]
primary_subset = primary_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
primary_subset = primary_subset[sex %in% c('female')]
primary_analysis_results = NULL
for (p in c('ukbb-f3','bcac-s-f3','bcac-l-f3','root-f3','latinas-f3','whi-aa-f3','whi-la-f3')) {
  for (a in c('eu','aa','la')) {
    for (s in c('female')) {
      working_subset = copy(primary_subset)
      working_subset = working_subset[sex == s & emerge_ancestry == a]
      working_subset = prs_norm(working_subset, prs_col = p, case_control_col = 'phekb_case_control')
      working_subset = working_subset[,.(case_control,prs_score,ancestry_specific_pc1,ancestry_specific_pc2,
                                         ancestry_specific_pc3,phekb_age,site,phekb_family_history,age_at_first_icd)]
      working_subset = clean_column(working_subset)
      working_subset[case_control == 'case',y:=age_at_first_icd][case_control != 'case',y:=phekb_age]
      working_subset = working_subset[!is.na(case_control)]
      glm_model = run_glm(working_subset)
      dt = clean_results(
        glm_model = glm_model,
        working_subset = working_subset,
        prs_col = p,
        sex = s,
        ancestry = a
      )
      primary_analysis_results = rbindlist(list(dt, primary_analysis_results), fill = TRUE)
    }
  }
}
primary_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/cox_analysis.csv')
