# Last updated: 2020-05-07
# Author: Cong Liu
# logistic regression analysis with interaction terms for PRS project

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
    dt = cbind(dt, t(coef(s)[2, ]) %>% data.table())
    dt[, prs_score_OR := or_display$table[1, 1]]
    dt[, prs_score_OR_lower95ci := or_display$table[1, 2]]
    dt[, prs_score_OR_upper95ci := or_display$table[1, 3]]
    dt[, prs_score_pvalue := coef(s)[2,4]]
    dt[, strata_OR := or_display$table[2, 1]]
    dt[, strata_OR_lower95ci := or_display$table[2, 2]]
    dt[, strata_OR_upper95ci := or_display$table[2, 3]]
    dt[, strata_OR_pvalue := coef(s)[3,4]]
    col_number = dim(coef(s))[1]
    dt[, interaction_OR := or_display$table[col_number - 1, 1]]
    dt[, interaction_OR_lower95ci := or_display$table[col_number - 1, 2]]
    dt[, interaction_OR_upper95ci := or_display$table[col_number - 1, 3]]
    dt[, interaction_pvalue := coef(s)[col_number,4]]
  }
  return(dt)
}

clean_lrt = function(lrt = NULL,
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
  
  
  
  if (!is.null(lrt$`Pr(>Chi)`)) {
    dt[, Df := lrt$Df[2]]
    dt[, Deviance := lrt$Deviance[2]]
    dt[, pvalue := lrt$`Pr(>Chi)`[2]]
  }
  return(dt)
}

run_glm = function(formula,dt){
  glm_model = tryCatch({
    glm(
      formula,
      family = binomial,
      data = dt
    )
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

# LRT for family history * PRS
interaction_subset_col = c(1, 4, 5, 17:27, 28, 32, 33, 50:52, 65:82,41,42)
prs_names = colnames(merged_table)[c(50:52, 65:82)]
interaction_subset = merged_table[, interaction_subset_col, with = FALSE]
interaction_subset = interaction_subset[phekb_case_control %in% c('case', 'control')]
interaction_subset = interaction_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
interaction_subset = interaction_subset[sex %in% 'female']
interaction_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu','aa','la')) {
    for (s in 'female') {
      working_subset = copy(interaction_subset)
      working_subset = working_subset[sex == s & emerge_ancestry == a]
      working_subset = prs_norm(working_subset, prs_col = p, case_control_col = 'phekb_case_control')
      working_subset = working_subset[,.(case_control,prs_score,ancestry_specific_pc1,ancestry_specific_pc2,
                                         ancestry_specific_pc3,phekb_age,site,phekb_family_history)]
      # working_subset = working_subset[,phekb_family_history := as.numeric(factor(phekb_family_history,levels = c('wo_fx','w_fx')))]
      working_subset = clean_column(working_subset)
      formula1 = as.formula("case_control ~ prs_score*phekb_family_history + ancestry_specific_pc1 + ancestry_specific_pc2 + ancestry_specific_pc3 + phekb_age")
      formula2 = as.formula("case_control ~ prs_score + phekb_family_history + ancestry_specific_pc1 + ancestry_specific_pc2 + ancestry_specific_pc3 + phekb_age")
      
      glm_model_1 = run_glm(formula = formula1, dt = working_subset)
      glm_model_2 = run_glm(formula = formula2, dt = working_subset)
      lrt = anova(glm_model_2, glm_model_1, test="LRT")
      
      dt = clean_lrt(
        lrt = lrt,
        working_subset = working_subset,
        prs_col = p,
        sex = s,
        ancestry = a
      )
      interaction_analysis_results = rbindlist(list(dt, interaction_analysis_results), fill = TRUE)
    }
  }
}
interaction_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/interaction_analysis_results_fx.csv')

# LRT for breast density * PRS
interaction_subset_col = c(1, 4, 5, 17:27, 28, 32, 33, 38, 50:52, 65:82)
prs_names = colnames(merged_table)[c(50:52, 65:82)]
interaction_subset = merged_table[, interaction_subset_col, with = FALSE]
interaction_subset = interaction_subset[phekb_case_control %in% c('case', 'control')]
interaction_subset = interaction_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
interaction_subset = interaction_subset[sex %in% 'female']
interaction_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu','aa','la')) {
    for (s in 'female') {
      working_subset = copy(interaction_subset)
      working_subset = working_subset[sex == s & emerge_ancestry == a]
      working_subset = prs_norm(working_subset, prs_col = p, case_control_col = 'phekb_case_control')
      working_subset = working_subset[,.(case_control,prs_score,ancestry_specific_pc1,ancestry_specific_pc2,
                                         ancestry_specific_pc3,phekb_age,site,breast_density)]
      # working_subset = working_subset[,breast_density := as.numeric(factor(breast_density,levels = c('entirely_fat','scattered_fibro_dens','hetero_dense','extremely_dense')))]
      working_subset = clean_column(working_subset)
      formula1 = as.formula("case_control ~ prs_score*breast_density + ancestry_specific_pc1 + ancestry_specific_pc2 + ancestry_specific_pc3 + phekb_age")
      formula2 = as.formula("case_control ~ prs_score + breast_density + ancestry_specific_pc1 + ancestry_specific_pc2 + ancestry_specific_pc3 + phekb_age")
      glm_model_1 = run_glm(formula = formula1, dt = working_subset)
      glm_model_2 = run_glm(formula = formula2, dt = working_subset)
      lrt = anova(glm_model_2, glm_model_1, test="LRT")
      dt = clean_lrt(
        lrt = lrt,
        working_subset = working_subset,
        prs_col = p,
        sex = s,
        ancestry = a
      )
      interaction_analysis_results = rbindlist(list(dt, interaction_analysis_results), fill = TRUE)
    }
  }
}
interaction_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/interaction_analysis_results_density.csv')

# LRT analysis for site.
interaction_subset_col = c(1, 4, 5, 17:27, 28, 32, 33, 38, 50:52, 65:82)
prs_names = colnames(merged_table)[c(50:52, 65:82)]
interaction_subset = merged_table[, interaction_subset_col, with = FALSE]
interaction_subset = interaction_subset[phekb_case_control %in% c('case', 'control')]
interaction_subset = interaction_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
interaction_subset = interaction_subset[sex %in% 'female']
interaction_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu','aa','la')) {
    for (s in 'female') {
      working_subset = copy(interaction_subset)
      working_subset = working_subset[sex == s & emerge_ancestry == a]
      working_subset = prs_norm(working_subset, prs_col = p, case_control_col = 'phekb_case_control')
      working_subset = working_subset[,.(case_control,prs_score,ancestry_specific_pc1,ancestry_specific_pc2,
                                         ancestry_specific_pc3,phekb_age,site,phekb_family_history)]
      # working_subset = working_subset[,breast_density := as.numeric(factor(breast_density,levels = c('entirely_fat','scattered_fibro_dens','hetero_dense','extremely_dense')))]
      working_subset = clean_column(working_subset)
      formula1 = as.formula("case_control ~ prs_score * site + ancestry_specific_pc1 + ancestry_specific_pc2 + ancestry_specific_pc3 + phekb_age+phekb_family_history")
      formula2 = as.formula("case_control ~ prs_score + ancestry_specific_pc1 + ancestry_specific_pc2 + ancestry_specific_pc3 + phekb_age+phekb_family_history")
      glm_model_1 = run_glm(formula = formula1, dt = working_subset)
      glm_model_2 = run_glm(formula = formula2, dt = working_subset)
      lrt = anova(glm_model_2, glm_model_1, test="LRT")
      dt = clean_lrt(
        lrt = lrt,
        working_subset = working_subset,
        prs_col = p,
        sex = s,
        ancestry = a
      )
      interaction_analysis_results = rbindlist(list(dt, interaction_analysis_results), fill = TRUE)
    }
  }
}
interaction_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/sensitivity_analysis_results_site.csv')

