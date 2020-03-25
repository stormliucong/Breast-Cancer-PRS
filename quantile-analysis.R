# Last updated: 2020-03-23
# Author: Cong Liu
# quantile analysis for PRS project


rm(list=ls())
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
colnames(merged_table)

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


# add quantile column into dt.
add_prs_quantile = function(dt,ancestry){
  results = tryCatch({
    dt = copy(dt)
    ref_labels = c('40-60%')
    if (ancestry == 'eu') {
      labels = c(
        '<1%',
        '1-5%',
        '5-10%',
        '10-20%',
        '20-40%',
        '40-60%',
        '60-80%',
        '80-90%',
        '90-95%',
        '95-99%',
        '>=99%'
      )
      dt[, quantile := cut(
        prs_score,
        breaks = quantile(prs_score, probs = c(
          0, 1, 5, 10, 20, 40, 60, 80, 90, 95, 99, 100
        ) / 100),
        labels = labels,
        right = TRUE
      )]
      dt[is.na(quantile), quantile := '<1%']
      
    } else{
      labels = c('<20%', '20-40%', '40-60%', '60-80%', '>=80%')
      dt[, quantile := cut(
        prs_score,
        breaks = quantile(prs_score, probs = c(0, 20, 40, 60, 80, 100) / 100),
        labels = labels,
        right = TRUE
      )]
      dt[is.na(quantile), quantile := '<20%']
    }
  }, warning = function(w) {
    # warning-handler-code
  }, error = function(e) {
    # error-handler-code
  }, finally = {
    NA
  })
  if(!is.null(results)){
    return(dt)
  }else{
    return(NULL)
  }
  
}

# clean glm results.
# clean glm results.
clean_results = function(glm_model = NULL,
                         working_subset_q,
                         quantile,
                         prs_col,
                         ancestry = 'ALL',
                         sex = 'ALL',
                         family_history = 'ALL',
                         site = 'ALL',
                         subtype = 'ALL',
                         density = 'ALL') {
  dt = data.table(NULL)
  dt[, prs_model := prs_col]
  dt[, quantile := quantile]
  dt[, ancestry := ancestry]
  dt[, sex := sex]
  dt[, family_history := family_history]
  dt[, site := site]
  dt[, subtype := subtype]
  dt[, density := density]
  dt[, case_number := sum(working_subset_q$case_control == 'case')]
  dt[, control_number := sum(working_subset_q$case_control == 'control')]
  
  
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
    dt[, OR := or_display$table[1, 1]]
    dt[, OR_lower95ci := or_display$table[1, 2]]
    dt[, OR_upper95ci := or_display$table[1, 3]]
    dt[, logOR := log(OR)]
    dt[, logOR_lower95ci := log(OR_lower95ci)]
    dt[, logOR_upper95ci := log(OR_upper95ci)]
    dt[, rsqr := NagelkerkeR2(glm_model)$R2]
  }
  return(dt)
}

run_glm = function(dt){
  glm_model = tryCatch({
    glm(
      case_control ~ .,
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


# primary analysis
primary_subset_col = c(1,4,5,17:27,28,32,33,50:52,65:82)
prs_names = colnames(merged_table)[c(50:52,65:82)]
primary_subset = merged_table[,primary_subset_col,with = FALSE]
primary_subset = primary_subset[phekb_case_control %in% c('case','control')]
primary_subset = primary_subset[emerge_ancestry %in% c('eu','aa','la')]
primary_subset = primary_subset[sex %in% c('male','female')]
primary_analysis_results = NULL
for(p in prs_names){
  for(a in c('eu','aa','la')){
    for(s in c('male','female')){
      working_subset = copy(primary_subset)
      working_subset = working_subset[sex==s & emerge_ancestry==a]
      working_subset = prs_norm(working_subset,prs_col = p,case_control_col = 'phekb_case_control')
      working_subset = add_prs_quantile(working_subset,ancestry = a)
      working_subset = working_subset[,.(case_control,quantile,ancestry_specific_pc1,ancestry_specific_pc2,
                                         ancestry_specific_pc3,
                                         phekb_age,site,phekb_family_history)]
      for(q in unique(working_subset$quantile)){
        if(q != '40-60%'){
          working_subset_q = working_subset[quantile == q | quantile == '40-60%']
          working_subset_q[,quantile:=factor(quantile,levels = c('40-60%',q))]
          working_subset_q = clean_column(working_subset_q)
          glm_model = run_glm(working_subset_q)
          dt = clean_results(glm_model = glm_model, working_subset_q = working_subset_q, prs_col = p,quantile = q,sex = s,ancestry = a)
          primary_analysis_results = rbindlist(list(dt, primary_analysis_results), fill=TRUE)
        }
      }
    }
  }
}
primary_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/primary_analysis_qt.csv')

# primary analysis - sensitivity analysis
# s1 - results in ICD defined case-control
sensitivity_subset_col = c(1, 4, 5, 17:27, 43, 41, 33, 50:52,65:82)
prs_names = colnames(merged_table)[c(50:52,65:82)]
sensitivity_subset = merged_table[, sensitivity_subset_col, with = FALSE]
sensitivity_subset = sensitivity_subset[icd_case_control %in% c('case', 'control')]
sensitivity_subset = sensitivity_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
sensitivity_subset = sensitivity_subset[sex %in% c('male', 'female')]
sensitivity_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu', 'aa', 'la')) {
    for (s in c('male', 'female')) {
      working_subset = copy(sensitivity_subset)
      working_subset = working_subset[sex == s & emerge_ancestry == a]
      working_subset = prs_norm(working_subset, prs_col = p, case_control_col = 'icd_case_control')
      working_subset = add_prs_quantile(working_subset,ancestry = a)
      working_subset = working_subset[,.(case_control,quantile,ancestry_specific_pc1,ancestry_specific_pc2,
                                         ancestry_specific_pc3,
                                         year_birth,site,phekb_family_history)]
      for(q in unique(working_subset$quantile)){
        if(q != '40-60%'){
          working_subset_q = working_subset[quantile == q | quantile == '40-60%']
          working_subset_q[,quantile:=factor(quantile,levels = c('40-60%',q))]
          working_subset_q = clean_column(working_subset_q)
          glm_model = run_glm(working_subset_q)
          dt = clean_results(glm_model = glm_model,working_subset_q = working_subset_q, prs_col = p,quantile = q,sex = s,ancestry = a)
          sensitivity_analysis_results = rbindlist(list(dt, sensitivity_analysis_results), fill = TRUE)
          
        }
      } 
    }
  }
}
sensitivity_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/sensitivity_analysis_icd_qt.csv')

# s2 - comparing w/ reminder population.
sensitivity_subset_col = c(1, 4, 5, 17:27, 28, 32, 33, 50:52,65:82)
prs_names = colnames(merged_table)[c(50:52,65:82)]
sensitivity_subset = merged_table[, sensitivity_subset_col, with = FALSE]
sensitivity_subset = sensitivity_subset[phekb_case_control %in% c('case', 'control')]
sensitivity_subset = sensitivity_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
sensitivity_subset = sensitivity_subset[sex %in% c('male', 'female')]
sensitivity_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu', 'aa', 'la')) {
    for (s in c('male', 'female')) {
      working_subset = copy(sensitivity_subset)
      working_subset = working_subset[sex == s & emerge_ancestry == a]
      working_subset = prs_norm(working_subset, prs_col = p, case_control_col = 'phekb_case_control')
      working_subset = add_prs_quantile(working_subset,ancestry = a)
      for(q in unique(working_subset$quantile)){
        working_subset_q = working_subset[quantile == q | quantile == '40-60%']
        working_subset_q[quantile != q,quantile_class:= 0]
        working_subset_q[quantile == q,quantile_class:= 1]
        working_subset_q = working_subset_q[,.(case_control,quantile_class,ancestry_specific_pc1,ancestry_specific_pc2,
                                           ancestry_specific_pc3,
                                           phekb_age,site,phekb_family_history)]
        working_subset_q = clean_column(working_subset_q)
        glm_model = run_glm(working_subset_q)
        dt = clean_results(glm_model = glm_model, working_subset_q = working_subset_q, prs_col = p,quantile = q,sex = s,ancestry = a)
        sensitivity_analysis_results = rbindlist(list(dt, sensitivity_analysis_results), fill = TRUE)
        
      } 
    }
  }
}
sensitivity_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/sensitivity_analysis_overall_qt.csv')


# stratified analysis.
# site.
stratify_subset_col = c(1,4,5,17:27,28,32,33,50:52,65:82)
prs_names = colnames(merged_table)[c(50:52,65:82)]
stratify_subset = merged_table[, stratify_subset_col, with = FALSE]
stratify_subset = stratify_subset[phekb_case_control %in% c('case', 'control')]
stratify_subset = stratify_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
stratify_subset = stratify_subset[sex %in% c('male', 'female')]
stratify_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu', 'aa', 'la')) {
    for (s in c('male', 'female')) {
      for (st in c("mrsh",
                   "vand",
                   "kpuw",
                   "colu",
                   "mayo",
                   "nwun",
                   "geis",
                   "harv")) {
        working_subset = copy(stratify_subset)
        working_subset = working_subset[sex == s & emerge_ancestry == a & site == st]
        working_subset = prs_norm(working_subset, prs_col = p, case_control_col = 'phekb_case_control')
        working_subset = add_prs_quantile(working_subset,ancestry = a)
        if(!is.null(working_subset)){
          working_subset = working_subset[,.(case_control,quantile,ancestry_specific_pc1,ancestry_specific_pc2,
                                             ancestry_specific_pc3,
                                             phekb_age,site,phekb_family_history)]
          for(q in unique(working_subset$quantile)){
            if(q != '40-60%'){
              working_subset_q = working_subset[quantile == q | quantile == '40-60%']
              working_subset_q[,quantile:=factor(quantile,levels = c('40-60%',q))]
              working_subset_q = clean_column(working_subset_q)
              glm_model = run_glm(working_subset_q)
              dt = clean_results(glm_model = glm_model, working_subset_q = working_subset_q, prs_col = p,quantile = q,sex = s,ancestry = a, site = st)
              stratify_analysis_results = rbindlist(list(dt, stratify_analysis_results), fill=TRUE)
              
            }
          }
        }
      }
    }
  }
}
stratify_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/stratify_analysis_site_qt.csv')

# family history
stratify_subset_col = c(1,4,5,17:27,28,32,33,50:52,65:82)
prs_names = colnames(merged_table)[c(50:52,65:82)]
stratify_subset = merged_table[, stratify_subset_col, with = FALSE]
stratify_subset = stratify_subset[phekb_case_control %in% c('case', 'control')]
stratify_subset = stratify_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
stratify_subset = stratify_subset[sex %in% c('male', 'female')]
stratify_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu', 'aa', 'la')) {
    for (s in c('male', 'female')) {
      for (st in c("wo_fx", "w_fx")) {
        working_subset = copy(stratify_subset)
        working_subset = working_subset[sex == s &
                                          emerge_ancestry == a
                                        &
                                          ((
                                            phekb_case_control == 'case' & phekb_family_history == st
                                          )
                                          |
                                            phekb_case_control == 'control'
                                          )]
        working_subset = prs_norm(working_subset,
                                  prs_col = p,
                                  case_control_col = 'phekb_case_control')
        working_subset = add_prs_quantile(working_subset,ancestry = a)
        
        if(!is.null(working_subset)){
          working_subset = working_subset[,.(case_control,quantile,ancestry_specific_pc1,ancestry_specific_pc2,
                                             ancestry_specific_pc3,
                                             phekb_age,site,phekb_family_history)]
          for(q in unique(working_subset$quantile)){
            if(q != '40-60%'){
              working_subset_q = working_subset[quantile == q | quantile == '40-60%']
              working_subset_q[,quantile:=factor(quantile,levels = c('40-60%',q))]
              working_subset_q = clean_column(working_subset_q)
              glm_model = run_glm(working_subset_q)
              dt = clean_results(glm_model = glm_model, working_subset_q = working_subset_q, prs_col = p,quantile = q,sex = s,ancestry = a, family_history = st)
              stratify_analysis_results = rbindlist(list(dt, stratify_analysis_results), fill=TRUE)
              
            }
          }
        }
      }
    }
  }
}
stratify_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/stratify_analysis_fx_qt.csv')

# subtype
stratify_subset_col = c(1,4,5,17:27,28,32,33,36,44:82)
prs_names = colnames(merged_table)[44:82]
stratify_subset = merged_table[, stratify_subset_col, with = FALSE]
stratify_subset = stratify_subset[phekb_case_control %in% c('case', 'control')]
stratify_subset = stratify_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
stratify_subset = stratify_subset[sex %in% c('male', 'female')]
stratify_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu', 'aa', 'la')) {
    for (s in c('male', 'female')) {
      for (st in c('positive', 'negative')) {
        working_subset = copy(stratify_subset)
        working_subset = working_subset[sex == s &
                                          emerge_ancestry == a
                                        &
                                          ((
                                            phekb_case_control == 'case' & hormone_receptor_string_value == st
                                          )
                                          |
                                            phekb_case_control == 'control'
                                          )]
        working_subset = prs_norm(working_subset,
                                  prs_col = p,
                                  case_control_col = 'phekb_case_control')
        working_subset = add_prs_quantile(working_subset,ancestry = a)
        
        if(!is.null(working_subset)){
          working_subset = working_subset[,.(case_control,quantile,ancestry_specific_pc1,ancestry_specific_pc2,
                                             ancestry_specific_pc3,
                                             phekb_age,site,phekb_family_history,hormone_receptor_string_value)]
          for(q in unique(working_subset$quantile)){
            if(q != '40-60%'){
              working_subset_q = working_subset[quantile == q | quantile == '40-60%']
              working_subset_q[,quantile:=factor(quantile,levels = c('40-60%',q))]
              working_subset_q = clean_column(working_subset_q)
              glm_model = run_glm(working_subset_q)
              dt = clean_results(glm_model = glm_model, working_subset_q = working_subset_q, prs_col = p,quantile = q,sex = s,ancestry = a, subtype = st)
              stratify_analysis_results = rbindlist(list(dt, stratify_analysis_results), fill=TRUE)
              
            }
          }
        }
      }
    }
  }
}
stratify_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/stratify_analysis_subtype_qt.csv')

# density
stratify_subset_col = c(1,4,5,17:27,28,32,33,50:52,65:82,38)
prs_names = colnames(merged_table)[c(50:52,65:82)]
stratify_subset = merged_table[, stratify_subset_col, with = FALSE]
stratify_subset = stratify_subset[phekb_case_control %in% c('case', 'control')]
stratify_subset = stratify_subset[emerge_ancestry %in% c('eu', 'aa', 'la')]
stratify_subset = stratify_subset[sex %in% c('male', 'female')]
stratify_analysis_results = NULL
for (p in prs_names) {
  for (a in c('eu', 'aa', 'la')) {
    for (s in c('male', 'female')) {
      for (st in c('scattered_fibro_dens',
                   'hetero_dense',
                   'extremely_dense',
                   'entirely_fat')) {
        working_subset = copy(stratify_subset)
        working_subset = working_subset[sex == s &
                                          emerge_ancestry == a
                                        &
                                          ((phekb_case_control == 'case' & breast_density == st)
                                           |
                                             phekb_case_control == 'control'
                                          )]
        working_subset = prs_norm(working_subset,
                                  prs_col = p,
                                  case_control_col = 'phekb_case_control')
        working_subset = add_prs_quantile(working_subset,ancestry = a)
        if(!is.null(working_subset)){
          for(q in unique(working_subset$quantile)){
            working_subset = working_subset[,.(case_control,quantile,ancestry_specific_pc1,ancestry_specific_pc2,
                                               ancestry_specific_pc3,
                                               phekb_age,site,phekb_family_history,breast_density)]
            if(q != '40-60%'){
              working_subset_q = working_subset[quantile == q | quantile == '40-60%']
              working_subset_q[,quantile:=factor(quantile,levels = c('40-60%',q))]
              working_subset_q = clean_column(working_subset_q)
              glm_model = run_glm(working_subset_q)
              dt = clean_results(glm_model = glm_model, working_subset_q = working_subset_q, prs_col = p,quantile = q,sex = s,ancestry = a, density = st)
              stratify_analysis_results = rbindlist(list(dt, stratify_analysis_results), fill=TRUE)
              
            }
          }
        }
      }
    }
  }
}
stratify_analysis_results %>% fwrite('~/2019-PRS-pipeline/workplace/stratify_analysis_density_qt.csv')


