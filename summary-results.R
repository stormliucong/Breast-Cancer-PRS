# Last updated: 2020-04-03
# Author: Cong Liu
# summarize tables for PRS project.


rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(stringi)
library(tidyr)

# Table 2
dt = fread('./primary_analysis_auc.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(AUC_95CI = paste0(round(AUC,digits = 2),
                           ' (',round(AUC_lower95ci,digits = 2), 
                           '-', round(AUC_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female'& filter == 'f3' & case_number > 10) %>% as_tibble() %>%
  mutate(rsqr = round(rsqr * 10^3,2))
tb %>%
  dplyr::select(ancestry,prs_model,rsqr,AUC_95CI) %>%
  arrange(ancestry,prs_model)

# Supplementary Table S5
dt = fread('./primary_analysis_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female'& case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(ancestry,filter,prs_model,OR_95CI) %>%
  arrange(ancestry,filter,prs_model)

dt = fread('./primary_analysis_qt.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3') %>% 
  mutate(quantile = factor(as.character(quantile),levels = 
                             c(
                               '<1%',
                               '1-5%',
                               '5-10%',
                               '10-20%',
                               '<20%',
                               '20-40%',
                               '40-60%',
                               '60-80%',
                               '>=80%',
                               '80-90%',
                               '90-95%',
                               '95-99%',
                               '>=99%'
                             )
  )) %>%
  as_tibble()

tb %>% 
  select(prs_model,quantile,ancestry,case_number,control_number,OR_95CI) %>% 
  arrange(prs_model,ancestry,quantile)

# Supplementary Table S7
dt = fread('./sensitivity_analysis_overall_qt.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3') %>% 
  mutate(quantile = factor(as.character(quantile),levels = 
                             c(
                               '>=80%','>=90%','>=95%',">=99%"
                             )
  )) %>%
  as_tibble()

tb %>% 
  select(prs_model,quantile,ancestry,case_number,control_number,OR_95CI) %>% 
  arrange(prs_model,ancestry,quantile)

# Supplementary Table S8
dt = fread('./stratify_analysis_subtype_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3') %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,subtype,OR_95CI,case_number,control_number) %>%
  arrange(ancestry,subtype,prs_model)

# Supplementary Table S9
dt = fread('./stratify_analysis_fx_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3') %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,family_history,OR_95CI,case_number,control_number) %>%
  arrange(ancestry,family_history,prs_model)

# Supplementary Table S10
dt = fread('./stratify_analysis_density_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3') %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,density,OR_95CI,case_number,control_number) %>%
  arrange(ancestry,density,prs_model)
  
# Supplementary Table S11
dt = fread('./sensitivity_analysis_icd_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female'& filter == 'f3' & case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(ancestry,prs_model,OR_95CI) %>%
  spread(ancestry,OR_95CI) %>% arrange(prs_model) 

# Supplementary Table S12
dt = fread('./stratify_analysis_site_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,site,OR_95CI,case_number,control_number) %>%
  filter(
    (prs_model %in% c('bcac-s','bcac-l','ukbb')) | 
      (prs_model %in% c('whi-la','latinas') & ancestry == 'la') |
      (prs_model %in% c('whi-aa','root') & ancestry == 'aa')
  ) %>%
  arrange(ancestry,site,prs_model)

# 
# dt = fread('./interaction_analysis_results_fx.csv')
# tb = dt %>% 
#   separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
#   mutate(prs_score_OR_95CI = paste0(round(prs_score_OR,digits = 2),
#                                     ' (',round(prs_score_OR_lower95ci,digits = 2), 
#                                     '-', round(prs_score_OR_upper95ci,digits = 2),')')) %>%
#   mutate(prs_score_pvalue = p.adjust(prs_score_pvalue,method = 'fdr')) %>%
#   mutate(strata_OR_pvalue = p.adjust(strata_OR_pvalue,method = 'fdr')) %>%
#   mutate(interaction_pvalue = p.adjust(interaction_pvalue,method = 'fdr')) %>%
#   mutate(prs_score_pvalue = if_else(prs_score_pvalue < 1e-3,'<0.001',as.character(round(prs_score_pvalue,3)))) %>%
#   mutate(strata_OR_pvalue = if_else(strata_OR_pvalue < 1e-3,'<0.001',as.character(round(strata_OR_pvalue,3)))) %>%
#   mutate(interaction_pvalue = if_else(interaction_pvalue < 1e-3,'<0.001',as.character(round(interaction_pvalue,3)))) %>%
#   
#   mutate(strata_OR_95CI = paste0(round(strata_OR,digits = 2),
#                                  ' (',round(strata_OR_lower95ci,digits = 2), 
#                                  '-', round(strata_OR_upper95ci,digits = 2),')')) %>%
#   mutate(interaction_OR_95CI = paste0(round(interaction_OR,digits = 2),
#                                       ' (',round(interaction_OR_lower95ci,digits = 2), 
#                                       '-', round(interaction_OR_upper95ci,digits = 2),')')) %>%
#   filter(sex == 'female' & filter == 'f3' & case_number > 10) %>% as_tibble()
# 
# tb %>% select(prs_model,ancestry,prs_score_OR_95CI,prs_score_pvalue,
#               strata_OR_95CI,strata_OR_pvalue,
#               interaction_OR_95CI,interaction_pvalue) %>% 
#   arrange(ancestry,prs_model) %>%
#   fwrite('./table_s12.csv')
# 
# 
# 
# dt = fread('./interaction_analysis_results_site.csv')
# tb = dt %>% 
#   separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
#   filter(sex == 'female' & filter == 'f3' & case_number > 10 & ancestry == 'eu' & prs_model %in% c('bcac-l','bcac-s','ukbb')) %>% 
#   as_tibble() %>% 
#   select(prs_model,case_number,control_number,pvalue) %>% 
#   arrange(prs_model)
# 
# dt = fread('./interaction_analysis_results_fx.csv')
# tb = dt %>% 
#   separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
#   filter(sex == 'female' & filter == 'f3' & case_number > 10 & ancestry == 'eu' & prs_model %in% c('bcac-l','bcac-s','ukbb')) %>% 
#   as_tibble() %>% 
#   select(prs_model,case_number,control_number,pvalue) %>% 
#   arrange(prs_model)
#   
# dt = fread('./interaction_analysis_results_density.csv')
# tb = dt %>% 
#   separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
#   filter(sex == 'female' & filter == 'f3' & case_number > 10 & ancestry == 'eu' & prs_model %in% c('bcac-l','bcac-s','ukbb')) %>% 
#   as_tibble() %>% 
#   select(prs_model,case_number,control_number,pvalue) %>% 
#   arrange(prs_model)
# 
# 
