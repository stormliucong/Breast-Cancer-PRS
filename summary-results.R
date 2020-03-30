# Last updated: 2020-03-30
# Author: Cong Liu
# summarize tables for PRS project.


rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(stringi)
library(tidyr)


# Table S6
dt = fread('./primary_analysis_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'male' & ancestry == 'eu' & filter == 'f3' & case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(prs_model,OR_95CI,case_number,control_number) %>% 
  arrange(prs_model) %>%
  fwrite('./table_s6.csv')

# Table S7
dt = fread('./cox_analysis.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(HR_95CI = paste0(round(`exp(coef)`,digits = 2),
                          ' (',round(`lower .95`,digits = 2 ), 
                          '-', round(`upper .95`,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,HR_95CI,case_number,control_number) %>% 
  arrange(ancestry,prs_model) %>%
  fwrite('./table_s7.csv')

# Table S8
dt = fread('./primary_analysis_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female'& case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(ancestry,filter,prs_model,OR_95CI) %>%
  arrange(ancestry,filter,prs_model) %>% fwrite('./table_s8.csv')

# Table S9
dt = fread('./sensitivity_analysis_icd_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female'& filter == 'f3' & case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(ancestry,prs_model,OR_95CI) %>%
  spread(ancestry,OR_95CI) %>% arrange(prs_model) %>% fwrite('./table_s9.csv')

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
  arrange(ancestry,prs_model) %>% fwrite('./table_2.csv')

# Supplementary Table S10
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
  arrange(prs_model,ancestry,quantile) %>%
  fwrite('./table_s10.csv')

# Table S11
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
  arrange(prs_model,ancestry,quantile) %>%
  fwrite('./table_s11.csv')

# Table S12
dt = fread('./stratify_analysis_subtype_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3') %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,subtype,OR_95CI,case_number,control_number) %>%
  arrange(ancestry,subtype,prs_model) %>%
  fwrite('./table_s12.csv')

# Table S13
dt = fread('./stratify_analysis_fx_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3') %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,family_history,OR_95CI,case_number,control_number) %>%
  arrange(ancestry,family_history,prs_model) %>%
  fwrite('./table_s13.csv')

# Table S14
dt = fread('./stratify_analysis_density_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3') %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,density,OR_95CI,case_number,control_number) %>%
  arrange(ancestry,density,prs_model) %>%
  fwrite('./table_s14.csv')

# Table S15
dt = fread('./stratify_analysis_site_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                          ' (',round(OR_lower95ci,digits = 2), 
                          '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(prs_model,ancestry,site,OR_95CI,case_number,control_number) %>%
  arrange(ancestry,site,prs_model) %>%
fwrite('./table_s15.csv')



