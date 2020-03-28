# Last updated: 2020-03-26
# Author: Cong Liu
# summarize tables for PRS project.


rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(stringi)


dt = fread('~/2019-PRS-pipeline/workplace/primary_analysis_or.csv')
tb = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  mutate(OR_95CI = paste0(round(OR,digits = 2),
                             ' (',round(OR_lower95ci,digits = 2), 
                             '-', round(OR_upper95ci,digits = 2),')')) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>% as_tibble()
tb %>%
  dplyr::select(ancestry,prs_model,OR_95CI) %>%
  spread(ancestry,OR_95CI) %>% fwrite(~)


