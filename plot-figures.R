# Last updated: 2020-03-23
# Author: Cong Liu
# plot figures for PRS project.

rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(stringi)

point_estimator_plot = function(dt, stratify=NULL){
  
    dt = dt %>% 
      separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
      filter(sex == 'female' & filter == 'f3' & case_number > 10)
    pe_plot = dt %>% 
      ggplot(aes(y=logOR, x=ancestry, color=prs_model)) +
      geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = prs_model),position = position_dodge(.7)) + 
      scale_shape_manual(values = rep(c(4,8,15,16,17,18),3)) +
      coord_flip() + 
      theme_bw() + 
      theme(legend.position="right") +
      xlab("")
  if(!is.null(stratify)){
    warning("stratified by : ", stratify)
    pe_plot = pe_plot + facet_wrap(~get(stratify),scales = "free")
  }
  return(pe_plot)
}

quantile_plot = function(dt){
  dt = dt %>% 
    separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
    filter(sex == 'female' & filter == 'f3') %>%
    mutate(ancestry = factor(ancestry,levels = c('eu','aa','la')))
  qt_plot =   dt %>% ggplot(aes(x=quantile,y = logOR,group=prs_model, color=prs_model)) + 
    geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = prs_model),position = position_dodge(.5)) + 
    scale_shape_manual(values = rep(c(4,8,15,16,17,18),3)) +
    theme_bw() + 
    theme(legend.position="right") +
    xlab("") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 13)) + 
    facet_wrap(~ancestry,scales = 'free',ncol = 3)
  return(qt_plot)
}



## primary analysis
or = fread('~/2019-PRS-pipeline/workplace/primary_analysis_or.csv')
point_estimator_plot(or)
## sensitivity analysis icd
or = fread('~/2019-PRS-pipeline/workplace/sensitivity_analysis_icd_or.csv')
point_estimator_plot(or)
# stratify analysis by site
or = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_site_or.csv')
point_estimator_plot(or,stratify = 'site')
# stratify analysis by family history
or = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_fx_or.csv')
point_estimator_plot(or,stratify = 'family_history')
# stratify analysis by subtype
or = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_subtype_or.csv')
point_estimator_plot(or,stratify = 'subtype')
# stratify analysis by density
or = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_density_or.csv')
point_estimator_plot(or,stratify = 'density')

# quantile plot.
dt = fread('~/2019-PRS-pipeline/workplace/primary_analysis_qt.csv')
quantile_plot(dt)
dt = fread('~/2019-PRS-pipeline/workplace/sensitivity_analysis_icd_qt.csv')
quantile_plot(dt)
dt = fread('~/2019-PRS-pipeline/workplace/sensitivity_analysis_overall_qt.csv')
quantile_plot(dt)
dt = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_site_qt.csv')
dt = dt %>% filter(site == 'harv')
quantile_plot(dt)
dt = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_fx_qt.csv')
dt = dt %>% filter(family_history == 'w_fx')
quantile_plot(dt)
dt = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_subtype_qt.csv')
dt = dt %>% filter(subtype == 'positive')
quantile_plot(dt)
dt = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_density_qt.csv')
dt = dt %>% filter(density == 'extremely_dense')
quantile_plot(dt)

# auc plot




# 


# stratified figure.
stratify_analysis_or = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_site_or.csv')
stratify_analysis_or = stratify_analysis_or %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10)
point_estimator <- stratify_analysis_or %>% ggplot(aes(y=logOR, x=prs_model, color=ancestry)) +
  # add error bars, parameterized by other columns of 'estimates'
  # add point estimate, colored according to the problem column of 'estimate'  
  # geom_point(size = 3, shape = 1, position = position_dodge(.5)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = ancestry),position = position_dodge(.5)) + 
  coord_flip() + 
  theme_bw() + 
  xlab("") +
  facet_wrap(~site) + 
  theme(legend.position="right")
point_estimator

stratify_analysis_or = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_fx_or.csv')
stratify_analysis_or = stratify_analysis_or %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10)
point_estimator <- stratify_analysis_or %>% ggplot(aes(y=logOR, x=prs_model, color=ancestry)) +
  # add error bars, parameterized by other columns of 'estimates'
  # add point estimate, colored according to the problem column of 'estimate'  
  # geom_point(size = 3, shape = 1, position = position_dodge(.5)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = ancestry),position = position_dodge(.5)) + 
  coord_flip() + 
  theme_bw() + 
  xlab("") +
  facet_wrap(~family_history) + 
  theme(legend.position="right")
point_estimator

stratify_analysis_or = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_subtype_or.csv')
stratify_analysis_or = stratify_analysis_or %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10)
point_estimator <- stratify_analysis_or %>% ggplot(aes(y=logOR, x=prs_model, color=ancestry)) +
  # add error bars, parameterized by other columns of 'estimates'
  # add point estimate, colored according to the problem column of 'estimate'  
  # geom_point(size = 3, shape = 1, position = position_dodge(.5)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = ancestry),position = position_dodge(.5)) + 
  coord_flip() + 
  theme_bw() + 
  xlab("") +
  facet_wrap(~subtype) + 
  theme(legend.position="right")
point_estimator


stratify_analysis_or = fread('~/2019-PRS-pipeline/workplace/stratify_analysis_density_or.csv')
stratify_analysis_or = stratify_analysis_or %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10)
point_estimator <- stratify_analysis_or %>% ggplot(aes(y=logOR, x=prs_model, color=ancestry)) +
  # add error bars, parameterized by other columns of 'estimates'
  # add point estimate, colored according to the problem column of 'estimate'  
  # geom_point(size = 3, shape = 1, position = position_dodge(.5)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = ancestry),position = position_dodge(.5)) + 
  coord_flip() + 
  theme_bw() + 
  xlab("") +
  facet_wrap(~density) + 
  theme(legend.position="right")
point_estimator
