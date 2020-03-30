# Last updated: 2020-03-24
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
    filter(sex == 'female' & filter == 'f3' & case_number > 10) %>%
    mutate(`PRS Model` = prs_model)
  pe_plot = dt %>% 
    ggplot(aes(y=logOR, x=ancestry, color=`PRS Model`)) +
    geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = `PRS Model`),position = position_dodge(.7)) + 
    scale_shape_manual(values = rep(c(8,15,16,17),4)) +
    coord_flip() + 
    theme_bw() + 
    theme(legend.position="right",text = element_text(size=20)) +
    xlab("") + 
    labs(color = 'PRS Model',shape = 'PRS Model') +
    scale_x_discrete(labels=c("eu" = "European", 
                              "la" = "Latina",
                              "aa" = "African American"))

  if(!is.null(stratify)){
    warning("stratified by : ", stratify)
    pe_plot = pe_plot + facet_wrap(~get(stratify),nrow = 2,scales = 'free_x')
  }
  
  return(pe_plot)
}

quantile_plot = function(dt){
  dt = dt %>% 
    separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
    filter(sex == 'female' & filter == 'f3' & ancestry == 'eu') %>%
    mutate(ancestry = factor(ancestry,levels = c('European','African American','Latina')))
  qt_plot =   dt %>% ggplot(aes(x=quantile,y = logOR,group=prs_model, color=prs_model)) + 
    geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = prs_model),position = position_dodge(.5)) + 
    scale_shape_manual(values = rep(c(8,15,16,17),3)) +
    theme_bw() + 
    theme(legend.position="right",text = element_text(size=20)) +
    xlab("") + 
    labs(color = 'PRS Model',shape = 'PRS Model') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 13)) 
    # + facet_wrap(~ancestry,scales = 'free',ncol = 3)
  return(qt_plot)
}



# Figure 1
or = fread('./primary_analysis_or.csv')
point_estimator_plot(or)
# Figure 2
dt = fread('./primary_analysis_qt.csv')
dt = dt %>% filter(prs_model %in% c('ukbb-f3','bcac-l-f3','bcac-s-f3')) %>%
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
))
quantile_plot(dt)

# Figure 3
dt = fread('./sensitivity_analysis_overall_qt.csv')
dt = dt %>% filter(prs_model %in% c('ukbb-f3','bcac-l-f3','bcac-s-f3')) %>%
  mutate(quantile = factor(as.character(quantile),levels =
                             c(
                               '>=80%',
                               '>=90%',
                               '>=95%',
                               '>=99%'
                             )
  ))
quantile_plot(dt)

# Figure 4
or = fread('./stratify_analysis_subtype_or.csv')
or = or %>% 
  filter(ancestry == 'eu' & prs_model %in% c('ukbb-f3',
                                                       'bcac-s-f3',
                                                       'bcac-l-f3',
                                                       'bcac-s-erp-f3',
                                                       'bcac-s-ern-f3',
                                                       'bcac-s-h-erp-f3',
                                                       'bcac-s-h-ern-f3',
                                                       'bcac-l-erp-f3',
                                                       'bcac-l-erp-f3')) %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>%
  mutate(`PRS Model` = prs_model)
pe_plot = or %>% 
  ggplot(aes(y=logOR, x=subtype, color=`PRS Model`)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = `PRS Model`),position = position_dodge(.7)) + 
  scale_shape_manual(values = rep(c(8,15,16,17),4)) +
  coord_flip() + 
  theme_bw() + 
  theme(legend.position="right",text = element_text(size=20)) +
  xlab("") + 
  labs(color = 'PRS Model',shape = 'PRS Model')
pe_plot

# Figure 5
dt = fread('./stratify_analysis_fx_or.csv')
dt = dt %>% mutate(family_history=case_when(family_history=='w_fx'~'with family history',
                                            TRUE~'without family history'))
point_estimator_plot(dt,stratify = 'family_history')

# Figure 6
dt = fread('./stratify_analysis_density_or.csv')
dt = dt %>% mutate(density=case_when(density=='entirely_fat'~'entirely fatty',
                                            density=='extremely_dense'~'extremely dense',
                                            density=='hetero_dense'~'hetero dense',
                                            TRUE~'scattered fibro dense'))
point_estimator_plot(dt,stratify = 'density')

# Figure 7
or = fread('./stratify_analysis_site_or.csv')
point_estimator_plot(or,stratify = 'site')


