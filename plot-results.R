# Last updated: 2020-05-07
# Author: Cong Liu
# plot figures for PRS project.

rm(list=ls())
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(stringi)

quantile_plot = function(dt){
  dt = dt %>% 
    separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
    filter(sex == 'female' & filter == 'f3' & ancestry == 'eu') %>%
    mutate(ancestry = factor(ancestry,levels = c('European','African American','Latina'))) %>%
    mutate(prs_model = toupper(prs_model))
    
  qt_plot =   dt %>% ggplot(aes(x=quantile,y = logOR,group=prs_model, color=prs_model)) + 
    geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci, shape = prs_model),position = position_dodge(.5)) + 
    scale_shape_manual(values = rep(c(8,15,16,17),3)) +
    theme_bw() + 
    theme(legend.position="right",text = element_text(size=20)) +
    xlab("") + 
    ylab("Odds Ratio (log scale)") + 
    labs(color = 'PRS Model',shape = 'PRS Model') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 13)) 
    # + facet_wrap(~ancestry,scales = 'free',ncol = 3)
  return(qt_plot)
}



# Figure 1
or = fread('./primary_analysis_or.csv')
dt = or %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>%
  mutate(`PRS Model` = prs_model)
dt = dt %>% filter(
  (prs_model %in% c('bcac-s','bcac-l','ukbb')) | 
    (prs_model %in% c('whi-la','latinas') & ancestry == 'la') |
    (prs_model %in% c('whi-aa','root') & ancestry == 'aa')
) %>% mutate(ancestry=case_when(
  ancestry == 'aa'~'(B) African American',
  ancestry == 'la'~'(C) LatinX',
  ancestry == 'eu'~'(A) European'
)) %>% 
  mutate(`PRS Model` = toupper(`PRS Model`))
pe_plot = dt %>% 
  ggplot(aes(y=logOR, x=`PRS Model`)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci),position = position_dodge(.7)) + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position="right",text = element_text(size=20)) +
  ylab("Odds Ratio (log scale)") +
  xlab('') + 
  labs(color = 'PRS Model',shape = 'PRS Model') +
  scale_x_discrete(labels=c("eu" = "European", 
                            "la" = "Latina",
                            "aa" = "African American")) + 
  facet_wrap(~ancestry,nrow = 3,scales = 'free_y')
pe_plot

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
                                                       'bcac-l-ern-f3')) %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>%
  mutate(`PRS Model` = prs_model)
pe_plot = or %>% 
  filter(
    (prs_model %in% c('ukbb','bcac-s','bcac-l','bcac-s-h-erp','bcac-s-erp','bcac-l-erp') & subtype == 'positive') | 
    (prs_model %in% c('ukbb','bcac-s','bcac-l','bcac-s-h-ern','bcac-s-ern','bcac-l-ern') & subtype == 'negative')
  ) %>% 
  mutate(subtype = case_when(
    subtype == 'negative'~'(B) Estrogen Receptor (ER)-Negative',
    subtype == 'positive'~'(A) Estrogen Receptor (ER)-Positive'
  )) %>%
  mutate(prs_model = toupper(prs_model)) %>%
  ggplot(aes(y=logOR, x=prs_model)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci),position = position_dodge(.7)) + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position="right",text = element_text(size=20)) +
  xlab("") + 
  ylab("Odds Ratio (log scale)") +
  facet_wrap(~subtype,scales = 'free_y',ncol = 1)
pe_plot

# Figure 5
dt = fread('./stratify_analysis_fx_or.csv')
dt = dt %>% mutate(family_history=case_when(family_history=='w_fx'~'with a family history',
                                            TRUE~'without a family history'))
dt = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>%
  mutate(`PRS Model` = prs_model)
dt = dt %>% filter(
  (prs_model %in% c('bcac-s','bcac-l','ukbb')) | 
    (prs_model %in% c('whi-la','latinas') & ancestry == 'la') |
    (prs_model %in% c('whi-aa','root') & ancestry == 'aa')
) %>% mutate(ancestry=case_when(
  ancestry == 'aa'~'(B) African American',
  ancestry == 'la'~'(C) LatinX',
  ancestry == 'eu'~'(A) European'
)) %>% 
  mutate(`PRS Model` = toupper(`PRS Model`))
pe_plot = dt %>% 
  ggplot(aes(y=logOR, x=`PRS Model`,color=family_history)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci,shape = family_history),position = position_dodge(.7)) + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position="top",text = element_text(size=20)) +
  ylab("Odds Ratio (log scale)") +
  xlab('') + 
  labs(color = '',shape = '') +
  facet_wrap(~ancestry,ncol = 1,scales = 'free_y')
pe_plot

# Figure 6
dt = fread('./stratify_analysis_density_or.csv')
dt = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>%
  mutate(`PRS Model` = prs_model)
dt = dt %>% filter(
  (prs_model %in% c('bcac-s','bcac-l','ukbb')) | 
    (prs_model %in% c('whi-la','latinas') & ancestry == 'la') |
    (prs_model %in% c('whi-aa','root') & ancestry == 'aa')
) %>% mutate(ancestry=case_when(
  ancestry == 'aa'~'(B) African American',
  ancestry == 'la'~'(C) LatinX',
  ancestry == 'eu'~'(A) European'
)) %>% 
  mutate(`PRS Model` = toupper(`PRS Model`))
pe_plot = dt %>% 
  ggplot(aes(y=logOR, x=`PRS Model`,color=density)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci,shape = density),position = position_dodge(.7)) + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position="top",text = element_text(size=20)) +
  ylab("Odds Ratio (log scale)") +
  xlab('') + 
  labs(color = '',shape = '') +
  facet_wrap(~ancestry,ncol = 1,scales = 'free_y')
pe_plot
# Figure 7
dt = fread('./stratify_analysis_site_or.csv')
dt = dt %>% 
  separate(prs_model, into = c("prs_model", "filter"), sep = "\\-(?!.*-)", remove = FALSE) %>%
  filter(sex == 'female' & filter == 'f3' & case_number > 10) %>%
  mutate(`PRS Model` = prs_model)
dt = dt %>% filter(
  (prs_model %in% c('bcac-s','bcac-l','ukbb')) | 
    (prs_model %in% c('whi-la','latinas') & ancestry == 'la') |
    (prs_model %in% c('whi-aa','root') & ancestry == 'aa')
) %>% mutate(ancestry=case_when(
  ancestry == 'aa'~'(B) African American',
  ancestry == 'la'~'(C) LatinX',
  ancestry == 'eu'~'(A) European'
)) %>% 
  mutate(`PRS Model` = toupper(`PRS Model`))
pe_plot = dt %>% 
  ggplot(aes(y=logOR, x=`PRS Model`,color=site)) +
  geom_pointrange(aes(ymin = logOR_lower95ci, ymax = logOR_upper95ci),size=.3,position = position_dodge(.7)) + 
  coord_flip() + 
  theme_bw() + 
  scale_shape_manual(values = rep(c(8,15,16),3)) +
  theme(legend.position="top",text = element_text(size=20)) +
  ylab("Odds Ratio (log scale)") +
  xlab('') + 
  labs(color = '',shape = '') +
  facet_wrap(~ancestry,ncol = 1,scales = 'free_y')
pe_plot


