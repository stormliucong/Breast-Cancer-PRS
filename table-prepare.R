# Last updated: 2020-03-24
# Author: Cong Liu
# data processing script for PRS project
# generate big table for statistic analysis

rm(list=ls())
sort( sapply(ls(),function(x){object.size(get(x))}))  # memory usage
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

# Get all eMERGE ID and overall PCs.
overall_PC = fread('~/2019-PRS-pipeline/raw/genomic-data/chr1-22.plink_pca.e123_imputation_sample_manifest.csv')
colnames(overall_PC)[1] = 'eMERGE_ID'
colnames(overall_PC)[3] = 'ethnicity'
colnames(overall_PC)[4] = 'sex'


# Get eMERGE ancestry and ancestry specific PCs.
PC_AA = fread('~/2019-PRS-pipeline/raw/genomic-data/ancestry_specific_pcas/chr1-22.african.merged.maf.05.LD_1000_50_.7_pruned.no_tri.plink.pca-approx.eigenvec',sep = '\t')
PC_AA[,eMERGE_ancestry:='AA']
PC_EU = fread('~/2019-PRS-pipeline/raw/genomic-data/ancestry_specific_pcas/chr1-22.european.merged.maf.05.LD_1000_50_.7_pruned.no_tri.plink.pca-approx.eigenvec',sep = '\t')
PC_EU[,eMERGE_ancestry:='EU']
colnames(PC_EU)[1] = c('FID')
PC_LA = fread('~/2019-PRS-pipeline/raw/genomic-data/ancestry_specific_pcas/chr1-22.hispanic.merged.maf.05.LD_1000_50_.7_pruned.no_tri.plink.pca-approx.eigenvec',sep = '\t')
PC_LA[,eMERGE_ancestry:='LA']
PC_AS = fread('~/2019-PRS-pipeline/raw/genomic-data/ancestry_specific_pcas/chr1-22.asian.merged.maf.05.LD_1000_50_.7_pruned.no_tri.plink.pca-approx.eigenvec',sep = '\t')
PC_AS[,eMERGE_ancestry:='AS']
ancestry_specific_PC = rbind(PC_AA,PC_EU,PC_LA,PC_AS)
colnames(ancestry_specific_PC)[2] = 'eMERGE_ID'
colnames(ancestry_specific_PC)[3:12] = paste0('ancestry_specific_',colnames(ancestry_specific_PC)[3:12])
ancestry_specific_PC[,FID:=NULL]

# Get PheKB phenotypes.
# demo
demo = fread('~/2019-PRS-pipeline/workplace/eMERGE_phekb_demo.csv')
colnames(demo)[2:8] = c('phekb_age','phekb_sex','phekb_race','phekb_ethnicity','phekb_case_control','phekb_hr_status','phekb_family_history')
demo = demo[,!('phekb_hr_status'), with=FALSE]
demo = unique(demo)
# subtype
subtype = fread('~/2019-PRS-pipeline/workplace/eMERGE_phekb_subtype.csv')
subtype = subtype[Hormone_receptor == 'ER']
subtype = subtype[, .SD[which.max(Age_at_hormone_receptor_measurement)], by=eMERGE_ID]
subtype = subtype[,'Hormone_receptor':=NULL]
# density
density = fread('~/2019-PRS-pipeline/workplace/eMERGE_phekb_density.csv')
density = density[,.SD[which.max(Age_at_mammogram_report)],by=eMERGE_ID]

# Get ICD defined case-control and ICD-age
icd_table = fread('~/2019-PRS-pipeline/raw/phenotypic-data/latest/EMERGE_201907_ICD_GWAS_0.csv', header = TRUE)
# handling NA in icd_table.
icd_flag_table = icd_table[!ICD_FLAG %>% is.na(),.(1),c('ICD_CODE','ICD_FLAG')]
icd_table_na = icd_table[ICD_FLAG %>% is.na()]
icd_table_na_added = merge(icd_table_na,icd_flag_table,by=c('ICD_CODE'))
icd_table_na_added = icd_table_na_added[,c("ICD_FLAG"):=ICD_FLAG.y]
icd_table_na_added = icd_table_na_added[,c("ICD_FLAG.x","ICD_FLAG.y","V1"):=NULL]
icd_nona = icd_table[!ICD_FLAG %>% is.na()]
new_icd_table = rbindlist(list(icd_nona,icd_table_na_added),use.names=TRUE)
new_icd_table[,"CODE":=gsub("\\.","",ICD_CODE)]
new_icd_table = new_icd_table[,.(SUBJECT_ID,AGE_AT_OBSERVATION,CODE,ICD_FLAG)]
# all breast cancer events.
bc_icd_code = fread(file = '~/2019-PRS-pipeline/workplace/bc-icd-codes.csv',header = TRUE)
bc_icd_code[,"TYPE":=ifelse(TYPE=='ICD10CM',10,TYPE)]
bc_icd_code[,"TYPE":=ifelse(TYPE=='ICD9CM',9,TYPE)]
bc_icd_code[,"CODE":=gsub("\\.","",CODE)]
bc_icd_code[,"ICD_FLAG" := as.numeric(TYPE)]
bc_icd_code = bc_icd_code[,.(CODE,ICD_FLAG)]
# find icd case
emerge_table_bc_icd = merge(new_icd_table,bc_icd_code,by=c('CODE','ICD_FLAG') )
first_bc_icd = emerge_table_bc_icd[,.(min(AGE_AT_OBSERVATION)),(SUBJECT_ID)]
colnames(first_bc_icd) = c('eMERGE_ID','age_at_first_icd')
first_bc_icd[,'icd_case_control':='case']

# add icd family hisotry
family_history = icd_table[(ICD_CODE == "V16.3" & ICD_FLAG == "9") |(ICD_CODE == "Z80.3" & ICD_FLAG == "10")]
family_history[,icd_family_history:="w_fx"]
family_history = family_history[,.(SUBJECT_ID,icd_family_history)]
colnames(family_history)[1] = 'eMERGE_ID'
family_history = unique(family_history)

# Read PRS related tables.
prs_dir = '~/2019-PRS-pipeline/workplace/prs-scores/'
prs_scores_dt = NULL
for(prs_name in c('BCAC-S','BCAC-L','WHI-LA','WHI-AA','UKBB','ROOT','LATINAS',
                  'BCAC-S-ERP','BCAC-S-ERN','BCAC-S-H-ERP','BCAC-S-H-ERN','BCAC-L-ERP','BCAC-L-ERN')){
  for(f in c(1,2,3)){
    # read prs models.
    prs_scores = fread(paste0(prs_dir,prs_name,'_f',f,'.plink.profile'),sep = " ",header = TRUE)
    prs_scores = prs_scores[,.(IID,SCORE)]
    prs_scores[,prs_model:=paste0(prs_name,'-f',f)]
    prs_scores_dt = rbind(prs_scores,prs_scores_dt)
  }
}
prs_scores_dt_wide = dcast(prs_scores_dt,IID~prs_model,value.var='SCORE',fun.aggregate = sum)
colnames(prs_scores_dt_wide)[1] = c('eMERGE_ID')

# read brith decade
age_demo = fread(file = '~/2019-PRS-pipeline/raw/phenotypic-data/latest/EMERGE_201907_DEMO_GWAS_1.csv',header = TRUE)
age_demo = age_demo[,.(SUBJECT_ID,YEAR_BIRTH)]
colnames(age_demo)[1] = c('eMERGE_ID')

# merge all tables together.
# overall_PC
# ancestry_specific_PC
# demo
# subtype
# density
# first_bc_icd
# prs_scores_dt_wide
dt1 = merge(overall_PC,ancestry_specific_PC,by = 'eMERGE_ID',all.x = TRUE)
dt2 = merge(dt1,demo,by = 'eMERGE_ID',all.x = TRUE)
dt3 = merge(dt2,subtype,by = 'eMERGE_ID',all.x = TRUE)
dt4 = merge(dt3,density,by = 'eMERGE_ID',all.x = TRUE)
dt5 = merge(dt4,first_bc_icd,by = 'eMERGE_ID',all.x = TRUE)
dt6 = merge(dt5,family_history,by = 'eMERGE_ID',all.x = TRUE)
dt7 = merge(dt6,age_demo,by = 'eMERGE_ID',all.x = TRUE)
dt8 = merge(dt7,prs_scores_dt_wide,by = 'eMERGE_ID',all.x = TRUE)

colnames(dt8) = tolower(colnames(dt8))
colnames(dt8)

merged_table = dt8
for(col in colnames(merged_table)){
  if(class(merged_table[,get(col)]) == 'character'){
    # to lower case.
    merged_table[,eval(quote(col)):= as.factor(tolower(merged_table[,get(col)]))]
  }
}

# change family history to w_fx if found in either ICD or phenotyping algorithm.
merged_table[is.na(icd_family_history)| icd_family_history =='',icd_family_history:='wo_fx']
merged_table[is.na(icd_case_control)|icd_case_control=='',icd_case_control:='control']
merged_table[is.na(phekb_family_history)|phekb_family_history=='.'|phekb_family_history=='',phekb_family_history:=icd_family_history]


merged_table %>% fwrite('~/2019-PRS-pipeline/workplace/merged_table.txt')
merged_table = fread('~/2019-PRS-pipeline/workplace/merged_table.txt')
