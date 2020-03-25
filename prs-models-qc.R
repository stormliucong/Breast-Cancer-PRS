# Last updated: 2020-03-20
# Author: Cong Liu
# QC script for PRS project


rm(list=ls())
library(data.table)

# Generate control list for each ancestry.
# We will calculate HWE for variants in each ancestry using control samples only.
demo = fread('~/2019-PRS-pipeline/raw/genomic-data/chr1-22.plink_pca.e123_imputation_sample_manifest.csv')
case_control = fread('~/2019-PRS-pipeline/workplace/eMERGE_phekb_demo.csv')
case_control[,ancestry:=case_when(
  race=='White' & ethnicity=='Not Hispanic or Latino'~'EU',
  race=='Black or African American' & ethnicity=='Not Hispanic or Latino'~'AA',
  ethnicity=='Hispanic or Latino'~'LA'
)]
control = case_control[case_control_status == 'control']
control[ancestry=='EU'][,.(0,eMERGE_ID)] %>% fwrite('~/2019-PRS-pipeline/workplace/prs-models/EU_control.fam',sep = '\t',quote = F,row.names = F)
control[ancestry=='AA'][,.(0,eMERGE_ID)] %>% fwrite('~/2019-PRS-pipeline/workplace/prs-models/AA_control.fam',sep = '\t',quote = F,row.names = F)
control[ancestry=='LA'][,.(0,eMERGE_ID)] %>% fwrite('~/2019-PRS-pipeline/workplace/prs-models/LA_control.fam',sep = '\t',quote = F,row.names = F)


# Read var imp results.
var_impt = fread('~/2019-PRS-pipeline/raw/genomic-data/imputed_versus_genotyped_chip_counts_Rsq_metrics_by_chr_pos',sep = "\t",header = FALSE)
var_impt_simple = var_impt[,.(V1,V2,V4,V6,V9)]
colnames(var_impt_simple) = c('SNP_NEW','allele','genotyped','imputed','mean_Rsq')

# Executre prs-models-qc-helper.sh to generate prs-models.maf.frq.
# Read MAF results.
maf_freq = fread('~/2019-PRS-pipeline/workplace/genomeic-subset/prs-models.maf.frq')
maf_freq = maf_freq[V1!='CHR']
maf_freq_simple = maf_freq[,.(V2,V5)]
maf_freq_simple = data.table(maf_freq_simple)
maf_freq_simple$V5 = as.numeric(maf_freq_simple$V5)
colnames(maf_freq_simple) = c('SNP_NEW','MAF')

# Executre prs-models-qc-helper.sh to generate prs-models.maf.hwe
# Read hwe results.
hwe = fread('~/2019-PRS-pipeline/workplace/genomeic-subset/prs-models.maf.hwe')
hwe = hwe[V1!='CHR']
hwe_simple = hwe[,.(V2,V9)]
hwe_simple = data.table(hwe_simple)
hwe_simple$V9 = as.numeric(hwe_simple$V9)
colnames(hwe_simple) = c('SNP_NEW','HWE')
hwe_long = NULL
for(ancestry in c('EU','AA','LA')){
  hwe = fread(paste0('~/2019-PRS-pipeline/workplace/prs-models/bed/prs-models/prs-models.',ancestry,'_control.hwe'))
  hwe = hwe[V1!='CHR']
  hwe_simple = hwe[,.(V2,V9)]
  hwe_simple = data.table(hwe_simple)
  hwe_simple$V9 = as.numeric(hwe_simple$V9)
  colnames(hwe_simple) = c('SNP_NEW','HWE')
  hwe_simple[,A:=ancestry]
  hwe_long = rbind(hwe_simple,hwe_long)
}
hwe_wide = dcast(hwe_long, SNP_NEW~A,value.var='HWE',fun.aggregate = sum)

# Executre prs-models-qc-helper.sh to generate Merge.missnp
# Read 3+ allele SNPs.
missnp = fread('~/2019-PRS-pipeline/workplace/genomeic-subset/Merge.missnp',header = FALSE)
colnames(missnp) = 'SNP_NEW'

# Executre prs-models-qc-helper.sh to generate Merge.exclude.missnp.*
# Read mismatch SNPs.
bim <- fread("~/2019-PRS-pipeline/workplace/genomeic-subset/Merge.exclude.missnp.bim")
setnames(bim, colnames(bim), c("CHR", "SNP_NEW", "CM", "BP", "B.A1", "B.A2"))

################################
# filter functions start here #
################################
# print var passed ambiguous filter 
ambiguous_filter = function(dt){
  dt = dt[!(
    (effect_allele=='A' & reference_allele == 'T') | 
      (effect_allele=='T' & reference_allele == 'A') |
      (effect_allele=='G' & reference_allele == 'C') |
      (effect_allele=='C' & reference_allele == 'G') )]
  cat('SNP left after applying ambiguous filter: ;;',dt$SNP_NEW %>% unique() %>% length(),'\n')
  return(dt)
}
# print var passed imputation filter.
imputation_filter = function(dt, rsq=0.8){
  info = merge(dt,var_impt_simple, by = 'SNP_NEW',all.x=TRUE)
  pass_snps = info[mean_Rsq > rsq,SNP_NEW]
  dt = dt[SNP_NEW %in%pass_snps]
  dt$SNP_NEW %>% unique() %>% length()
  cat('SNP left after applying imputation filter: with Rsq ',rsq, ':;;',dt$SNP_NEW %>% unique() %>% length(),'\n')
  return(dt)
}
# print var passed MAF filter 
# MAF won't affect the SNP filtering.
maf_filter = function(dt,maf=0.005){
  info = merge(dt,maf_freq_simple, by = 'SNP_NEW',all.x=TRUE)
  pass_snps = info[MAF >= maf,SNP_NEW]
  dt = dt[SNP_NEW %in%pass_snps]
  dt$SNP_NEW %>% unique() %>% length()
  cat('SNP left after applying MAF filter: with MAF ',maf, ':;;',dt$SNP_NEW %>% unique() %>% length(),'\n')
  return(dt)
}
# HWE filter
hwe_filter = function(dt,hwe=1e-6,ancestry='ALL'){
  info = merge(dt,hwe_wide, by = 'SNP_NEW',all.x=TRUE)
  if(ancestry=='ALL'){
    pass_snps = info[AA > hwe & EU > hwe & LA > hwe,SNP_NEW]
    dt = dt[SNP_NEW %in%pass_snps]
    dt$SNP_NEW %>% unique() %>% length()
    cat('SNP left after applying HWE filter (ALL): with HWE ',hwe, ':;;',dt$SNP_NEW %>% unique() %>% length(),'\n')
  }
  if(ancestry=='ANY'){
    pass_snps = info[AA > hwe | EU > hwe | LA > hwe,SNP_NEW]
    dt = dt[SNP_NEW %in%pass_snps]
    dt$SNP_NEW %>% unique() %>% length()
    cat('SNP left after applying HWE filter (ANY): with HWE ',hwe, ':;;',dt$SNP_NEW %>% unique() %>% length(),'\n')
  }
  return(dt)
}
# multi allele filter
multi_allele_filter = function(dt){
  dt = dt[!SNP_NEW %in% missnp$SNP_NEW]
  cat('SNP left after applying multi-allele filter: ;;',dt$SNP_NEW %>% unique() %>% length(),'\n')
  return(dt)
}
# mismatch SNP filter and update target bim by flip strand.
mismatch_filter = function(dt){
  info = merge(dt, bim, by = c("SNP_NEW"),all.x=TRUE)
  complement = function(x) {
    switch (
      x,
      "A" = "T",
      "C" = "G",
      "T" = "A",
      "G" = "C",
      return(NA)
    )
  }
  # Get SNPs that have the same alleles across base and target but with complementary
  com_snps = info[sapply(effect_allele, complement) == B.A1 & sapply(reference_allele, complement) == B.A2, SNP_NEW]
  com_recode = info[sapply(effect_allele, complement) == B.A2 & sapply(reference_allele, complement) == B.A1, SNP_NEW]
  # Now update the dt file
  dt[SNP_NEW %in% c(com_snps,com_recode), c("effect_allele", "reference_allele") :=
       list(sapply(effect_allele, complement),
            sapply(reference_allele, complement))]
  # And update the info structure
  info[SNP_NEW %in% c(com_snps,com_recode), c("effect_allele", "reference_allele") :=
         list(sapply(effect_allele, complement),
              sapply(reference_allele, complement))]
  matched = info[(effect_allele == B.A1 & reference_allele == B.A2) |
                   (reference_allele == B.A1 & effect_allele == B.A2)]
  dt = dt[SNP_NEW%in%matched$SNP_NEW]
  cat('SNP left after applying allele mismatch filter: ;;',dt$SNP_NEW %>% unique() %>% length(),'\n')
  return(dt)
}
# chain filter to evalutate single filter effect.
filter_each = function(dt){
  dt2 = ambiguous_filter(dt)
  dt4 = mismatch_filter(dt)
  dt3 = multi_allele_filter(dt)
  dt1 = imputation_filter(dt,rsq = 0) # 0, 0.3, 0.8
  dt1 = imputation_filter(dt,rsq = 0.3) # 0, 0.3, 0.8
  dt1 = imputation_filter(dt,rsq = 0.8) # 0, 0.3, 0.8
  dt6 = hwe_filter(dt,ancestry = 'ANY')
  dt6 = hwe_filter(dt,ancestry = 'ALL')
  dt5 = maf_filter(dt,maf = 0)
  dt6 = maf_filter(dt,maf = 0.005)
  dt6 = maf_filter(dt,maf = 0.01)
  dt6 = maf_filter(dt,maf = 0.05)
}
# stack filter together for different filter pipeline.
filter_stacks = function(dt,filter=1){
  if(filter == 1){
    dt = ambiguous_filter(dt)
    dt = mismatch_filter(dt)
    dt = multi_allele_filter(dt)
    dt = imputation_filter(dt,rsq = 0) # 0, 0.3, 0.8
  }
  if(filter == 2){
    dt = ambiguous_filter(dt)
    dt = mismatch_filter(dt)
    dt = multi_allele_filter(dt)
    dt = imputation_filter(dt,rsq = 0.3) # 0, 0.3, 0.8
    dt = maf_filter(dt,maf = 0.005)
    dt = hwe_filter(dt,ancestry = 'ANY')
  }
  if(filter == 3){
    dt = ambiguous_filter(dt)
    dt = mismatch_filter(dt)
    dt = multi_allele_filter(dt)
    dt = imputation_filter(dt,rsq = 0.8) # 0, 0.3, 0.8
    dt = maf_filter(dt,maf = 0.01)
    dt = hwe_filter(dt,ancestry = 'ALL')
  }
  return(dt)
}
# generate plink input
plink_format_dt = function(dt,subtype = 'ALL'){
  if(subtype == 'ALL'){
    dt = dt[,.(SNP_NEW,CHR,BP,effect_allele,reference_allele,ALL)]
  }
  if(subtype == 'ERP'){
    dt = dt[,.(SNP_NEW,CHR,BP,effect_allele,reference_allele,`ER-positive`)]
  }
  if(subtype == 'ERN'){
    dt = dt[,.(SNP_NEW,CHR,BP,effect_allele,reference_allele,`ER-negative`)]
  }
  if(subtype == 'H-ERP'){
    dt = dt[,.(SNP_NEW,CHR,BP,effect_allele,reference_allele,`hybrid ER-positive`)]
  }
  if(subtype == 'H-ERN'){
    dt = dt[,.(SNP_NEW,CHR,BP,effect_allele,reference_allele,`hybrid ER-negative`)]
  }
  return(dt)
}

# Generate PRS plink input files after differet QC filters.
dir = '~/2019-PRS-pipeline/workplace/prs-models/'
for(prs_name in c('BCAC-S','BCAC-L','WHI-LA','WHI-AA','UKBB','ROOT','LATINAS')){
  cat('Perform QC filter for prs model: ;;', prs_name)
  # effect of single model.
  # filter_each(bcac_l) # report each filter effect.
  for(f in c(1,2,3)){
    # read prs models.
    prs = fread(paste0(dir,prs_name,'.txt'),sep = "\t",header = TRUE)
    if(prs_name == 'BCAC-S'){
      for(s in c('ERP','ERN','H-ERP','H-ERN')){
        filter_stacks(prs,filter = f) %>% plink_format_dt(subtype = s) %>%
          fwrite(file = paste0(dir,prs_name,'-',s,'_f',f,'.txt'),row.names = FALSE,sep = '\t')
      }
    }
    if(prs_name == 'BCAC-L'){
      for(s in c('ERP','ERN')){
        filter_stacks(prs,filter = f) %>% plink_format_dt(subtype = s) %>%
          fwrite(file = paste0(dir,prs_name,'-',s,'_f',f,'.txt'),row.names = FALSE,sep = '\t')
      }
    }
    filter_stacks(prs,filter = f) %>% plink_format_dt() %>%
      fwrite(file = paste0(dir,prs_name,'_f',f,'.txt'),row.names = FALSE,sep = '\t')
  }
}







