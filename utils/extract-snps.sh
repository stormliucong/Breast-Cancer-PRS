RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/";
TARGET_DIR=${RAW_DIR}/genomic-data;
PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink";
# RESULTS_DIR='/home/cl3720/2019-PRS-pipeline/results/breast_cancer';
RESULTS_DIR='/home/cl3720/2019-PRS-pipeline/sandbox'
SCRIPT_DIR="/home/cl3720/2019-PRS-pipeline/script";
SOFTWARE_DIR="/home/cl3720/2019-PRS-pipeline/software";

for chr in $(seq 1 22); do
       ${PLINK_DIR}/plink \
           --bfile ${TARGET_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged \
           --extract ${RESULTS_DIR}/PRS_313_full.valid.snprange --range \
           --make-bed \
           --out ${RESULTS_DIR}/chr${chr}.dose.emerge_ids.consented.merged.prs313 ;
    done