RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/"
TARGET_DIR=${RAW_DIR}/genomic-data
PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink"

# chr=22
# ${PLINK_DIR}/plink \
#   --bfile ${TARGET_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged \
#   --maf 0.10 --indep 50 5 1.5 \
#   --out ${TARGET_DIR}/pruned/chr${chr}.dose.emerge_ids.consented.merged ;

# ${PLINK_DIR}/plink \
#   --bfile ${TARGET_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged \
#   --extract ${TARGET_DIR}/pruned/chr${chr}.dose.emerge_ids.consented.merged.prune.in \
#   --make-bed \
#   --out ${TARGET_DIR}/pruned/chr${chr}.dose.emerge_ids.consented.merged ;


for chr in $(seq 1 22); do
    ${PLINK_DIR}/plink \
      --bfile ${TARGET_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged \
      --maf 0.10 --indep 50 5 1.5 \
      --out ${TARGET_DIR}/pruned/chr${chr}.dose.emerge_ids.consented.merged ;

    ${PLINK_DIR}/plink \
      --bfile ${TARGET_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged \
      --extract ${TARGET_DIR}/pruned/chr${chr}.dose.emerge_ids.consented.merged.prune.in \
      --make-bed \
      --out ${TARGET_DIR}/pruned/chr${chr}.dose.emerge_ids.consented.merged ;
done