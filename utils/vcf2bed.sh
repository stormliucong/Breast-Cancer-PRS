RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/"
TARGET_DIR=${RAW_DIR}/genomic-data/merged_vcf_gz/
PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink"

for chr in $(seq 1 22); do
        ${PLINK_DIR}/plink \
        --vcf ${TARGET_DIR}/chr${chr}.dose.emerge_ids.consented.merged.vcf.gz \
        --out ${TARGET_DIR}/../bed/chr${chr}.dose.emerge_ids.consented.merged.plink \
        --make-bed &
done
wait

# for chr in $(seq 2 22); do
#         ${PLINK_DIR}/plink \
#         --vcf ${TARGET_DIR}/chr${chr}.dose.emerge_ids.consented.merged.vcf.gz \
#         --out ${TARGET_DIR}/chr${chr}.dose.emerge_ids.consented.merged.plink \
#         --make-just-bim &
# done
# wait

# for chr in $(seq 9 9); do
#         ${PLINK_DIR}/plink \
#         --vcf ${TARGET_DIR}/chr${chr}.dose.emerge_ids.consented.merged.vcf.gz \
#         --out ${TARGET_DIR}/../bed/chr${chr}.dose.emerge_ids.consented.merged.plink \
#         --make-bed
# done
# wait
