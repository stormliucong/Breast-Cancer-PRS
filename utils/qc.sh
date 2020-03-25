RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/"
TARGET_DIR=${RAW_DIR}/genomic-data
PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink"
QC_DIR=${TARGET_DIR}/qc
cd ${QC_DIR}

# Standard GWAS QC
for chr in $(seq 1 22); do
    ${PLINK_DIR}/plink \
    --bfile ${TARGET_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-bed \
    --out ${QC_DIR}/chr${chr}.dose.emerge_ids.consented.merged.qc ;
done
