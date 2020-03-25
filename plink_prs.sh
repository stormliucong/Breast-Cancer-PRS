# Use Plink to calculate PRS score for each prs models.

PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink";
PRS_MODELS_DIR='/home/cl3720/2019-PRS-pipeline/workplace/prs-models';
GSUBSET_DIR='/home/cl3720/2019-PRS-pipeline/workplace/genomeic-subset'
SCORE_DIR='/home/cl3720/2019-PRS-pipeline/workplace/prs-scores'
cd ${PRS_MODELS_DIR};

for BASE_FILE in ${PRS_MODELS_DIR}/*_f*.txt;do
    BASE=$(basename ${BASE_FILE} .txt);
    # 1,4,9 - SNP_ID, effective_allele, effective size.
    ${PLINK_DIR}/plink \
        --bfile ${GSUBSET_DIR}/Merge.exclude.missnp \
        --score ${BASE_FILE} 1 4 6 header \
        --out ${SCORE_DIR}/${BASE}.plink;
done