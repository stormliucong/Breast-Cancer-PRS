RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/";
TARGET_DIR=${RAW_DIR}/genomic-data;
PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink";
RESULTS_DIR='/home/cl3720/2019-PRS-pipeline/results/breast_cancer';
SCRIPT_DIR="/home/cl3720/2019-PRS-pipeline/script";
SOFTWARE_DIR="/home/cl3720/2019-PRS-pipeline/software";


cd ${RESULTS_DIR};
PHENOTYPE_FILE=${RESULTS_DIR}/'phenotype.AllBreastCancer'; # different for primary cancer and seconadry cancer
# TARGET_FILE=${TARGET_DIR}/qc/Merge;
TARGET_FILE=${TARGET_DIR}/bed/Merge;
for KEEP_FILE in ${RESULTS_DIR}/keep_*.txt;do
    KEEP=$(basename ${KEEP_FILE} .txt);
    # ${PLINK_DIR}/plink \
    # --bfile ${TARGET_DIR}/qc/Merge \
    # --pheno ${PHENOTYPE_FILE} \
    # --keep ${KEEP_FILE} \
    # --make-bed \
    # --out ${RESULTS_DIR}/${KEEP}.ldpred;
    for BASE_FILE in ${RESULTS_DIR}/PRS_*.txt;do
        BASE=$(basename ${BASE_FILE} .txt);
        ldpred coord \
        --rs SNP_NEW \
        --A1 effect_allele \
        --A2 reference_allele \
        --pos BP \
        --chr CHR \
        --pval P \
        --eff beta \
        --eff_type LOGOR \
        --ssf-format CUSTOM \
        --N 100000 \
        --ssf ${BASE_FILE} \
        --out ${RESULTS_DIR}/${BASE}.${KEEP}.coord \
        --gf ${RESULTS_DIR}/${KEEP}.ldpred;
    done;
done;