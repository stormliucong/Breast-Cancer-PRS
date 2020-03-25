RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/";
TARGET_DIR=${RAW_DIR}/genomic-data;
PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink";
RESULTS_DIR='/home/cl3720/2019-PRS-pipeline/results/breast_cancer';
# RESULTS_DIR='/home/cl3720/2019-PRS-pipeline/sandbox'
SCRIPT_DIR="/home/cl3720/2019-PRS-pipeline/script";
SOFTWARE_DIR="/home/cl3720/2019-PRS-pipeline/software";
BASE_DIR=${RAW_DIR}/base/breast_cancer;

cd ${RESULTS_DIR};
PHENOTYPE_FILE=${RESULTS_DIR}/'phenotype.AllBreastCancer'; # different for primary cancer and seconadry cancer
# TARGET_FILE=${TARGET_DIR}/qc/Merge;
# TARGET_FILE=${TARGET_DIR}/bed/Merge;

# remove ambigous SNPs from base file




for BASE_FILE in ${RESULTS_DIR}/PRS_*.txt;do
    BASE=$(basename ${BASE_FILE} .txt);
    less ${BASE_FILE} | awk '{print $1}' > ${RESULTS_DIR}/${BASE}.valid.snp;
    less ${BASE_FILE} | awk 'NR>1 {print $2,$3,$3,$1}' > ${RESULTS_DIR}/${BASE}.valid.snprange;
    mkdir -p ${RESULTS_DIR}/bed/${BASE};
    cd ${RESULTS_DIR}/bed/${BASE};
    for chr in $(seq 1 22); do
       ${PLINK_DIR}/plink \
           --bfile ${TARGET_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged \
           --extract ${RESULTS_DIR}/${BASE}.valid.snprange --range \
           --make-bed \
           --out ${RESULTS_DIR}/bed/${BASE}/chr${chr}.dose.emerge_ids.consented.merged ;
    done

    find . -name "*.bim" > ForMerge.list ;
    sed -i 's/.bim//g' ForMerge.list ;
    ${PLINK_DIR}/plink --merge-list ForMerge.list --out ${RESULTS_DIR}/bed/${BASE}/Merge ;
    
    # cd ${RESULTS_DIR};
    # # 1,4,9 - SNP_ID, effective_allele, effective size.
    # ${PLINK_DIR}/plink \
    #     --bfile ${RESULTS_DIR}/bed/Merge \
    #     --score ${BASE_FILE} 1 4 9 header \
    #     --out ${RESULTS_DIR}/${BASE}.erp.plink.out;

    # # 1,4,11 - SNP_ID, effective_allele, effective size.
    # ${PLINK_DIR}/plink \
    #     --bfile ${RESULTS_DIR}/bed/Merge \
    #     --score ${BASE_FILE} 1 4 11 header \
    #     --out ${RESULTS_DIR}/${BASE}.erp_hybrid.plink.out;

    # # 1,4,10 - SNP_ID, effective_allele, effective size.
    # ${PLINK_DIR}/plink \
    #     --bfile ${RESULTS_DIR}/bed/Merge \
    #     --score ${BASE_FILE} 1 4 10 header \
    #     --out ${RESULTS_DIR}/${BASE}.ern.plink.out;

    # # 1,4,12 - SNP_ID, effective_allele, effective size.
    # ${PLINK_DIR}/plink \
    #     --bfile ${RESULTS_DIR}/bed/Merge \
    #     --score ${BASE_FILE} 1 4 12 header \
    #     --out ${RESULTS_DIR}/${BASE}.ern_hybrid.plink.out;
done