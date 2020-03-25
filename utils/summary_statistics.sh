RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/"
BASE_DIR=${RAW_DIR}/base/


less ${BASE_DIR}/breast_cancer/Mavaddat_supp_table_7.txt | \
    awk '{print $2":"$3,$2,$3,$5,$4,$7,"0"}' | \
    sed '1d' | \
    sed '1iSNP CHR BP A1 A2 BETA P' | \
    wc -l

less ${BASE_DIR}/breast_cancer/Mavaddat_supp_table_7.txt | \
    awk '{print $2":"$3,$2,$3,$5,$4,$7,"0"}' | \
    sed '1d' | \
    sed '1iSNP CHR BP A1 A2 BETA P' \
    > Mavaddat_2019_SNPS-313.assoc




