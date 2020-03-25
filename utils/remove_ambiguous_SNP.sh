RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/"
BASE_DIR=${RAW_DIR}/base/breast_cancer
QC_DIR=${BASE_DIR}/qc
cd ${BASE_DIR}

for file in PRS_*.txt;do
    less ${file} |\
    awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' \
    > ${QC_DIR}/${file} ;
done
