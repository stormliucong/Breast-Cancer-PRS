RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw";
TARGET_DIR=${RAW_DIR}/genomic-data;
PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink";
RESULTS_DIR='/home/cl3720/2019-PRS-pipeline/results/breast_cancer';
SCRIPT_DIR="/home/cl3720/2019-PRS-pipeline/script";
SOFTWARE_DIR="/home/cl3720/2019-PRS-pipeline/software";

cd ${RESULTS_DIR};
PHENOTYPE_FILE=${RESULTS_DIR}/'phenotype.AllBreastCancer'; # different for primary cancer and seconadry cancer




for KEEP_FILE in ${RESULTS_DIR}/keep_*.txt;do
    KEEP=$(basename ${KEEP_FILE} .txt);
    for BASE_FILE in ${RESULTS_DIR}/PRS_*.txt;do
        BASE=$(basename ${BASE_FILE} .txt);
        less ${BASE_FILE} | awk '{print $1}' > ${RESULTS_DIR}/${BASE}.valid.snp;
        less ${BASE_FILE} | awk 'NR>1 {print $2,$3,$3,$1}' > ${RESULTS_DIR}/${BASE}.valid.snprange;

        cd ${RESULTS_DIR}/bed/;
        rm *
        for chr in $(seq 1 22); do
        ${PLINK_DIR}/plink \
            --bfile ${TARGET_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged.exclude.missnp \
            --extract ${RESULTS_DIR}/${BASE}.valid.snprange --range \
            --make-bed \
            --out ${RESULTS_DIR}/bed/chr${chr}.dose.emerge_ids.consented.merged.exclude.missnp.bc ;
        done

        find . -name "*.missnp.bc.bim" > ForMerge.list ;
        sed -i 's/.bim//g' ForMerge.list ;
        ${PLINK_DIR}/plink --merge-list ForMerge.list --out ${RESULTS_DIR}/bed/Merge;

        TARGET_FILE=${RESULTS_DIR}/bed/Merge;


        Rscript ${SOFTWARE_DIR}/PRSice/PRSice.R --dir ${SOFTWARE_DIR}/PRSice \
            --prsice ${SOFTWARE_DIR}/PRSice/PRSice_linux \
            --A1 effect_allele \
            --A2 reference_allele \
            --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
            --base ${BASE_FILE} \
            --binary-target T \
            --bp BP \
            --chr CHR \
            --clump-kb 250 \
            --clump-p 1.000000 \
            --clump-r2 0.100000 \
            --interval 5e-05 \
            --missing MEAN_IMPUTE \
            --model add \
            --beta \
            --out ${RESULTS_DIR}/${BASE}.${KEEP}.PRSice.out \
            --score avg \
            --seed 2307739121 \
            --snp SNP_NEW \
            --stat beta \
            --target ${TARGET_FILE} \
            --thread 10 \
            --upper 0.5 \
            --print-snp \
            --keep-ambig \
            --keep ${KEEP_FILE} \
            --pheno ${PHENOTYPE_FILE} ;
        
        if [ -f ${RESULTS_DIR}/${BASE}.${KEEP}.PRSice.out.valid ]; then
            Rscript ${SOFTWARE_DIR}/PRSice/PRSice.R --dir ${SOFTWARE_DIR}/PRSice \
            --prsice ${SOFTWARE_DIR}/PRSice/PRSice_linux \
            --A1 effect_allele \
            --A2 reference_allele \
            --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
            --base ${BASE_FILE} \
            --binary-target T \
            --bp BP \
            --chr CHR \
            --clump-kb 250 \
            --clump-p 1.000000 \
            --clump-r2 0.100000 \
            --interval 5e-05 \
            --missing MEAN_IMPUTE \
            --model add \
            --beta \
            --out ${RESULTS_DIR}/${BASE}.${KEEP}.PRSice.out \
            --score avg \
            --seed 2307739121 \
            --snp SNP_NEW \
            --stat beta \
            --target ${TARGET_FILE} \
            --thread 10 \
            --upper 0.5 \
            --print-snp \
            --keep-ambig \
            --keep ${KEEP_FILE} \
            --pheno ${PHENOTYPE_FILE} \
            --extract ${RESULTS_DIR}/${BASE}.${KEEP}.PRSice.out.valid;
        else
            echo "PRS pipeline finished"
        fi
        
    done;
done;
