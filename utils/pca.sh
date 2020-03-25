RAW_DIR="/home/cl3720/2019-PRS-pipeline/raw/"
TARGET_DIR=${RAW_DIR}/genomic-data
PLINK_DIR="/home/cl3720/2019-PRS-pipeline/software/plink"
cd ${TARGET_DIR}/pruned
# Get a list of all PLINK files
find . -name "*.bim" > ForMerge.list ;
sed -i 's/.bim//g' ForMerge.list ;

# Merge all projects into a single PLINK file
${PLINK_DIR}/plink --merge-list ForMerge.list --out Merge ;

if test -f "Merge.missnp"; then
   for chr in $(seq 1 22); do
       ${PLINK_DIR}/plink \
           --bfile ${TARGET_DIR}/pruned/chr${chr}.dose.emerge_ids.consented.merged \
           --exclude Merge.missnp \
           --make-bed \
           --out ${TARGET_DIR}/pruned/chr${chr}.dose.emerge_ids.consented.merged.exclude.missnp ;
   done
   find . -name "*.missnp.bim" > ForMerge.exclude.missnp.list ;
   sed -i 's/.bim//g' ForMerge.exclude.missnp.list ;
   ${PLINK_DIR}/plink --merge-list ForMerge.exclude.missnp.list --out Merge ;
fi
${PLINK_DIR}/plink --merge-list ForMerge.exclude.missnp.list --out Merge.exclude.missnp ;
${PLINK_DIR}/plink --bfile Merge.exclude.missnp --pca 10 --out ${TARGET_DIR}/pca/ALL.dose.emerge_ids.consented.merged.prune.in;
