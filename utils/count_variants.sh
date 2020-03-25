a=0
PWD='/home/cl3720/2019-PRS-pipeline/raw/genomic-data/bed'
for file in ${PWD}/*.bim;do
    tmp=`less ${file} | wc -l`;
    echo $tmp;

done
