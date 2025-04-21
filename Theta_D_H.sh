#!/bin/bash

k=$1

CHROM_LENGTHS=(248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189 90338345 83257441 80373285 58617616 64444167 46709983 50818468)
echo -e $k"\t"chr$k"\t"1"\t"${CHROM_LENGTHS[$(($k-1))]} >region.chr$k.txt

for group in `ls -d ../00.input/*/|awk -F '/' '{print $(NF-1)}'`
do
    mkdir -p $group
    vcf="../00.input/$group/${group}.chr${k}.vcf.gz"
    ./Theta_D_H.Est \
        --gzvcf $vcf \
        --region region.chr$k.txt \
        --samples ../00.input/$group/$group.list \
        --window_shift 50000@10000 \
        --out $group/$group.chr$k
done
