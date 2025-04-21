#!/bin/bash

k=$1
predictGMAP_d='/home/sunyumeng/software/predictGMAP'
selscan_d='/home/sunyumeng/software/selscan-linux-2.0.0'
map_d='/home/sunyumeng/data/PLINK_format_genetic_map'


for group in `ls -d ../00.input/*/|awk -F '/' '{print $(NF-1)}'`
do
    mkdir -p $group
    vcf="../00.input/$group/${group}.chr${k}.vcf.gz"
    if [ ! -e "chr${k}.query" ]
    then
        bcftools view -H $vcf|cut -f2 > chr${k}.query
        ${predictGMAP_d}/src/predictGMAP --ref $map_d/plink.chr${k}.GRCh38.map --query chr${k}.query --out chr${k}.predicted.map --max-gap 3000000000
    fi 1>$group/${group}.chr${k}.log 2>&1
    ${selscan_d}/selscan --ihs \
        --vcf $vcf \
        --map chr${k}.predicted.map \
        --threads 10 \
        --out $group/${group}.chr${k}
    ${selscan_d}/norm --ihs --files $group/${group}.chr${k}.ihs.out --bp-win --winsize 50000 1>>$group/${group}.chr${k}.log 2>&1
done 