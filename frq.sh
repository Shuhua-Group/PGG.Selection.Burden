#!/bin/bash

k=$1

for group in `ls -d ../00.input/*/|awk -F '/' '{print $(NF-1)}'`
do
    mkdir -p $group
    vcf="../00.input/$group/${group}.chr${k}.vcf.gz"
    plink2 --vcf $vcf --double-id --set-all-var-ids @:# --freq --out $group/$group.chr$k
done