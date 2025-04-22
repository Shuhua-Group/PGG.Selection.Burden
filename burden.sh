#!/bin/bash

k=$1

cwd=`pwd`
for group in `ls -d ../00.input/*/|awk -F '/' '{print $(NF-1)}'`
do
    mkdir -p $cwd/$group;cd $cwd/$group
    vcf="../../00.input/$group/${group}.chr${k}.vcf.gz"
    /usr/bin/python3.6 ../burden.py $vcf ../chr${k}.DamageSnp.txt ../kegg.input.txt $group.chr$k 1>${group}.chr${k}.log 2>&1
done 