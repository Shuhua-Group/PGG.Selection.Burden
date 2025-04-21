#!/bin/bash

k=$1
SweeD_d='/home/sunyumeng/software/SweeD_v3.2.1_Linux'

## step1
CHROM_LENGTHS=(248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189 90338345 83257441 80373285 58617616 64444167 46709983 50818468)
grid=`echo ${CHROM_LENGTHS[$(($k-1))]} 50000 |awk '{print int($1/$2+0.5)}'`

cwd=`pwd`
## step2
for group in `ls -d ../00.input/*/|awk -F '/' '{print $(NF-1)}'`
do
  mkdir -p $cwd/$group; cd $cwd/$group
  vcf="../../00.input/$group/${group}.chr${k}.vcf"
  ${SweeD_d}/SweeD -name $group.chr${k}.50kb  \
    -input $vcf \
    -grid ${grid}
done