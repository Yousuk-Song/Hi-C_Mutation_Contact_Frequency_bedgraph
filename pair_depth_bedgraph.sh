#!/usr/bin/bash

HNT_DIR=/data2/home/ujkim/HNT_pipeline
sort -k1,2n ${1} > ${1}.sort.txt
wait
python ${HNT_DIR}/to_bedgraph.py ${1}.sort.txt ${2}
wait
rm ${1}.sort.txt

