#!/usr/bin/bash

Pipeline_DIR=
sort -k1,2n ${1} > ${1}.sort.txt
wait
python ${Pipeline_DIR}/to_bedgraph.py ${1}.sort.txt ${2}
wait
rm ${1}.sort.txt

