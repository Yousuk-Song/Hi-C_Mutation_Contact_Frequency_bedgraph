#!/usr/bin/bash

./rename.py sample_ID_list.csv 
wait
./filter.py sample_ID_list.csv.output.txt pointmutation.csv
