#!/bin/bash
./brenda_kcat_km_val.pl > kcat_km.txt;
./brenda_kcat.pl > kcat.txt;
./brenda_km.pl > km.txt;
./brenda_ki.pl > ki.txt;
./brenda_sub_pdt.pl > sub_pdt.txt;
./brenda_ptms.pl > ptm.txt

#cwd = '.'
#for file in "$dir"/*.txt; do
#    while read LINE; do
#      sed -i 's/#ecNumber/&\n/g'
#  done
