#!/bin/bash
./brenda_ec.pl > ec.txt;
./brenda_enzymenam.pl > enzymenam.txt;
./brenda_ic50.pl > ic50.txt;
./brenda_kcat_km_val.pl > kcat_km.txt;
./brenda_kcat.pl > kcat.txt;
./brenda_km.pl > km.txt;
./brenda_ki.pl > ki.txt;
./brenda_pathways.pl > pathways.txt;
./brenda_pdb.pl > pdb.txt;
./brenda_sub_pdt.pl > sub_pdt.txt;
./brenda_ptms.pl > ptm.txt
./brenda_reactions.pl > reactopms.txt
#cwd = '.'
#for file in "$dir"/*.txt; do
#    while read LINE; do
#      sed -i 's/#ecNumber/\n/g'
#  done
