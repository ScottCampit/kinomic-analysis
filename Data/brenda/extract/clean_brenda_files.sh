#!/bin/bash
########################## CLEAN BRENDA FILES ##################################
# kcat
sed -i 's/\r$//' kcat.txt &
sed -i 's/ecNumber/\n&/g' kcat.txt &
sed -i 's/[*]/\t/g' kcat.txt &
sed -i 's/[#]/\t/g' kcat.txt &
sed -i 's/[!]//g' kcat.txt &

# Km
sed -i 's/\r$//' km.txt &
sed -i 's/ecNumber/\n&/g' km.txt &
sed -i 's/[*]/\t/g' km.txt &
sed -i 's/[#]/\t/g' km.txt &
sed -i 's/[!]//g' km.txt &

# Ki
sed -i 's/\r$//' ki.txt &
sed -i 's/ecNumber/\n&/g' ki.txt &
sed -i 's/[*]/\t/g' ki.txt &
sed -i 's/[#]/\t/g' ki.txt &
sed -i 's/[!]//g' ki.txt &

# kcat/km
sed -i 's/\r$//' kcat_km.txt &
sed -i 's/ecNumber/\n&/g' kcat_km.txt &
sed -i 's/[*]/\t/g' kcat_km.txt &
sed -i 's/[#]/\t/g' kcat_km.txt &
sed -i 's/[!]//g' kcat_km.txt &

# ptms
sed -i 's/\r$//' ptms.txt &
sed -i 's/ecNumber/\n&/g' ptms.txt &
sed -i 's/[*]/\t/g' ptms.txt &
sed -i 's/[#]/\t/g' ptms.txt &
sed -i 's/[!]//g' ptms.txt &

# Substrates and Products
sed -i 's/\r$//' sub_pdt.txt &
sed -i 's/ecNumber/\n&/g' sub_pdt.txt &
sed -i 's/[*]/\t/g' sub_pdt.txt &
sed -i 's/[#]/\t/g' sub_pdt.txt &
sed -i 's/[!]//g' sub_pdt.txt

# EC Numbers
sed -i 's/\r$//' ec.txt &
sed -i 's/ecNumber/\n&/g' ec.txt &
sed -i 's/[*]/\t/g' ec.txt &
sed -i 's/[#]/\t/g' ec.txt &
sed -i 's/[!]//g' ec.txt
