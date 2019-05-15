#!/bin/bash

# Clean up uniprot accession
sed -i 's/"\t\t\t\t\t\t\t\t\t\t"\/n&/g' uniprot_kinetics.tab
