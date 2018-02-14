#!/bin/bash

clear

#Script that does the over all analysis
#Takes as input the mutant number and the name of the directory
#nikolaos@mutant"$1":/linuxhome/tmp/nikolas/"$2"/ contains 10 directories from 1 to 10 one for each specififity 

echo "Hello $USER"

ssh mutant"$1" <<EOF
cd /linuxhome/tmp/nikolas/"$2"/
rm -f *.pl
rm -f *.sh
cp ~/Documents/essential/PeptideAnalysis/*.pl .
cp ~/Documents/essential/PeptideAnalysis/SIscript.sh .
cp ~/Documents/essential/PeptideAnalysis/Speckscript.sh .
chmod +x SIscript.sh
chmod +x Speckscript.sh
echo "Files copied"
echo "Collecting Population"
./FullPopCalc.pl -years 150000 -out SortedAverage"$1".txt
echo "Finished"
echo "Calculating Diversity indexes"
./SIscript.sh -kir
echo "KIR SRI Ok"
./SIscript.sh -haplotypes
echo "Haplotype SRI Ok"
echo "Calculating Specificity"
./Speckscript.sh
echo "Done"
EOF
