#!/bin/bash


for i in {1..10}
do
echo "Processing directory" $i
./CheckEvol.pl -dir $i -number_loci_mhc 1 -loci_kir 5
echo "Average (non)self peptide recognition Ok"
./Contraction.pl -dir $i -number_loci_mhc 1 -loci_kir 5
echo "Haplotype contraction checked"
./PepHaLicense.pl -dir $i -loci_kir 5 -number_loci_mhc 1
echo "Haplotype average licensing Ok"
./PepKIRLicense.pl -dir $i -loci_kir 5 -number_loci_mhc 1
echo "KIR average licensing Ok"
echo "Done"
done

