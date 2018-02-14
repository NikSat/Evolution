#!/bin/bash

for i in {1..10}
do 
./SimpIndexHaplotypes.pl -dir $i $1 -number_loci_mhc 1 -loci_kir 5
done

echo "Goodbye"
