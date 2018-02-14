#!/bin/bash


clear

echo "Processing Population Files"

./PAverage.pl -dir PepPresentation
./Rprep.pl -dir PepPresentation

echo "Plotting"

R CMD BATCH plotPop.R
gnuplot44 plotDynamics.gnu

echo "Calculating Licensing Indexes"

for i in {1..10}; do ./AvHapLicense.pl -dir PepPresentation -number "$i"; done

for i in {1..10}; do ./AvKIRLicense.pl -dir PepPresentation -number "$i"; done

echo "Plotting"

gnuplot44 plotHap.gnu

gnuplot44 plotKIR.gnu

echo "Calculating Diversity Indexes"

for i in {1..10}; do ./SICalc.pl -dir PepPresentation -kir -number "$i"; done
./SSICalc.pl -dir PepPresentation -years 150000

for i in {1..10}; do ./SICalc.pl -dir PepPresentation -hap -number "$i"; done
./SSICalc.pl -dir PepPresentation -years 150000 -hap

./RSIprep.pl -dir PepPresentation -kir

./RSIprep.pl -dir PepPresentation -hap

echo "Plotting"

gnuplot44 plotSRIkir.gnu
gnuplot44 plotSRIhap.gnu
R CMD BATCH plotSRIhap.R
R CMD BATCH plotSRI.R

echo "Calculating Diversity Indexes"

for i in {1..10}; do ./PresCalc.pl -dir PepPresentation -number "$i"; done

echo "Plotting"

gnuplot44 plotSelf.gnu
gnuplot44 plotNonSelf.gnu

echo "Estimating Haplotype Contraction"

for i in {1..10}; do ./ContractionCalc.pl -dir PepPresentation -number "$i"; done

echo "Plotting"

gnuplot44 plotCON.gnu

echo "Finished"
