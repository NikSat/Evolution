
data2<-read.table("RSIprephapPepPresentation.txt")

png("SIhapBox.png")

lines(boxplot(data2$V2~data2$V1, col=4, xlab="Percentage of MHC-KIR recognition ", ylab="Simpson's Reciprocal Index", xaxt='n') $stats[c(3),])

axis (1, at=1:10, lab=c(100, 96, 70, 40, 20, 9.8, 4.1, 2.4, 0.1, 0.06))

dev.off()

