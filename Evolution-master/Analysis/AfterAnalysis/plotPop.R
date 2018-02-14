
data<-read.table("RprepTotal.txt")

png("PopBox.png")

lines( boxplot(data$V2~data$V1, col=2, xlab="Percentage of MHC-KIR recognition ", ylab="Population Size", xaxt='n', title="Total Population") $stats[c(3),])

axis (1, at=1:10, lab=c(100, 96, 70, 40, 20, 9.8, 4.1, 2.4, 0.1, 0.06))

dev.off()
