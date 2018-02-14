#Dark colours
set style line 1 linecolor rgb '#B0171F' linetype 1 linewidth 4  # Indian Red
set style line 2 linecolor rgb '#0000FF' linetype 1 linewidth 4  # Dark Blue
set style line 3 linecolor rgb '#008000' linetype 1 linewidth 4  # Green
set style line 4 linecolor rgb '#000000' linetype 1 linewidth 4  # Black
set style line 5 linecolor rgb '#8B2252' linetype 1 linewidth 4  # VioletRed4
set style line 6 linecolor rgb '#FF7000' linetype 1 linewidth 4  # Orange
set style line 7 linecolor rgb '#DAA520' linetype 1 linewidth 4  # Golden Rod
set style line 8 linecolor rgb '#4B0052' linetype 1 linewidth 4  # Indigo
set style line 9 linecolor rgb '#008b45' linetype 1 linewidth 4  # SpringGreen3
set style line 10 linecolor rgb '#FF0000' linetype 1 linewidth 4 #pt 13 #Red


 
set key out
set terminal pngcairo size 1200,800
set output "AverageContraction.png"
set ylabel "Average Loci Number" font ",20"
set xlabel "Years" font ",20"
set xtics font ",15"
set ytics font ",15"
set yr [0:6]
#set label "Decoy Presenting Virus" font ",20" at 70000,9
plot "LociAverage_L=1.txt" every 100 using 1:2 w l ls 1 title "Rec (%)=100" ,\
"LociAverage_L=2.txt" every 100 using 1:2 w l ls 2 title "Rec (%)=96" ,\
"LociAverage_L=3.txt" every 100 using 1:2 w l ls 3 title "Rec (%)=70" ,\
"LociAverage_L=4.txt" every 100 using 1:2 w l ls 4 title "Rec (%)=40" ,\
"LociAverage_L=5.txt" every 100 using 1:2 w l ls 5 title "Rec (%)=20" ,\
"LociAverage_L=6.txt" every 100 using 1:2 w l ls 6 title "Rec (%)=9.8" ,\
"LociAverage_L=7.txt" every 100 using 1:2 w l ls 7 title "Rec (%)=4.1" ,\
"LociAverage_L=8.txt" every 100 using 1:2 w l ls 8 title "Rec (%)=2.4" ,\
"LociAverage_L=9.txt" every 100 using 1:2 w l ls 9 title "Rec (%)=0.1" ,\
"LociAverage_L=10.txt" every 100 using 1:2 w l ls 10 title "Rec (%)=0.06" 
# for [i=4:10] "../DownRegulation/AverageHapSpec_L=".i.".txt" every 100 using 1:($2-$3) w l lt 0 lw 2 lc i notitle,\
# for [i=4:10] "../DownRegulation/AverageHapSpec_L=".i.".txt" every 100 using 1:($2+$3) w l lt 0 lw 2 lc i notitle


# ../DownRegulation/
#Down Regulating Virus
#AVKIRDownregulation
#set yr [0:10]
#set label "Down Regulating Virus" font ",20" at 70000,15
#set yr [2:16]
