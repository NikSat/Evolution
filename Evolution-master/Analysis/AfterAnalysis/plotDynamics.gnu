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
set output "PopulationDynamics.png"
set ylabel "Population Size" font ",20"
set xlabel "KIR-MHC recognition" font ",20"
set xtics ("100" 1,"96" 2,"70" 3,"40" 4,"20" 5,"9.8" 6,"4.1" 7,"2.4" 8,"0.1" 9,"0.06" 10)
set xtics font ",15"
set ytics font ",15"

set yr [0:5000]
#set label "Decoy Presenting Virus" font ",20" at 70000,9
plot "AveragePopulation.txt" using 1:2 w l lc 1 linewidth 4 title "Total Population" ,\
"AveragePopulation.txt" using 1:4 w l ls 2 title "Acute Infected" ,\
"AveragePopulation.txt" using 1:6 w l ls 3 title "Chronically Infected" ,\
"AveragePopulation.txt" using 1:8 w l ls 7 title "Immune"

# ../DownRegulation/
#Down Regulating Virus
#AVKIRDownregulation
#set yr [0:10]
#set label "Down Regulating Virus" font ",20" at 70000,15
#set yr [2:16]
