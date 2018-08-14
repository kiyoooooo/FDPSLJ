reset
set terminal pdfcairo color enhanced
set output "DrawEnergy.pdf"
set tics font "Arial,15"   # 目盛りのフォントの変更
set xlabel font "Arial,15"
set ylabel font "Arial,15"
set xlabel "step"
set ylabel "Energy"
set format x "%3.1f"
set format y "%3.1f"
#set xrange [0:7]
#et yrange [0:3]
#set key font "Arial,24"
#unset key # 凡例を表示しない

#set key right bottom
set xtics rotate by -90

plot "DrawEnergy.dat" using 1:2 w l lw 2 title "total",\
"DrawEnergy.dat" using 1:3 w l lw 2 title "kinetic",\
"DrawEnergy.dat" using 1:4 w l lw 2 title "potential"

set terminal aqua
set output