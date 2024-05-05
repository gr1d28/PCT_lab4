#!/usr/bin/gnuplot
set termoption enhanced
set terminal svg size 1300,600 font "Arial, 16"
set output "N-body.svg"
set style line 1 lc rgb "0xDC143C" lt 1 lw 4 pt 9 ps 1
set style line 2 lc rgb "green" lt 1 lw 4 pt 9 ps 1
set style line 3 lc rgb "blue" lt 1 lw 4 pt 9 ps 1
set style line 4 lc rgb "black" lt 1 lw 4 pt 9 ps 1
set style line 5 lc rgb "yellow" lt 1 lw 4 pt 9 ps 1
set style line 6 lc rgb "violet" lt 1 lw 4 pt 9 ps 1
set border lw 2
set grid
set key top left
set xlabel "Количество потоков"
set ylabel "Коэффициент ускорения" rotate by 90
set xtics 1
set mxtics
set xrange [0:8]
set yrange [0:8]
set format x "%6.0f"
set format y "%.6f"
plot "line.dat" using 1:2 title "line" with linespoints ls 1 , \
"critical.dat" using 1:2 title "critical" with linespoints ls 2, \
"atomic.dat" using 1:2 title "atomic" with linespoints ls 3, \
"n_blocked.dat" using 1:2 title "n blocked" with linespoints ls 4, \
"redundant_calc.dat" using 1:2 title "redundant calc" with linespoints ls 5, \
"add_memory.dat" using 1:2 title "add memory" with linespoints ls 6