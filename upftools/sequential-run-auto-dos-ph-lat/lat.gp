#!/usr/bin/gnuplot -persist

# plot
set terminal png
set output "lat.png"
set xlabel "X"
set ylabel "Lattice constant, a / nm"
set xrange [0:1]
set yzeroaxis
plot "lat.dat" u 1:($2/1.4142) w p, "lat_line.dat" u 1:($2/1.4142) w l
