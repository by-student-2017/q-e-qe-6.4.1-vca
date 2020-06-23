#!/usr/bin/gnuplot -persist

# get fermi energy
stats "ef.txt" nooutput
ef = STATS_mean_x

# plot
set terminal png
set output "band.png"
set yrange [-14:6]
set ylabel "Energy / eV"
set xzeroaxis
set noxtics
plot "bands.out.gnu" u 1:($2-ef) w l t ""
