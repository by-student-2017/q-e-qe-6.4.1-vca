#!/usr/bin/gnuplot -persist

# get fermi energy
stats "ef.txt" nooutput
ef = STATS_mean_x

# plot
set terminal png
set output "tdos.png"
set xlabel "Energy / eV"
set ylabel "Density of State / eV"
set xrange [-14:6]
set yzeroaxis
plot "pwscf.pdos_tot" u ($1-ef):2 w l t "tdos"
