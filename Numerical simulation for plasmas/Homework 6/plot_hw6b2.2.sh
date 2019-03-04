# plot file for hw6b.1

set style data lines
set logscale xy

plot 'hw6_errors.dat' u 1:2

set xlabel "N"
set ylabel "Error"
set term png
set output "hw6b_2.2.png"
replot
set term x11
