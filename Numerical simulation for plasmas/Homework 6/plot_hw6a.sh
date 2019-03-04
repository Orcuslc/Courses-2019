# plot file for hydro

set style data lines
plot 'sample1.out.0000' u 4:7
replot 'sample1.out.0002' u 4:7
replot 'sample1.out.0004' u 4:7
replot 'sample1.out.0006' u 4:7
replot 'sample1.out.0008' u 4:7
replot 'sample1.out.0010' u 4:7

set xlabel "x"
set ylabel "u_x"
set term png
set output "sample.png"
replot
set term x11
