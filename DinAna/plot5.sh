#!/bin/sh

gnuplot << TOEND

set term epslatex color solid size 10cm,10cm
set output "ljp-scan.tex"
set view map
set xlabel '\$k\$'
set ylabel '\$k^\\prime\$'
set cblabel '\$\\lambda\$
set xrange [0.1:0.975]
set yrange [0.1:0.975]
set tics out
set tics nomirror
set palette rgbformulae 30,31,32
splot "ljp-scan.txt" u 1:2:3 title '' w pm3d

TOEND
