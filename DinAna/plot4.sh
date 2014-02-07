#!/bin/sh

gnuplot << TOEND
# funkcija, ki smo jo dobili z linearno regresijo
#########################
set term epslatex color solid
set output "ljp-determ.tex"
f(x) = 1.023596753407e+00 * x + -1.298349919359e+02

set title 'Iskanje Ljapunovega eksponenta za marginalno kaoti\\v cen sistem'
set xlabel '\$t_i\$'
set ylabel '\$\\sum_{t_i < t} y_i\$'

set key left
set grid
set label '\$\\lambda \\approx 1.02,\\ C \\approx -129.83\$' at 400,10

plot "ljp-determ.txt" u 1:2 w l title 'eksperiment', f(x) lt -1 title 'regresijska premica'
unset output

TOEND

exit 0
