reset

set term epslatex color solid size 12cm,10cm
set output "sila1.tex"
set xrange[0.0001:1]
set xlabel '$x$'
set ylabel '$\rho^2(x)$'
set grid
set title 'Spreminjanje sile z raztezkom'

f(x) = x + exp((-1)*x)

plot f(x) title ''
