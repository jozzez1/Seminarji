reset

set term epslatex color solid size 12cm,10cm
set output "sila.tex"
set log xy
set xrange[0.0001:1]
set xlabel '$x$'
set ylabel '$\mathscr{F}(x)$'
set grid
set title 'Spreminjanje sile z raztezkom

f(x) = x - 0.25 + 0.25/((1 - x)**2)

plot f(x) title ''
