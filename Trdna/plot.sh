#!/bin/sh

echo "reset;
	set term epslatex color solid;
	set output 'graf1.tex'
	set xlabel '\$\omega\$';
	set ylabel '\$\varepsilon\$';
	set title '\$\beta = 14\$';
	set log y;
	plot \"graf1.txt\" u 1:(abs(\$3)) w l title '\$\varepsilon\$';" | gnuplot -p

epstopdf graf1.eps

echo "reset;
	set term epslatex color solid;
	set output 'graf2.tex';
	set view 62,31;
	set xlabel '\$\omega\$';
	set ylabel '\$\beta\$';
	set zlabel '\$\varepsilon\$';
	splot \"graf2.txt\" w l title '\$\varepsilon\$';" | gnuplot -p

epstopdf graf2.eps

exit 0;

