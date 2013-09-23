#!/bin/sh

echo "reset;
	set term epslatex color solid;
	set output 'graf1.tex'
	set xlabel '\$\omega\$';
	set ylabel '\$\varepsilon - 1\$';
	set title '\$\beta = 12\$';
	set grid;
	plot \"graf1.txt\" u 1:3 w l title '\$\varepsilon - 1\$';" | gnuplot -p

epstopdf graf1.eps

echo "reset;
	set term epslatex color solid;
	set output 'graf2.tex';
	set xlabel '\$\omega\$';
	set ylabel '\$\beta\$';
	set zlabel '\$\varepsilon - 1\$';
	set ztics 0.5
	splot \"graf2.txt\" w l title '\$\varepsilon - 1\$';" | gnuplot -p

epstopdf graf2.eps

echo "reset;
	set term epslatex color solid;
	set output 'graf3.tex';
	set xlabel '\$\beta\$';
	set ylabel '\$\varepsilon - 1\$';
	set log x;
	set title '\$\omega = 0\$';
	set grid;
	plot \"graf3.txt\" u 1:3 w l title '\$\varepsilon - 1\$';" | gnuplot -p

epstopdf graf3.eps

pdflatex dn-joze.tex

exit 0;

