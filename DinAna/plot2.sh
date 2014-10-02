#!/bin/sh

tau=$1
theta=$2
k1=$3
k2=$4

filename="new-t$1-q$2-k$3-K$4"
output=$filename.png
input=$filename.txt

echo "reset;
	set term pngcairo enhanced truecolor size 4200, 2400 fontscale 3.5;
	set output \"$filename.png\";
	set tics scale 2.5;
	set tics out;
	set tics nomirror;
	set encoding utf8;
	set title '{/Symbol t} = $tau, {/Symbol q} = $theta, k = $k1, k\´ = $k2';
	set xrange [0:2*pi];
	set yrange [-5:5];
	set xlabel \"x(t)\";
	set xtics ('0' 0, '{/Symbol p}/2' pi/2, '{/Symbol p}' pi, \
		'3{/Symbol p}/2' 3*pi/2, '2{/Symbol p}' 2*pi);
	set ylabel \"p(t)\";
	set multiplot layout 1,3;
	plot \"$filename.txt\" u 1:2 w p pt 7 ps 0.1 lt rgb \"\#6599E5\" title '';
	plot \"$filename-1.txt\" u 1:2 w p pt 7 ps 0.1 lt rgb \"\#00B389\" title '';
	plot \"$filename-2.txt\" u 1:2 w p pt 7 ps 0.1 lt rgb \"\#FC2BD4\" title '';
	unset multiplot;
	unset output" | gnuplot

echo "Done!"

convert $filename.png -resize 1800x900 $filename-res.png
feh $filename.png &

exit 0
