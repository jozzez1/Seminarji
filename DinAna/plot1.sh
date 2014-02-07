#!/bin/sh

tau=$1
theta=$2
k1=$3
k2=$4

filename="zoom-t$1-q$2-k$3-K$4"
output=$filename.png
input=$filename.txt

echo "reset;
	set term pngcairo enhanced truecolor size 1600, 1400 fontscale 3.5;
	set tics scale 2.5;
	set tics out;
	set tics nomirror;
	set encoding utf8;
	set title '{/Symbol t} = $tau, {/Symbol q} = $theta, k = $k1, k\´ = $k2';
	set output \"$output\";
	set xrange [0.9*pi:1.1*pi];
	set xtics ('0.90{/Symbol p}' 0.9*pi, \
		'0.95{/Symbol p}' 0.95*pi, \
		'{/Symbol p}' pi, \
		'1.05{/Symbol p}' 1.05*pi, \
		'1.10{/Symbol p}' 1.10*pi);
	set ytics 0.1
	set xlabel \"x(t)\";
	set ylabel \"p(t)\";
	plot \"$input\" u 1:2 w p pt 7 ps 0.1 lt rgb \"\#6599E5\" title '';
	unset output" | gnuplot

echo "Done!"

convert $filename.png -resize 800x700 $filename-res.png
feh $filename.png &

exit 0

