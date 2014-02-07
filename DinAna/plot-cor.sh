#!/bin/sh

tau=$1
theta=$2
k1=$3
k2=$4

filename="corr-$1-$2-$3-$4"
output=$filename.png
input=$filename.txt

gnuplot -p << TOEND

reset
set term pngcairo enhanced color size 1000,1600

set output "$output"

set tics out
set tics nomirror

set multiplot layout 2,1

plot "$input" u 1:3 w l title 'C<x>'

plot "$input" u 1:4 w l title 'C<p>'

unset multiplot
unset output

TOEND

echo "Done!"

feh $output &

exit 0

