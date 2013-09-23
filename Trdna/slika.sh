#!/bin/sh

echo 'splot "output.txt" u 1:2:3 w l title "rezultat"' | gnuplot -p

exit 0
