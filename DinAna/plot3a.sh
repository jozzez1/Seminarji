#!/bin/sh

gnuplot << TOEND

set xlabel '\$t_{\\text{max}}\$'

set term epslatex color solid size 16cm,7cm
set xrange [0:50]
set tics out
set tics nomirror

set output "new-corr-1.tex"
set multiplot layout 1,2 title '\v Casovni korelaciji, \$\\tau = \\sqrt{2}-1,\\ k=0.1,\\ k^\\prime = 0.15\$'

set key top left
set title '\$C^{(x)}_t\$'
set ylabel '\$\\langle x_0 x_t \\rangle - \\langle x_0 \\rangle\\langle x_t\\rangle\$'
plot "corr-0.41-0.00-0.10-0.15.txt" u 1:3 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.10-0.15.txt" u 1:3 w l lt 3 title '\$\\theta = 0.25\$'

set title '\$C^{(p)}_t\$'
set ylabel '\$\\langle p_0 p_t \\rangle - \\langle p_0 \\rangle\\langle p_t\\rangle\$'
plot "corr-0.41-0.00-0.10-0.15.txt" u 1:4 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.10-0.15.txt" u 1:4 w l lt 3 title '\$\\theta = 0.25\$'

unset multiplot
unset output

set output "new-corr-2.tex"
set multiplot layout 1,2 title '\v Casovni korelaciji, \$\\tau = \\sqrt{2}-1,\\ k=0.5,\\ k^\\prime = 0.3\$'

set key top left
set title '\$C^{(x)}_t\$'
set ylabel '\$\\langle x_0 x_t \\rangle - \\langle x_0 \\rangle\\langle x_t\\rangle\$'
plot "corr-0.41-0.00-0.50-0.30.txt" u 1:3 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.50-0.30.txt" u 1:3 w l lt 3 title '\$\\theta = 0.25\$'

set title '\$C^{(p)}_t\$'
set ylabel '\$\\langle p_0 p_t \\rangle - \\langle p_0 \\rangle\\langle p_t\\rangle\$'
plot "corr-0.41-0.00-0.50-0.30.txt" u 1:4 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.50-0.30.txt" u 1:4 w l lt 3 title '\$\\theta = 0.25\$'

unset multiplot
unset output

set output "new-corr-3.tex"
set multiplot layout 1,2 title '\v Casovni korelaciji, \$\\tau = \\sqrt{2}-1,\\ k=0.9,\\ k^\\prime = 10^{-6}\$'

set key top left
set title '\$C^{(x)}_t\$'
set ylabel '\$\\langle x_0 x_t \\rangle - \\langle x_0 \\rangle\\langle x_t\\rangle\$'
plot "corr-0.41-0.00-0.90-0.00.txt" u 1:3 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.90-0.00.txt" u 1:3 w l lt 3 title '\$\\theta = 0.25\$'

set title '\$C^{(p)}_t\$'
set ylabel '\$\\langle p_0 p_t \\rangle - \\langle p_0 \\rangle\\langle p_t\\rangle\$'
plot "corr-0.41-0.00-0.90-0.00.txt" u 1:4 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.90-0.00.txt" u 1:4 w l lt 3 title '\$\\theta = 0.25\$'

unset multiplot
unset output

set output "new-corr-4.tex"
set multiplot layout 1,2 title '\v Casovni korelaciji, \$\\tau = \\sqrt{2}-1,\\ k = 10^{-6},\\ k^\\prime = 0.9\$'

set key top left
set title '\$C^{(x)}_t\$'
set ylabel '\$\\langle x_0 x_t \\rangle - \\langle x_0 \\rangle\\langle x_t\\rangle\$'
plot "corr-0.41-0.00-0.00-0.90.txt" u 1:3 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.00-0.90.txt" u 1:3 w l lt 3 title '\$\\theta = 0.25\$'

set title '\$C^{(p)}_t\$'
set ylabel '\$\\langle p_0 p_t \\rangle - \\langle p_0 \\rangle\\langle p_t\\rangle\$'
plot "corr-0.41-0.00-0.00-0.90.txt" u 1:4 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.00-0.90.txt" u 1:4 w l lt 3 title '\$\\theta = 0.25\$'

unset multiplot
unset output

TOEND
