#!/bin/sh

gnuplot << TOEND
# ok, let's make this script a bit different :D
set xlabel '\$t_{\\text{max}}\$'

set term epslatex color solid size 16cm,7cm
set xrange [0:50]
set tics out
set tics nomirror
set output "corr.tex"

set multiplot layout 1,2 title '\v Casovni korelaciji, \$\\tau = 1,\\ k=0.3,\\ k^\\prime = 0.6\$'

set key top left
set title '\$C^{(x)}_t\$'
set ylabel '\$\\langle x_0 x_t \\rangle - \\langle x_0 \\rangle\\langle x_t\\rangle\$'
plot "corr-1.00-0.00-0.30-0.60.txt" u 1:3 w l title '\$\\theta = 0\$', \
	"corr-1.00-0.25-0.30-0.60.txt" u 1:3 w l lt 3 title '\$\\theta = 0.25\$'

set title '\$C^{(p)}_t\$'
set ylabel '\$\\langle p_0 p_t \\rangle - \\langle p_0 \\rangle\\langle p_t\\rangle\$'
plot "corr-1.00-0.00-0.30-0.60.txt" u 1:4 w l title '\$\\theta = 0\$', \
	"corr-1.00-0.25-0.30-0.60.txt" u 1:4 w l lt 3 title '\$\\theta = 0.25\$'

unset multiplot
unset output

set output "corr2.tex"
set multiplot layout 1,2 title '\v Casovni korelaciji, \$\\tau = 1,\\ k=0.3,\\ k^\\prime = 5.0\$'

set title '\$C^{(x)}_t\$'
set ylabel '\$\\langle x_0 x_t \\rangle - \\langle x_0 \\rangle\\langle x_t\\rangle\$'
plot "corr-1.00-0.00-0.30-5.00.txt" u 1:3 w l title '\$\\theta = 0\$', \
	"corr-1.00-0.25-0.30-5.00.txt" u 1:3 w l lt 3 title '\$\\theta = 0.25\$'

set title '\$C^{(p)}_t\$'
set ylabel '\$\\langle p_0 p_t \\rangle - \\langle p_0 \\rangle\\langle p_t\\rangle\$'
plot "corr-1.00-0.00-0.30-5.00.txt" u 1:4 w l title '\$\\theta = 0\$', \
	"corr-1.00-0.25-0.30-5.00.txt" u 1:4 w l lt 3 title '\$\\theta = 0.25\$'

unset multiplot
unset output

set output "corr3.tex"
set multiplot layout 1,2 title '\v Casovni korelaciji, \$\\tau = 1.5,\\ k=0.3,\\ k^\\prime = 0.6\$'

set title '\$C^{(x)}_t\$'
set ylabel '\$\\langle x_0 x_t \\rangle - \\langle x_0 \\rangle\\langle x_t\\rangle\$'
plot "corr-1.50-0.23-0.30-0.60.txt" u 1:3 w l lt 3 title '\$\\theta = 0.23\$', \
	"corr-1.50-0.00-0.30-0.60.txt" u 1:3 w l lt 1 title '\$\\theta = 0\$'

set title '\$C^{(p)}_t\$'
set ylabel '\$\\langle p_0 p_t \\rangle - \\langle p_0 \\rangle\\langle p_t\\rangle\$'
plot "corr-1.50-0.23-0.30-0.60.txt" u 1:4 w l lt 3 title '\$\\theta = 0.23\$', \
	"corr-1.50-0.00-0.30-0.60.txt" u 1:4 w l lt 1 title '\$\\theta = 0\$'

unset multiplot
unset output

set output "corr4.tex"
set multiplot layout 1,2 title '\v Casovni korelaciji, \$\\tau = \\sqrt{2}-1,\\ k=0.3,\\ k^\\prime = 0.6\$'

set title '\$C^{(x)}_t\$'
set ylabel '\$\\langle x_0 x_t \\rangle - \\langle x_0 \\rangle\\langle x_t\\rangle\$'
plot "corr-0.41-0.00-0.30-0.60.txt" u 1:3 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.30-0.60.txt" u 1:3 w l lt 3 title '\$\\theta = 0.25\$'

set title '\$C^{(p)}_t\$'
set ylabel '\$\\langle p_0 p_t \\rangle - \\langle p_0 \\rangle\\langle p_t\\rangle\$'
plot "corr-0.41-0.00-0.30-0.60.txt" u 1:4 w l title '\$\\theta = 0\$', \
	"corr-0.41-0.25-0.30-0.60.txt" u 1:4 w l lt 3 title '\$\\theta = 0.25\$'

unset multiplot
unset output

#set ylabel '\$D_t = \\langle (p_0 - p_t)^2\\rangle / 2t\$'
#set output "diff.tex"
#set xrange [1:10000]
#set log xy
#set multiplot layout 1,2 title 'Difuzija'
#set xtics 10
#
#set title '\$\\theta = 0\$'
#plot "diff-1.00-0.00-0.30-0.60.txt" u 1:(0.5 * \$2/\$1) w l title '\$\\tau = 1\$', \
#	"diff-1.00-0.00-0.30-5.00.txt" u 1:(0.5 * \$2/\$1) w l lt 3 title '\$\\tau = 1, k^\\prime = 5\$', \
#	"diff-1.50-0.00-0.30-0.60.txt" u 1:(0.5 * \$2/\$1) w l lt 4 title '\$\\tau = 1.5\$', \
#	"diff-0.41-0.00-0.30-0.60.txt" u 1:(0.5 * \$2/\$1) w l lt 5 title '\$\\tau = \\sqrt{2}-1\$', \
#	"diff-8.00-0.00-0.30-0.60.txt" u 1:(0.5 * \$2/\$1) w l lt 7 title '\$\\tau = 8\$'
#
#set title '\$\\theta = 0.25\$'
#plot "diff-1.00-0.25-0.30-0.60.txt" u 1:(0.5 * \$2/\$1) w l title '\$\\tau = 1\$', \
#	"diff-1.00-0.25-0.30-5.00.txt" u 1:(0.5 * \$2/\$1) w l lt 3 title '\$\\tau = 1, k^\\prime = 5\$', \
#	"diff-1.50-0.23-0.30-0.60.txt" u 1:(0.5 * \$2/\$1) w l lt 4 title '\$\\tau = 1.5\$', \
#	"diff-0.41-0.25-0.30-0.60.txt" u 1:(0.5 * \$2/\$1) w l lt 5 title '\$\\tau = \\sqrt{2}-1\$', \
#	"diff-8.00-0.23-0.30-0.60.txt" u 1:(0.5 * \$2/\$1) w l lt 7 title '\$\\tau = 8\$'
#
TOEND

exit 0

