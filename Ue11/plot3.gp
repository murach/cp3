#!/usr/bin/gnuplot
# kate: encoding utf8
set terminal png lw 2 size 800, 800
set size ratio 1
#set encoding utf8
set grid xtics mytics ytics
set key box
set key left top
#set mxtics 5
#set mytics 5
#set logscale x
#set logscale y
#set pointsize 1
set xlabel "L"
#set xrange [20:100]
set fit errorvariables

set ylabel "<|M|^2>"
set output "M2_kappa1_3d.png"
f(x) = chi*x**(-3)
fit f(x) 'output.dat' using 1:2:3 via chi
f2(x) = M02
fit f2(x) 'output.dat' using 1:2:3 via M02
set title sprintf('chi: %6.4f +- %6.4f, |M0|: %6.4f +- %6.4f', chi, chi_err, sqrt(M02), 1/sqrt(M02)*0.5*M02_err)
plot 'output.dat' using 1:2:3 w errorbars, f(x), f2(x)

set ylabel "U4"
set output "U4_kappa1_3d.png"
plot 'output.dat' using 1:4 w lp