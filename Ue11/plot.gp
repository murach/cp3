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
set output "M2_kappa1_2d.png"
f(x) = chi*x**(-2)
fit f(x) 'output.dat' using 1:2:3 via chi
f2(x) = a*x**(-eta)
fit f2(x) 'output.dat' using 1:2:3 via a,eta
set title sprintf('chi: %6.4f +- %6.4f, eta: %6.4f +- %6.4f', chi, chi_err, eta, eta_err)
plot 'output.dat' using 1:2:3 w errorbars, f(x), f2(x)

set ylabel "U4"
set output "U4_kappa1_2d.png"
plot 'output.dat' using 1:4 w lp