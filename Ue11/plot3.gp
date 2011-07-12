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

set ylabel "<|M|^2>"
set output "M2_kappa1_3d.png"
plot 'output.dat' using 1:2 w lp

set ylabel "U4"
set output "U4_kappa1_3d.png"
plot 'output.dat' using 1:4 w lp