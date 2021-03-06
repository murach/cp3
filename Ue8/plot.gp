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
set xlabel "B"
#set xrange [20:100]

set ylabel "M"
set output "mag.png"
plot 'hysterese.dat' using 1:2:3 with yerrorbars

set ylabel "M"
set output "mag2.png"
plot 'hysterese.dat' using 1:2 w lp