reset
set table 'color.dat'
splot 'u.dat'
unset table

set contour base
unset surface
set cntrparam levels 10
set table 'table.dat'
splot 'u.dat'
unset table

reset
set term x11
unset key
set palette rgbformulae 33,13,10
set xran[0:1]
set yran[0:1]
set size square
plot 'color.dat' with image, 'table.dat' u 1:2 w l lt -1 lw 1.5

set term postscript enhanced color
set out 'contour.eps'
replot
set out
set term x11
