reset
unset key
set grid
set size ratio -1
filename(n) = sprintf("grid-%d.gnuplot",n)
title(n) = sprintf("Grid level = %d",n)
outfile(n) = sprintf("grid-%d.svg",n)
N=system("ls -1 grid-*.gnuplot | wc -l")
do for [i=0:N-1] {
   set title title(i)
   plot filename(i) u 1:2 w l
   pause 2.0
   set term push
   set term svg
   set out outfile(i)
   replot
   set term pop
}
