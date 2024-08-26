reset
unset key
set grid
set yran[-1.1:1.1]
filename(n) = sprintf("sol_%d.gpl",n)
avgfilename(n) = sprintf("avg_%d.gpl",n)
N=system("ls -1 sol*.gpl | wc -l")
do for [i=0:N-1] {
   plot filename(i) u 1:2 w l lw 2,avgfilename(i) u 1:2 w p pt 6
   pause 0.2
}
