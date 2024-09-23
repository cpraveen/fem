reset
set term qt font "Times New Roman,16"
set grid
#set yran[-1.1:1.1]
filename(n) = sprintf("sol_%d.gpl",n)
avgfilename(n) = sprintf("avg_%d.gpl",n)
N=system("ls -1 sol*.gpl | wc -l")
do for [i=0:N-1] {
   plot filename(i) u 1:2 t 'Pressure' w l lw 2, \
        avgfilename(i) u 1:2 t 'Avg Pressure' w p pt 6, \
        filename(i) u 1:3 t 'Velocity' w l lw 3, \
        avgfilename(i) u 1:3 t 'Avg Velocity' w p pt 7
   pause 0.5
}
