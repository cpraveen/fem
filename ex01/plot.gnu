set xlabel 'x'
set ylabel 'u'
p 'solution.gnuplot' t 'FEM' w l lw 2,\
  x+sin(4*pi*x) t 'Exact'
