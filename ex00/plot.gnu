set xlabel 'x'
set ylabel 'u'
p 'sol.dat' t 'Numerical' w lp lw 2, \
  x+sin(2*pi*x) t 'Exact' w l lw 2
