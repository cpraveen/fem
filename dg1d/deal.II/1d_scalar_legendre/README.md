# 1-d linear scalar PDE using DG

To solve constant linear advection equation

```shell
cp linadv/*.h .
```

To solve Burgers equation

```shell
cp burgers/*.h .
```

Create Makefile (this only needs to be done once)

```shell
cmake .
make release
```

Compile

```shell
make
```

Set parameters in file `input.prm` and run

```shell
rm -f *.gpl
./dg input.prm
gnuplot anim.gnu
```

If you modify `dg.cc`, then you must compile and run again.

## Grid convergence test

Set `nrefine > 1` in `input.prm` and run the code to perform error convergence study. A tex file `error.tex` with convergence table will also be generated. You can compile and see the pdf.

```shell
pdflatex error.tex
open error.pdf
```
