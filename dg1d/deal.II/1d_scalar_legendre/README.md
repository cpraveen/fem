# 1-d linear scalar advection using DG

Create Makefile (this only needs to be done once)

```shell
cmake .
make release
```

Compile and run

```shell
make
rm -f *.gpl
./dg
gnuplot anim.gnu
```

If you modify `dg.cc`, then you must compile and run again.
