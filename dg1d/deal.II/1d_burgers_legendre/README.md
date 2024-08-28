# 1-d Burgers using DG

Copy file `dg.cc` from `1d_scalar_legendre`

```shell
cp ../1d_scalar_legendre/dg.cc .
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
