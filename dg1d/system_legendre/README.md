# 1-d system PDE using DG

This uses Legendre polynomials as basis functions.

* `acoustics`: Linear acoustic in homogeneous medium
* `shallow`: Shallow water with flat bottom
* `euler`: compressible Euler equation for idea gas

## acoustics: solves linear acoustics

```shell
cp acoustics/* .
cmake .
make release
make
rm -f *.gpl
./main input.prm
gnuplot anim.gnu
```

## Extension

The code in file `dg.h` works for any PDE system and should not need any modification, unless you want to implement, say, your own limiter scheme.

To implement your own PDE and problem, you need to write two files

* `pde.h`: contains PDE specific things
* `problem_data.h`: contains problem specific things

See `acoustics` directory for an example.
