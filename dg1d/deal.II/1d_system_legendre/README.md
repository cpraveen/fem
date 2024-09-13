# 1-d system PDE using DG

This uses Legendre polynomials as basis functions.

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
