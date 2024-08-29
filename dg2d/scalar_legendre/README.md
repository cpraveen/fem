# 2d scalar conservation law

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
rm -f *.vtu
./dg input.prm
visit -o sol*.vtu
```
