# 2d scalar conservation law

There are two test cases available.

* `test_linear.h`: advection along straight line, period = 2
* `test_rotate.h`: advection along circles, period = 2*pi

Copy one of these files as `test_data.h`.

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
