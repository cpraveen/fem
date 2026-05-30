# 1-D DG code for Euler equations

This code solves 1-D Euler equations using a DG scheme and requires Python3, Numpy and Matplotlib.

* Basis functions: Legendre polynomials
* Numerical flux: Rusanov scheme
* Limiter: TVD and TVB (conserved and characteristic)

By default, the code is set to run the Sod test case.

Run with conserved limiter

```
python euler.py 
```

Run with characteristic limiter

```
python euler.py -char_lim 1
```

For other test cases, pass `-ic` argument, see the help

```
python euler.py -h
```

for available options.

## Sod test case

```shell
python euler.py -ncell 100 -degree 2 -Tf 0.2 -ic sod -char_lim 1 -tvbM 100
```

## Shu-Osher test

```shell
python euler.py -ncell 400 -degree 2 -Tf 1.8 -ic shuosher -char_lim 1 \
                -tvbM 100 -plot_freq 100
```
