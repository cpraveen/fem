# 1-D DG code for Euler equations

This code solves 1-D Euler equations using a DG scheme.

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
