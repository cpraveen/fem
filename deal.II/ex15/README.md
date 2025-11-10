# Convection-diffusion equation

This is trying to solve a problem similar to 

> Quarteroni and Valli, Figure 8.3.1

But it does not give full details. Here is what we solve.

```text
-eps * Laplace(u) + du/dy = 0  in  (0,1) x (0,1)
u = 0   on  x = 0
u = 1   on  x = 1
u = 0.5*(1 + tanh(100*(x-0.5))) on y=0 and y=1
```

Compile and run

```shell
make
./demo
python plot.py
```

Solution does not look so bad, need to try with triangles.
