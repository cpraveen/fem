# Neumann Problem

This folder contains the implementation of finding solution pure Neumann problem which is described in Chapter 4 of the book. The problem solved is:

$$
  -\Delta u = f \quad \text{in} \quad \Omega \qquad\textrm{and}\qquad
  n \cdot \nabla u = g \quad \text{on} \quad \partial\Omega
$$

We take $f,g$ such that

$$
\int_\Omega f + \int_{\partial\Omega} g = 0
$$

so that the problem has a solution. To get unique solution, the mean value of solution is set to zero.

$$
\int_\Omega u = 0
$$

Steps to run

```bash
cmake .
make release
make
./neumann
```

After running, the numerical error will be written to a LaTeX file named error.tex.

To compile it to PDF:
```bash
pdflatex error.tex
```
