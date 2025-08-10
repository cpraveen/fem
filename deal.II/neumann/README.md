# Neumann Problem

This folder contains the implementation of finding solution pure Neumann problem which is described in Chapter 4 of the book. The problem solved is:

$$
  -\Delta u &= f \quad \text{in } \Omega \\
  n \cdot \nabla u &= g_N \quad \text{on } \partial\Omega
$$

where:
- $\Omega$ is the computational domain (typically 2D),
- $f$ is a given source term.
- $g_N$ is a given function.

To get unique solution, the mean value of solution is set to zero.

Steps to run

```bash
cmake .
make
./neumann
```

After running, the numerical error will be written to a LaTeX file named error.tex.

To compile it to PDF:
```bash
pdflatex error.tex
```
