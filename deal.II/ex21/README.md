# Unsteady heat equation

This is based on step-26 from deal.II tutorial. I have removed grid adaptation to make this simplified test case. See [step-26](https://dealii.org/current/doxygen/deal.II/step_26.html) for explanation of the method.

We solve heat equation with theta scheme. `theta=0.5` is the Crank-Nicholson scheme, whose value is set in the constructor.
