# Notes on implementation

- In some parts of the code we forego vectorization of `VectorizedArray` and loop over its elements to do certain calculations, I have listed the location which have a potential to be optimized. 
    - `compute_dt`
    - `initialize`
    - `assemble_mass_matrix`
- 
