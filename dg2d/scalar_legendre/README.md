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

## Dirichlet bc: test_rotate.h

Copy `test_rotate.h` as `test_data.h`.

Set `periodic = false` in input file or in main function and run the code. The boundary condition has been set to zero. You should get very same solution as with periodic bc.

Now solve the rotating problem in one quarter domain with boundary conditions. Modify the main function

```c++
   param.xmin = 0.0; param.xmax = XMAX;
   param.ymin = 0.0; param.ymax = YMAX;
   parse_parameters(ph, param);

   param.final_time = 0.5 * M_PI;
   param.periodic = false;
   auto initial_condition = Solution<2>();
   auto boundary_condition = Solution<2>();
   auto exact_solution = Solution<2>();
```

We are using the exact solution of rotating gaussian as boundary condition which now depends on time. The gaussian will rotate by 90 degree at the final time of pi/2. You can also set the final time and periodic information in the `input.prm` file.

Also try with a time independent boundary condition

```c++
#include <deal.II/base/function_parser.h>

   auto boundary_condition = FunctionParser<2>("exp(-50*((x-0.5)^2 + y^2))");
```

For t > pi/2, we get stationary solution.
