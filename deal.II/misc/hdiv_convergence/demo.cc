// Convergence of RT interpolation on Cartesian and quadrilateral grids
// Initial code generated with Google Gemini and then modified
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_raviart_thomas.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

// Define exact vector field u(x, y) = [ sin(pi * x), sin(pi * y) ]^T
template <int dim>
class ExactVectorField : public Function<dim>
{
public:
  ExactVectorField() : Function<dim>(dim) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const override
  {
    return std::sin(numbers::PI * p[component]);
  }

  virtual void vector_value(const Point<dim> &p,
                            Vector<double>   &values) const override
  {
    for (unsigned int d = 0; d < dim; ++d)
      values[d] = std::sin(numbers::PI * p[d]);
  }

  virtual Tensor<1, dim> gradient(const Point<dim> &p,
                                  const unsigned int component = 0) const override
  {
    Tensor<1, dim> grad;
    grad[component] = numbers::PI * std::cos(numbers::PI * p[component]);
    return grad;
  }

  virtual void vector_gradient(const Point<dim> &p,
                               std::vector<Tensor<1, dim>> &gradients) const override
  {
    for (unsigned int d = 0; d < dim; ++d)
      {
        gradients[d] = 0;
        gradients[d][d] = numbers::PI * std::cos(numbers::PI * p[d]);
      }
  }
};

template <int dim>
class RTInterpolationTest
{
public:
  RTInterpolationTest(const unsigned int degree);
  void run(const unsigned int initial_refine,
           const unsigned int n_cycles,
           const double perturb);

private:
  void perturb_grid(const double perturb, const unsigned refine);
  void setup_and_interpolate();
  void compute_and_record_errors(const unsigned int cycle);

  const unsigned int degree;

  Triangulation<dim> triangulation;
  FE_RaviartThomas<dim> fe;
  DoFHandler<dim> dof_handler;
  Vector<double> interpolated_solution;

  ConvergenceTable convergence_table;
  ExactVectorField<dim> exact_function;
};

template <int dim>
RTInterpolationTest<dim>::RTInterpolationTest(const unsigned int degree)
  : degree(degree)
  , fe(degree)
  , dof_handler(triangulation)
{}

template <int dim>
void RTInterpolationTest<dim>::setup_and_interpolate()
{
  dof_handler.distribute_dofs(fe);
  interpolated_solution.reinit(dof_handler.n_dofs());

  // Interpolate the exact vector field onto the RT space
  VectorTools::interpolate(dof_handler, exact_function, interpolated_solution);
}

template <int dim>
void RTInterpolationTest<dim>::compute_and_record_errors(const unsigned int cycle)
{
  const QGauss<dim> quadrature(degree + 3);

  // 1. Compute L2 Error
  Vector<double> cell_l2_errors(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    interpolated_solution,
                                    exact_function,
                                    cell_l2_errors,
                                    quadrature,
                                    VectorTools::L2_norm);
  const double l2_error = VectorTools::compute_global_error(triangulation,
                                                            cell_l2_errors,
                                                            VectorTools::L2_norm);

  // 2. Compute H(div) Semi-Norm Error (Divergence Error)
  Vector<double> cell_hdiv_semi_errors(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    interpolated_solution,
                                    exact_function,
                                    cell_hdiv_semi_errors,
                                    quadrature,
                                    VectorTools::Hdiv_seminorm);
  const double hdiv_semi_error = VectorTools::compute_global_error(triangulation,
                                                                   cell_hdiv_semi_errors,
                                                                   VectorTools::Hdiv_seminorm);

  // Record metrics in ConvergenceTable
  convergence_table.add_value("cycle", cycle);
  convergence_table.add_value("cells", triangulation.n_active_cells());
  convergence_table.add_value("dofs", dof_handler.n_dofs());
  convergence_table.add_value("L2_error", l2_error);
  convergence_table.add_value("Hdiv_error", hdiv_semi_error);
}

template <int dim>
void RTInterpolationTest<dim>::perturb_grid(const double perturb,
                                            const unsigned int refine)
{
  if(perturb <= 0.0) return;

  // Randomly perturb the grid
  const auto h = 1.0 / pow(2, refine);
  auto v = triangulation.begin_vertex();
  for (; v < triangulation.end_vertex(); ++v)
  {
    auto &p = v->vertex();
    if (p[0] > 0 && p[0] < 1 && p[1] > 0 && p[1] < 1)
    {
      double a = ((double)rand()) / RAND_MAX;
      double b = ((double)rand()) / RAND_MAX;
      p[0] += perturb * (2 * a - 1) * h;
      p[1] += perturb * (2 * b - 1) * h;
    }
  }
}

template <int dim>
void RTInterpolationTest<dim>::run(const unsigned int initial_refine,
                                   const unsigned int n_cycles,
                                   const double perturb)
{
  for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
    {
      GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
      const auto refine = initial_refine + cycle;
      triangulation.refine_global(refine);
      perturb_grid(perturb, refine);
      setup_and_interpolate();
      compute_and_record_errors(cycle);
    }

  // Formatting output table
  convergence_table.set_precision("L2_error", 3);
  convergence_table.set_scientific("L2_error", true);
  convergence_table.set_precision("Hdiv_error", 3);
  convergence_table.set_scientific("Hdiv_error", true);

  // Evaluate rates relative to mesh refinement
  convergence_table.evaluate_convergence_rates("L2_error", ConvergenceTable::reduction_rate_log2);
  convergence_table.evaluate_convergence_rates("Hdiv_error", ConvergenceTable::reduction_rate_log2);

  std::cout << "\n==============================================================\n";
  std::cout << " Convergence Results for RT(" << degree << ") Space in " << dim << "D ";
  if(perturb > 0.0)
     cout << " Quadrilateral grid\n";
  else
     cout << " Cartesian grid\n";
  std::cout << "==============================================================\n\n";
  convergence_table.write_text(std::cout);
}

int main()
{
  try
    {
      const unsigned int dim = 2;
      const unsigned int initial_refine = 3;
      const unsigned int n_cycles = 5;
      const double perturb = 0.1;

      for(unsigned int degree = 0; degree <= 1; ++degree)
      {
         RTInterpolationTest<dim> test(degree);
         test.run(initial_refine, n_cycles, 0.0);
      }

      for(unsigned int degree = 0; degree <= 1; ++degree)
      {
         RTInterpolationTest<dim> test(degree);
         test.run(initial_refine, n_cycles, perturb);
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << "\n\n----------------------------------------------------\n"
                << "Exception on processing: \n"
                << exc.what() << "\nAborting!\n"
                << "----------------------------------------------------\n";
      return 1;
    }
  catch (...)
    {
      std::cerr << "\n\n----------------------------------------------------\n"
                << "Unknown exception!\n"
                << "----------------------------------------------------\n";
      return 1;
    }

  return 0;
}
