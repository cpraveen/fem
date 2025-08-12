/*
Solves pure Neumann problem
   -Laplace(u) = f in Omega
         du/dn = g on dOmega
and we have compatibility
   int(Omega) f + int(dOmega) g = 0
To get unique solution, we put constraint
   int(Omega) u = 0
Author : Devansh Sonigra
Modifed: Praveen C
*/
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_mic.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

template<int dim>
class ExactSolution : public Function<dim>
{
public:
  ExactSolution() : Function<dim>() {}

  double value(const Point<dim>&   p,
               const unsigned int  component = 0) const override;
  Tensor<1, dim> gradient(const Point<dim>&   p,
                          const unsigned int  component = 0) const override;
};

template<>
double
ExactSolution<2>::value(const Point<2>& p,
                        const unsigned int /*component*/) const
{
  return p[0] * p[0] + p[1] * p[1] - 2.0 / 3.0;
}


template<>
Tensor<1, 2>
ExactSolution<2>::gradient(const Point<2>&   p,
                           const unsigned int) const
{
  Tensor<1, 2> values;
  values[0] = 2 * p[0];
  values[1] = 2 * p[1];
  return values;
}

template<int dim>
class RHS_function : public Function<dim>
{
public:
  RHS_function() : Function<dim>() {}

  double value(const Point<dim>&   p,
               const unsigned int  component = 0) const override;
};

template<>
double
RHS_function<2>::value(const Point<2>& /*p*/,
                       const unsigned int /*component*/) const
{
  return -4;
}

template<int dim>
class G_Neumann: public Function<dim>
{
public:
  G_Neumann() : Function<dim>() {}
  double value(const Point<dim>& p,
               const unsigned int component = 0) const override;
};

template<>
double
G_Neumann<2>::value(const Point<2>& p, const unsigned int) const
{
  const auto x = p[0];
  const auto y = p[1];
  const auto EPS = 1.0e-13;

  if(fabs(x) < EPS) // left
  {
    return 0.0;
  }
  else if(fabs(x-1.0) < EPS) // right
  {
    return 2.0;
  }
  else if(fabs(y) < EPS) // bottom
  {
    return 0.0;
  }
  else if(fabs(y-1.0) < EPS) // top
  {
    return 2.0;
  }
  else
  {
    DEAL_II_NOT_IMPLEMENTED();
  }
}

template<int dim>
class NeumannSolver
{
public:
  NeumannSolver(unsigned int nrefine, unsigned int degree);
  void run(std::vector<int>& ncell,
           std::vector<int>& ndofs,
           std::vector<double>& L2_error,
           std::vector<double>& H1_error,
           std::vector<int>& niterations);

private:
  void make_grid_and_dofs();
  void setup_system();
  void assemble_system();
  void solve(int& niteration);
  void compute_error(double& L2_error, double& H1_error);
  void refine_grid();

  unsigned int            nrefine;
  Triangulation<dim>      triangulation;
  const FE_SimplexP<dim>  fe;
  MappingFE<dim>          mapping;
  DoFHandler<dim>         dof_handler;

  SparsityPattern         sparsity_pattern;
  SparseMatrix<double>    system_matrix;

  Vector<double>          system_rhs;
  Vector<double>         system_solution;
};

template<int dim>
NeumannSolver<dim>::NeumannSolver(unsigned int nrefine, unsigned int degree):
  nrefine(nrefine),
  fe(degree),
  mapping(FE_SimplexP<dim>(1)),
  dof_handler(triangulation)
{}

template<int dim>
void
NeumannSolver<dim>::make_grid_and_dofs()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  GridGenerator::convert_hypercube_to_simplex_mesh(tria, triangulation);
  triangulation.refine_global(2);
}

template<int dim>
void
NeumannSolver<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  // Sparsity pattern WITHOUT Lagrange multiplier
  DynamicSparsityPattern dsp1(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp1);

  // Sparsity pattern WITH Lagrange multiplier
  DynamicSparsityPattern dsp2(dof_handler.n_dofs() + 1);

  // First copy dsp1 into dsp2
  for(types::global_dof_index i=0; i<dsp1.n_rows(); ++i)
    for(auto j = dsp1.begin(i); j < dsp1.end(i); ++j)
      dsp2.add(i, j->column());

  // Add last column/row
  for(unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
  {
    dsp2.add(i, dof_handler.n_dofs());
    dsp2.add(dof_handler.n_dofs(), i);
  }

  sparsity_pattern.copy_from(dsp2);

  //  Initializing the all the matrices
  system_matrix.reinit(sparsity_pattern);
  system_rhs.reinit(dof_handler.n_dofs() + 1);
  system_solution.reinit(dof_handler.n_dofs() + 1);
}


template<int dim>
void
NeumannSolver<dim>::assemble_system()
{
  // This takes input of what polynomial degree to be integrated exactly
  QGaussSimplex<dim> quadrature_formula(fe.degree + 1);
  FEValues<dim> fe_values(mapping, fe, quadrature_formula,
                          update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);

  QGaussSimplex < dim - 1 > face_quadrature_formula(fe.degree + 1);
  FEFaceValues<dim> face_fe_values(mapping, fe, face_quadrature_formula, 
                                   update_values | update_quadrature_points |
                                   update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  const unsigned int   face_n_q_points = face_quadrature_formula.size();

  RHS_function<dim> rhs_function;
  G_Neumann<dim> g_neu;

  FullMatrix<double>   cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs(dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  system_matrix = 0;
  system_rhs = 0;

  for(const auto&  cell : dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    std::vector<double>  basis_integrals(dofs_per_cell, 0.0);
    cell_matrix = 0;
    cell_rhs = 0;

    for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
    {
      double temp = rhs_function.value(fe_values.quadrature_point(q_point));
      for(unsigned int i = 0; i < dofs_per_cell ; ++i)
      {
        for(unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
                               fe_values.shape_grad(j, q_point) *
                               fe_values.JxW(q_point);
        }
        cell_rhs(i) += temp * 
                       fe_values.shape_value(i, q_point) *
                       fe_values.JxW(q_point);
        basis_integrals[i] += fe_values.shape_value(i, q_point) *
                              fe_values.JxW(q_point);
      }
    }

    for(unsigned int f = 0; f < cell->n_faces(); ++f)
    {
      if(cell->face(f)->at_boundary())
      {
        face_fe_values.reinit(cell, f);
        for(unsigned int q_point = 0; q_point < face_n_q_points; ++q_point)
        {
          double g_neu_value = g_neu.value(face_fe_values.quadrature_point(q_point));
          for(unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            cell_rhs(i) += g_neu_value *
                           face_fe_values.shape_value(i, q_point) *
                           face_fe_values.JxW(q_point);
          }
        }
      }
    }

    cell->get_dof_indices(local_dof_indices);
    for(unsigned int i = 0; i < dofs_per_cell ; ++i)
    {
      for(unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        system_matrix.add(local_dof_indices[i],
                          local_dof_indices[j],
                          cell_matrix(i, j));
      }
      system_rhs(local_dof_indices[i]) += cell_rhs(i);
      system_matrix.add(local_dof_indices[i],
                        dof_handler.n_dofs(),
                        basis_integrals[i]);
      system_matrix.add(dof_handler.n_dofs(),
                        local_dof_indices[i],
                        basis_integrals[i]);
    }

  }
}

template<int dim>
void
NeumannSolver<dim>::solve(int& niteration)
{
  SolverControl           solver_control(10000, 1e-12 * system_rhs.l2_norm());
  SolverGMRES<Vector<double>> cg(solver_control);

  SparseILU<double> preconditioner;
  preconditioner.initialize(system_matrix);

  // cg.solve takes solution as initial guess
  cg.solve(system_matrix, system_solution, system_rhs, preconditioner);
  niteration = solver_control.last_step();
}

template<int dim>
void
NeumannSolver<dim>::compute_error(double& L2_error, double& H1_error)
{
  ExactSolution<dim> exact_solution;
  Vector<double> solution(dof_handler.n_dofs());

  // Copy into solution excluding Lagrange multiplier, which is last component
  for(unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
  {
    solution(i) = system_solution(i);
  }

  Vector<double> difference_per_cell(triangulation.n_active_cells());
  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    difference_per_cell,
                                    QGaussSimplex<dim>(2 * fe.degree + 1),
                                    VectorTools::L2_norm);

  L2_error = difference_per_cell.l2_norm();

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    difference_per_cell,
                                    QGaussSimplex<dim>(2 * fe.degree + 1),
                                    VectorTools::H1_seminorm);

  H1_error = difference_per_cell.l2_norm();
}

template<int dim>
void
NeumannSolver<dim>::refine_grid()
{
  triangulation.refine_global(1);
}

template<int dim>
void
NeumannSolver<dim>::run(std::vector<int>& ncell,
                        std::vector<int>& ndofs,
                        std::vector<double>& L2_error,
                        std::vector<double>& H1_error,
                        std::vector<int>& niterations)
{
  for(unsigned int n = 0; n < nrefine; ++n)
  {
    if(n == 0)
    {
      make_grid_and_dofs();
    }
    else
    {
      refine_grid();
    }

    setup_system();
    assemble_system();
    solve(niterations[n]);
    compute_error(L2_error[n], H1_error[n]);

    ncell[n] = triangulation.n_active_cells();
    ndofs[n] = dof_handler.n_dofs();
  }
}

int
main()
{
  deallog.depth_console(0);
  unsigned int nrefine = 6;
  unsigned int degree = 1;

  NeumannSolver<2> problem(nrefine, degree);
  std::vector<int> ncell(nrefine), ndofs(nrefine), niterations(nrefine);
  std::vector<double> L2_error(nrefine), H1_error(nrefine);
  problem.run(ncell, ndofs, L2_error, H1_error, niterations);

  ConvergenceTable  convergence_table;
  for(unsigned int n = 0; n < nrefine; ++n)
  {
    convergence_table.add_value("cells", ncell[n]);
    convergence_table.add_value("dofs",  ndofs[n]);
    convergence_table.add_value("iterations",  niterations[n]);
    convergence_table.add_value("L2",    L2_error[n]);
    convergence_table.add_value("H1",    H1_error[n]);
  }

  convergence_table.set_precision("L2", 3);
  convergence_table.set_scientific("L2", true);

  convergence_table.set_precision("H1", 3);
  convergence_table.set_scientific("H1", true);

  convergence_table.set_tex_caption("cells", "\\# cells");
  convergence_table.set_tex_caption("dofs", "\\# dofs");
  convergence_table.set_tex_caption("iterations", "\\# iterations");
  convergence_table.set_tex_caption("L2", "$L^2$-error");
  convergence_table.set_tex_caption("H1", "$H^1$-error");

  convergence_table.set_tex_format("cells", "r");
  convergence_table.set_tex_format("dofs",  "r");
  convergence_table.set_tex_format("iterations",  "r");

  convergence_table.evaluate_convergence_rates
  ("L2", ConvergenceTable::reduction_rate_log2);
  convergence_table.evaluate_convergence_rates
  ("H1", ConvergenceTable::reduction_rate_log2);

  std::cout << std::endl;
  convergence_table.write_text(std::cout);

  std::ofstream error_table_file("error.tex");
  convergence_table.write_tex(error_table_file);

  return 0;
}
