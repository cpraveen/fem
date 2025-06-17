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
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/base/convergence_table.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;


class ExactSolution : public Function<2>
{
public:
   ExactSolution() : Function<2>() {}

   double value(const Point<2>&   p,
                const unsigned int  component = 0) const override;
   Tensor<1, 2> gradient(const Point<2>&   p,
                         const unsigned int  component = 0) const override;
};

double
ExactSolution::value(const Point<2>& p, const unsigned int /*component*/) const
{
   return sin(p[0] * p[1]);
}

Tensor<1, 2>
ExactSolution::gradient(const Point<2>&   p, const unsigned int) const
{
   Tensor<1, 2> values;
   values[0] = p[1] * cos(p[0] * p[1]);
   values[1] = p[0] * cos(p[0] * p[1]);
   return values;
}

class Projection
{
public:
   Projection(unsigned int nrefine);
   void run(std::vector<int>& ncell,
            std::vector<int>& ndofs,
            std::vector<double>& L2_error,
            std::vector<double>& H1_error);

private:
   void make_grid_and_dofs();
   void setup_system();
   void assemble_system();
   void solve();
   void compute_error(double& L2_error, double& H1_error);
   void refine_grid();

   unsigned int            nrefine;
   Triangulation<2>        triangulation;
   const FE_SimplexP<2>    fe;
   MappingFE<2>            mapping;
   DoFHandler<2>           dof_handler;

   SparsityPattern         sparsity_pattern;
   SparseMatrix<double>    mass_matrix;

   Vector<double>          system_rhs;
   Vector<double>          solution;
};

Projection::Projection(unsigned int nrefine):
   nrefine(nrefine),
   fe(1),
   mapping(FE_SimplexP<2>(1)),
   dof_handler(triangulation)
{}


void
Projection::make_grid_and_dofs()
{
   Triangulation<2> tria;
   GridGenerator::hyper_cube(tria);
   GridGenerator::convert_hypercube_to_simplex_mesh(tria, triangulation);
   triangulation.refine_global(2);
}

void
Projection::setup_system()
{
   dof_handler.distribute_dofs(fe);

   DynamicSparsityPattern dsp(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern(dof_handler, dsp);
   sparsity_pattern.copy_from(dsp);

//  Initializing the all the matrices
   mass_matrix.reinit(sparsity_pattern);
   system_rhs.reinit(dof_handler.n_dofs());
   solution.reinit(dof_handler.n_dofs());
}

void
Projection::assemble_system()
{
   mass_matrix = 0;
   system_rhs = 0;

   QGaussSimplex<2> quadrature_formula(fe.degree + 1);
   FEValues<2> fe_values(mapping, fe, quadrature_formula,
                         update_values |
                         update_quadrature_points |
                         update_JxW_values);

   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   ExactSolution exact_solution;

   FullMatrix<double>   cell_matrix(dofs_per_cell, dofs_per_cell);
   Vector<double>       cell_rhs(dofs_per_cell);
   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   for(const auto  &cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs = 0;

      for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
         const auto value = exact_solution.value(fe_values.quadrature_point(q_point));
         for(unsigned int i = 0; i < dofs_per_cell ; ++i)
         {
            for(unsigned int j = 0; j < dofs_per_cell; ++j)
            {
               cell_matrix(i, j) +=   fe_values.shape_value(i, q_point)
                                      * fe_values.shape_value(j, q_point)
                                      * fe_values.JxW(q_point);
            }
            cell_rhs(i) +=   value
                             * fe_values.shape_value(i, q_point)
                             * fe_values.JxW(q_point);
         }
      }

      cell->get_dof_indices(local_dof_indices);
      for(unsigned int i = 0; i < dofs_per_cell ; ++i)
      {
         for(unsigned int j = 0; j < dofs_per_cell; ++j)
         {
            mass_matrix.add(local_dof_indices[i], 
                            local_dof_indices[j], 
                            cell_matrix(i, j));
         }
         system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }

   }
}

void
Projection::solve()
{
   SolverControl            solver_control(1000, 1e-12 * system_rhs.l2_norm());
   SolverCG<Vector<double>> cg(solver_control);
   cg.solve(mass_matrix, solution, system_rhs, PreconditionIdentity());
}

void Projection::compute_error(double& L2_error, double& H1_error)
{
   ExactSolution exact_solution;

   Vector<double> difference_per_cell(triangulation.n_active_cells());
   VectorTools::integrate_difference(mapping,
                                     dof_handler,
                                     solution,
                                     exact_solution,
                                     difference_per_cell,
                                     QGaussSimplex<2>(fe.degree + 2),
                                     VectorTools::L2_norm);

   L2_error = difference_per_cell.l2_norm();

   VectorTools::integrate_difference(mapping,
                                     dof_handler,
                                     solution,
                                     exact_solution,
                                     difference_per_cell,
                                     QGaussSimplex<2>(fe.degree + 2),
                                     VectorTools::H1_seminorm);

   H1_error = difference_per_cell.l2_norm();
}

void
Projection::refine_grid()
{
   triangulation.refine_global(1);
}

void
Projection::run(std::vector<int>& ncell,
                std::vector<int>& ndofs,
                std::vector<double>& L2_error,
                std::vector<double>& H1_error)
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
      solve();
      compute_error(L2_error[n], H1_error[n]);

      ncell[n] = triangulation.n_active_cells();
      ndofs[n] = dof_handler.n_dofs();
   }
}

int
main()
{
   deallog.depth_console(0);
   unsigned int nrefine = 8;
   Projection problem(nrefine);
   std::vector<int> ncell(nrefine), ndofs(nrefine);
   std::vector<double> L2_error(nrefine), H1_error(nrefine);
   problem.run(ncell, ndofs, L2_error, H1_error);
   ConvergenceTable  convergence_table;
   for(unsigned int n = 0; n < nrefine; ++n)
   {
      convergence_table.add_value("cells", ncell[n]);
      convergence_table.add_value("dofs",  ndofs[n]);
      convergence_table.add_value("L2",    L2_error[n]);
      convergence_table.add_value("H1",    H1_error[n]);
   }

   convergence_table.set_precision("L2", 3);
   convergence_table.set_scientific("L2", true);

   convergence_table.set_precision("H1", 3);
   convergence_table.set_scientific("H1", true);

   convergence_table.set_tex_caption("cells", "\\# cells");
   convergence_table.set_tex_caption("dofs", "\\# dofs");
   convergence_table.set_tex_caption("L2", "$L^2$-error");
   convergence_table.set_tex_caption("H1", "$H^1$-error");

   convergence_table.set_tex_format("cells", "r");
   convergence_table.set_tex_format("dofs",  "r");

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
