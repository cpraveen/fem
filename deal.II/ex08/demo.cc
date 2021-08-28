/*
 Solve 2d laplace equation
 -Laplace(u) = 0  in Gamma shaped domain
 Exact solution is u = r^(2/3) * sin(2*theta/3)
 Boundary condition is dirichlet and taken from exact solution.
*/
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <fstream>
#include <iostream>


using namespace dealii;

//------------------------------------------------------------------------------
template <int dim>
class ExactSolution : public Function<dim>
{
public:
   ExactSolution () : Function<dim>() {}

   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
   virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                   const unsigned int  component = 0) const;
};

template <>
double ExactSolution<2>::value (const Point<2> &p,
                                const unsigned int /*component*/) const
{
   return (1.0 / M_PI) * atan2(p[1], p[0]);
}

template<>
Tensor<1,2> ExactSolution<2>::gradient (const Point<2>   &p,
                                        const unsigned int) const
{
   double r2 = p.norm_square();
   Tensor<1,2> values;
   values[0] = -p[1] / (M_PI * r2);
   values[1] =  p[0] / (M_PI * r2);
   return values;
}

//------------------------------------------------------------------------------
template <int dim>
class LaplaceProblem
{
public:
   LaplaceProblem (int degree, unsigned int nrefine);
   void run (std::vector<unsigned int> &ncell,
             std::vector<unsigned int> &ndofs,
             std::vector<double>       &L2_error,
             std::vector<double>       &H1_error);

private:
   void setup_system ();
   void assemble_system ();
   void solve ();
   void output_results ();
   void compute_error (double &L2_error, double &H1_error) const;
   void refine_grid ();

   unsigned int              nrefine;
   Triangulation<dim>        triangulation;
   FE_Q<dim>                 fe;
   DoFHandler<dim>           dof_handler;
   AffineConstraints<double> constraints;

   SparsityPattern           sparsity_pattern;
   SparseMatrix<double>      system_matrix;

   Vector<double>            solution;
   Vector<double>            system_rhs;
};


//------------------------------------------------------------------------------
template <int dim>
LaplaceProblem<dim>::LaplaceProblem (int degree, unsigned int nrefine) :
nrefine (nrefine),
fe (degree),
dof_handler (triangulation)
{}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
   std::cout
   << "   Number of active cells: "
   << triangulation.n_active_cells()
   << std::endl
   << "   Total number of cells: "
   << triangulation.n_cells()
   << std::endl;

   dof_handler.distribute_dofs (fe);

   std::cout
   << "   Number of degrees of freedom: "
   << dof_handler.n_dofs()
   << std::endl;

   constraints.clear ();
   DoFTools::make_hanging_node_constraints (dof_handler,
                                            constraints);
   VectorTools::interpolate_boundary_values (dof_handler,
                                             0,
                                             ExactSolution<dim>(),
                                             constraints);
   constraints.close ();

   DynamicSparsityPattern dsp(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints);
   sparsity_pattern.copy_from(dsp);

   system_matrix.reinit (sparsity_pattern);
   solution.reinit (dof_handler.n_dofs());
   system_rhs.reinit (dof_handler.n_dofs());
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
   system_matrix = 0;
   system_rhs    = 0;

   QGauss<dim>  quadrature_formula(2*fe.degree);
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values   | update_gradients |
                            update_quadrature_points | update_JxW_values);

   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();

   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
   Vector<double>       cell_rhs (dofs_per_cell);
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
               cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                    fe_values.shape_grad (j, q_point) *
                                    fe_values.JxW (q_point));
         }

      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (cell_matrix,
                                              cell_rhs,
                                              local_dof_indices,
                                              system_matrix,
                                              system_rhs);
   }

}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::solve ()
{
   SolverControl           solver_control (1000, 1e-12);
   SolverCG<>              cg (solver_control);

   PreconditionSSOR<SparseMatrix<double>> preconditioner;
   preconditioner.initialize(system_matrix, 1.2);

   cg.solve(system_matrix, solution, system_rhs,
            preconditioner);
   constraints.distribute (solution);

   std::cout
   << "   " << solver_control.last_step()
   << " CG iterations needed to obtain convergence."
   << std::endl;
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::output_results ()
{
   static int counter = 0;
   // compute error into system_rhs
   VectorTools::interpolate(dof_handler, ExactSolution<dim>(), system_rhs);
   constraints.distribute (system_rhs);
   system_rhs -= solution;

   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.add_data_vector (system_rhs, "error");
   data_out.build_patches (fe.degree);
   std::string fname = "solution-" + Utilities::int_to_string(counter,2)+".vtk";
   std::ofstream output (fname);
   data_out.write_vtk (output);
   ++counter;
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::compute_error (double &L2_error, double &H1_error) const
{
   std::cout << "Computing error norms\n";
   // compute error in solution
   ExactSolution<dim> exact_solution;
   Vector<double> difference_per_cell (triangulation.n_active_cells());
   VectorTools::integrate_difference (dof_handler,
                                      solution,
                                      exact_solution,
                                      difference_per_cell,
                                      QGauss<dim>(2*fe.degree+1),
                                      VectorTools::L2_norm);
   L2_error = difference_per_cell.l2_norm();

   // compute error in gradient
   VectorTools::integrate_difference (dof_handler,
                                      solution,
                                      exact_solution,
                                      difference_per_cell,
                                      QGauss<dim>(2*fe.degree+1),
                                      VectorTools::H1_seminorm);
   H1_error = difference_per_cell.l2_norm();
}
//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
   std::cout << "Refining grid based on error estimator\n";
   Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
   KellyErrorEstimator<dim>::estimate (dof_handler,
                                       QGauss<dim-1>(fe.degree+2),
                                       std::map<types::boundary_id, const Function<dim> *>(),
                                       solution,
                                       estimated_error_per_cell);
   GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
                                                      estimated_error_per_cell,
                                                      0.3, 0.03);
   triangulation.execute_coarsening_and_refinement ();
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::run (std::vector<unsigned int> &ncell,
                               std::vector<unsigned int> &ndofs,
                               std::vector<double>       &L2_error,
                               std::vector<double>       &H1_error)
{
   for(unsigned int n=0; n<nrefine; ++n)
   {
      if(n==0)
      {
         Point<dim> p1(-0.5, 0.0);
         Point<dim> p2(+0.5, 1.0);
         GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                   {20, 20},
                                                   p1,
                                                   p2);
      }
      else
         refine_grid ();

      setup_system();
      assemble_system ();
      solve ();
      output_results ();
      compute_error (L2_error[n], H1_error[n]);

      ncell[n] = triangulation.n_active_cells ();
      ndofs[n] = dof_handler.n_dofs ();
   }
}

//------------------------------------------------------------------------------
int main ()
{
   deallog.depth_console (0);
   int degree = 1;
   unsigned int nrefine = 7;
   LaplaceProblem<2> problem (degree, nrefine);
   std::vector<unsigned int> ncell(nrefine), ndofs(nrefine);
   std::vector<double> L2_error(nrefine), H1_error(nrefine);
   problem.run (ncell, ndofs, L2_error, H1_error);
   ConvergenceTable  convergence_table;
   for(unsigned int n=0; n<nrefine; ++n)
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

   std::cout << std::endl;
   convergence_table.write_text(std::cout);

   std::ofstream error_table_file("error.tex");
   convergence_table.write_tex(error_table_file);

   return 0;
}
