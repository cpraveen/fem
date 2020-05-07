/* 
 Solve 2d laplace equation
 -Laplace(u) = f(x) in (0,1)x(0,1)
 Exact solution is u = sin(2*pi*x) * sin(2*pi*y)
 f is obtained from exact solution.
 Boundary condition is dirichlet and taken from exact solution.
*/
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
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
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

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
   return sin(2*M_PI*p[0])*sin(2*M_PI*p[1]);
}

template<>
Tensor<1,2> ExactSolution<2>::gradient (const Point<2>   &p,
                                        const unsigned int) const
{
   Tensor<1,2> values;
   values[0] = 2*M_PI*cos(2*M_PI*p[0])*sin(2*M_PI*p[1]);
   values[1] = 2*M_PI*sin(2*M_PI*p[0])*cos(2*M_PI*p[1]);
   return values;
}
//------------------------------------------------------------------------------
template <int dim>
class RightHandSide : public Function<dim>
{
public:
   RightHandSide () : Function<dim>() {}
   
   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
};

template <>
double RightHandSide<2>::value (const Point<2> &p,
                                  const unsigned int /*component*/) const
{
   return 8*M_PI*M_PI*sin(2*M_PI*p[0])*sin(2*M_PI*p[1]);
}
//------------------------------------------------------------------------------
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
   BoundaryValues () : Function<dim>() {}
   
   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
};

template <int dim>
double BoundaryValues<dim>::value (const Point<dim>& /*p*/,
                                   const unsigned int /*component*/) const
{
   return 0.0;
}

//------------------------------------------------------------------------------
template <int dim>
class LaplaceProblem
{
public:
   LaplaceProblem (int degree, unsigned int nrefine);
   void run (unsigned int &ncell, unsigned int &ndofs, 
             double &L2_error, double &H1_error);
   
private:
   void make_grid_and_dofs ();
   void assemble_system ();
   void solve ();
   void output_results () const;
   void compute_error (double &L2_error, double &H1_error) const;
   
   unsigned int           nrefine;
   Triangulation<dim>     triangulation;
   FE_Q<dim>              fe;
   DoFHandler<dim>        dof_handler;
   
   SparsityPattern        sparsity_pattern;
   SparseMatrix<double>   system_matrix;
   
   Vector<double>         solution;
   Vector<double>         system_rhs;
};


//------------------------------------------------------------------------------
template <int dim>
LaplaceProblem<dim>::LaplaceProblem (int degree, unsigned int nrefine)
:
nrefine (nrefine),
fe (degree),
dof_handler (triangulation)
{}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
   GridGenerator::hyper_cube (triangulation, 0, 1);
   triangulation.refine_global (nrefine);
   
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
   
   DynamicSparsityPattern dsp(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, dsp);
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
   system_rhs = 0;

   QGauss<dim>  quadrature_formula(2*fe.degree);
   const RightHandSide<dim> right_hand_side;
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values   | update_gradients |
                            update_quadrature_points | update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   
   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
   Vector<double>       cell_rhs (dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   std::vector<double>  rhs_values (n_q_points);
   
   for (const auto &cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;
      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
               cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                    fe_values.shape_grad (j, q_point) *
                                    fe_values.JxW (q_point));
            
            cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                            rhs_values[q_point] *
                            fe_values.JxW (q_point));
         }
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));
         
         system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }
   }
   
   // boundary condition
   std::map<types::global_dof_index,double> boundary_values;
   VectorTools::interpolate_boundary_values (dof_handler,
                                             0,
                                             BoundaryValues<dim>(),
                                             boundary_values);
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::solve ()
{
   SolverControl           solver_control (1000, 1e-12);
   SolverCG<>              cg (solver_control);
   cg.solve (system_matrix, solution, system_rhs,
             PreconditionIdentity());
   
   std::cout
   << "   " << solver_control.last_step()
   << " CG iterations needed to obtain convergence."
   << std::endl;
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::output_results () const
{
   Vector<double> solution_error;
   solution_error.reinit(dof_handler.n_dofs());
   VectorTools::interpolate(dof_handler,
                            ExactSolution<dim>(),
                            solution_error);
   solution_error -= solution;

   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution_error, "error");
   data_out.build_patches (fe.degree);
   std::string fname = "error-" + Utilities::int_to_string(nrefine,2)+".vtk";
   std::ofstream output (fname);
   data_out.write_vtk (output);
   std::cout << "   Wrote to file " << fname << std::endl;
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::compute_error (double &L2_error, double &H1_error) const
{
   // compute L2 error in solution
   ExactSolution<dim> exact_solution;
   Vector<double> difference_per_cell (triangulation.n_active_cells());
   VectorTools::integrate_difference (dof_handler,
                                      solution,
                                      exact_solution,
                                      difference_per_cell,
                                      QGauss<dim>(2*fe.degree+1),
                                      VectorTools::L2_norm);
   L2_error = difference_per_cell.l2_norm();

   // compute L2 error in gradient
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
void LaplaceProblem<dim>::run (unsigned int &ncell, 
                               unsigned int &ndofs, 
                               double       &L2_error, 
                               double       &H1_error)
{
   make_grid_and_dofs();
   assemble_system ();
   solve ();
   output_results ();
   compute_error (L2_error, H1_error);

   ncell = triangulation.n_active_cells ();
   ndofs = dof_handler.n_dofs ();
}

//------------------------------------------------------------------------------
int main ()
{
   deallog.depth_console (0);
   int degree = 1;
   ConvergenceTable  convergence_table;   
   for(unsigned int n=5; n<10; ++n)
   {
      LaplaceProblem<2> problem (degree, n);
      unsigned int ncell, ndofs;
      double L2_error, H1_error;
      problem.run (ncell, ndofs, L2_error, H1_error);

      convergence_table.add_value("cells", ncell);
      convergence_table.add_value("dofs",  ndofs);
      convergence_table.add_value("L2",    L2_error);
      convergence_table.add_value("H1",    H1_error);
      std::cout << "--------------------------------------------------------\n";
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

   convergence_table.evaluate_convergence_rates("L2",
                                                ConvergenceTable::reduction_rate_log2);
   convergence_table.evaluate_convergence_rates("H1",
                                                ConvergenceTable::reduction_rate_log2);

   std::cout << std::endl;
   convergence_table.write_text(std::cout);

   std::ofstream error_table_file("error.tex");   
   convergence_table.write_tex(error_table_file);
   
   return 0;
}
