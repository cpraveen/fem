/* 
 Solve 1d laplace equation
 -u_xx = f(x)  for x in (0,1)
 Exact solution is u = x + sin(4 pi x)
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
class RightHandSide : public Function<1>
{
public:
   RightHandSide () : Function<1>() {}
   
   virtual double value (const Point<1>   &p,
                         const unsigned int  component = 0) const;
};

double RightHandSide::value (const Point<1> &p,
                                  const unsigned int /*component*/) const
{
   return 16*M_PI*M_PI*sin(4*M_PI*p[0]);
}

//------------------------------------------------------------------------------
class BoundaryValues : public Function<1>
{
public:
   BoundaryValues () : Function<1>() {}
   
   virtual double value (const Point<1>   &p,
                         const unsigned int  component = 0) const;
};

double BoundaryValues::value (const Point<1> &p,
                                   const unsigned int /*component*/) const
{
   return p[0] + sin(4*M_PI*p[0]);
}

//------------------------------------------------------------------------------
class LaplaceProblem
{
public:
   LaplaceProblem (int degree);
   void run ();
   
private:
   void make_grid_and_dofs ();
   void assemble_system ();
   void solve ();
   void output_results () const;
   
   Triangulation<1>     triangulation;
   FE_Q<1>              fe;
   DoFHandler<1>        dof_handler;
   
   SparsityPattern      sparsity_pattern;
   SparseMatrix<double> system_matrix;
   
   Vector<double>       solution;
   Vector<double>       system_rhs;
};


//------------------------------------------------------------------------------
LaplaceProblem::LaplaceProblem (int degree) :
fe (degree),
dof_handler (triangulation)
{}

//------------------------------------------------------------------------------
void LaplaceProblem::make_grid_and_dofs ()
{
   GridGenerator::hyper_cube (triangulation, 0, 1);
   triangulation.refine_global (8);
   
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
void LaplaceProblem::assemble_system ()
{
   QGauss<1>  quadrature_formula(2*fe.degree);
   
   const RightHandSide right_hand_side;
   
   FEValues<1> fe_values (fe, quadrature_formula,
                            update_values   | update_gradients |
                            update_quadrature_points | update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   
   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
   Vector<double>       cell_rhs (dofs_per_cell);
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   std::vector<double> rhs_values (n_q_points);
   
   typename DoFHandler<1>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
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
   
   // left boundary condition
   std::map<unsigned int,double> boundary_values;
   VectorTools::interpolate_boundary_values (dof_handler,
                                             0,
                                             BoundaryValues(),
                                             boundary_values);
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
   // right boundary condition
   VectorTools::interpolate_boundary_values (dof_handler,
                                             1,
                                             BoundaryValues(),
                                             boundary_values);
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}

//------------------------------------------------------------------------------
void LaplaceProblem::solve ()
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
void LaplaceProblem::output_results () const
{
   DataOut<1> data_out;
   
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches (fe.degree);
   std::ofstream output ("solution.gnuplot");
   data_out.write_gnuplot (output);
}

//------------------------------------------------------------------------------
void LaplaceProblem::run ()
{
   make_grid_and_dofs();
   assemble_system ();
   solve ();
   output_results ();
}

//------------------------------------------------------------------------------
int main ()
{
   deallog.depth_console (0);
   int degree = 2;
   LaplaceProblem laplace_problem_1d (degree);
   laplace_problem_1d.run ();
   
   return 0;
}
