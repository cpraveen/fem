/*
 Solve 2d/3d laplace equation
 -Laplace(u) = f(x,y) in (0,1)x(0,1)
 Exact solution is u = sin(2*pi*x) * sin(2*pi*y)
 f is obtained from exact solution.
 Boundary condition is dirichlet and taken from exact solution.
*/
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
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
// RHS function f in Poisson equation
template <int dim>
class RightHandSide : public Function<dim>
{
public:
   RightHandSide () : Function<dim>() {}

   double value (const Point<dim>   &p,
                 const unsigned int  component = 0) const override;
};

// RHS in 2-D
template <>
double RightHandSide<2>::value (const Point<2> &p,
                                const unsigned int /*component*/) const
{
   return 8*M_PI*M_PI*sin(2*M_PI*p[0])*sin(2*M_PI*p[1]);
}

// RHS in 3-D
template <>
double RightHandSide<3>::value (const Point<3> &p,
                                const unsigned int /*component*/) const
{
   return 12*M_PI*M_PI*sin(2*M_PI*p[0])*sin(2*M_PI*p[1])*sin(2*M_PI*p[2]);
}

//------------------------------------------------------------------------------
// Boundary condition
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
   BoundaryValues () : Function<dim>() {}

   double value (const Point<dim>   &p,
                 const unsigned int  component = 0) const override;
};

template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &/*p*/,
                                   const unsigned int /*component*/) const
{
   return 0.0;
}

//------------------------------------------------------------------------------
template <int dim>
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

   Triangulation<dim>     triangulation;
   MappingFE<dim>         mapping;
   FE_SimplexP<dim>       fe;
   DoFHandler<dim>        dof_handler;

   SparsityPattern        sparsity_pattern;
   SparseMatrix<double>   system_matrix;

   Vector<double>         solution;
   Vector<double>         system_rhs;
};


//------------------------------------------------------------------------------
template <int dim>
LaplaceProblem<dim>::LaplaceProblem (int degree) :
mapping(FE_SimplexP<dim>(1)),
fe (degree),
dof_handler (triangulation)
{}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
   Triangulation<dim> tria;
   GridGenerator::subdivided_hyper_cube (tria, 5, 0, 1);
   GridGenerator::convert_hypercube_to_simplex_mesh(tria, triangulation);

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
   system_rhs    = 0;

   QGaussSimplex<dim>  quadrature_formula(fe.degree + 1);
   const RightHandSide<dim> right_hand_side;
   FEValues<dim> fe_values (mapping, fe, quadrature_formula,
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
               cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *  // grad(phi_i)
                                    fe_values.shape_grad (j, q_point) *  // grad(phi_j)
                                    fe_values.JxW (q_point));            // det(J) * w

            cell_rhs(i) += (fe_values.shape_value (i, q_point) *  // phi_i
                            rhs_values[q_point] *                 // f
                            fe_values.JxW (q_point));             // det(J) * w
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
   VectorTools::interpolate_boundary_values (mapping,
                                             dof_handler,
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
   DataOut<dim> data_out;

   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches (mapping, fe.degree);
   std::ofstream output ("solution.vtk");
   data_out.write_vtk (output);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::run ()
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
   int degree = 1;
   LaplaceProblem<2> problem (degree);
   problem.run ();

   return 0;
}
