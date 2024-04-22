/*
   Solves
      -Laplace(u) = 1 in (0,1) x (0,1)
   with zero Dirichlet boundary conditions.

   All assembly operations are done using library functions.

   Use Cartesian/quad grid.
*/
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>

using namespace dealii;
const int dim = 2;

int main()
{
   const unsigned int degree = 1;

   // Create mesh
   Triangulation<dim> triangulation;
   GridGenerator::hyper_cube (triangulation);
   triangulation.refine_global (4);

   const MappingQ1<dim> mapping;

   // Declare FE space
   const FE_Q<dim> fe(degree);

   // Number the dofs
   DoFHandler<dim> dof_handler (triangulation);
   dof_handler.distribute_dofs (fe);

   // Create sparsity pattern
   DynamicSparsityPattern dsp(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, dsp);
   SparsityPattern sparsity_pattern;
   sparsity_pattern.copy_from(dsp);

   // Compute matrix
   SparseMatrix<double> system_matrix;
   system_matrix.reinit (sparsity_pattern);
   MatrixCreator::create_laplace_matrix (mapping,
                                         dof_handler, 
                                         QGauss<dim>(2*degree),
                                         system_matrix);

   // Compute rhs vector
   Vector<double> system_rhs;
   system_rhs.reinit (dof_handler.n_dofs());
   VectorTools::create_right_hand_side (mapping,
                                        dof_handler,
                                        QGauss<dim>(2*degree),
                                        Functions::ConstantFunction<dim>(1.0),
                                        system_rhs);

   // Solution vector
   Vector<double> solution;
   solution.reinit (dof_handler.n_dofs());

   // Apply boundary conditions
   std::map<types::global_dof_index,double> boundary_values;
   VectorTools::interpolate_boundary_values (mapping,
                                             dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);

   // Solve the problem
   SolverControl solver_control (1000, 1.0e-12);
   SolverCG<>    solver (solver_control);
   solver.solve (system_matrix,
                 solution,
                 system_rhs,
                 PreconditionIdentity());

   // Save solution to file
   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches (mapping, degree);

   std::ofstream output ("solution.vtk");
   data_out.write_vtk (output);
}
