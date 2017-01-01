/*
   Solves
      -Laplace(u) = 1 in (0,1) x (0,1)
   with zero Dirichlet boundary conditions.

   All assesmbly operations are done using library functions.
*/
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
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

int main()
{
   const unsigned int degree = 1;
   Triangulation<2> triangulation;
   GridGenerator::hyper_cube (triangulation);
   triangulation.refine_global (4);

   const FE_Q<2> fe(degree);

   DoFHandler<2> dof_handler (triangulation);
   dof_handler.distribute_dofs (fe);

   DynamicSparsityPattern dsp(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, dsp);
   SparsityPattern sparsity_pattern;
   sparsity_pattern.copy_from(dsp);

   SparseMatrix<double> system_matrix;
   system_matrix.reinit (sparsity_pattern);
   MatrixCreator::create_laplace_matrix (dof_handler, 
                                         QGauss<2>(2*degree), 
                                         system_matrix);

   Vector<double> system_rhs;
   system_rhs.reinit (dof_handler.n_dofs());
   VectorTools::create_right_hand_side (dof_handler,
                                        QGauss<2>(2*degree),
                                        ConstantFunction<2>(1.0),
                                        system_rhs);

   std::map<unsigned int,double> boundary_values;
   VectorTools::interpolate_boundary_values (dof_handler,
                                             0,
                                             ZeroFunction<2>(),
                                             boundary_values);

   Vector<double> solution;
   solution.reinit (dof_handler.n_dofs());

   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);

   SolverControl solver_control (1000, 1.0e-12);
   SolverCG<>    solver (solver_control);
   solver.solve (system_matrix,
                 solution,
                 system_rhs,
                 PreconditionIdentity());

   DataOut<2> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches (degree);

   std::ofstream output ("solution.vtk");
   data_out.write_vtk (output);
}
