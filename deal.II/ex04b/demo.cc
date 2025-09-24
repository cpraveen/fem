/* 
 Solve 2d Poisson equation
    -Laplace(u) = f(x) in (0,1)x(0,1)
 Exact solution is 
    u = sin(2*pi*x) * sin(2*pi*y)
 f is obtained from exact solution.
 Boundary condition is dirichlet and taken from exact solution.
 Perform grid convergence test.
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
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <fstream>
#include <iostream>

using namespace dealii;

//------------------------------------------------------------------------------
template <int dim>
class ExactSolution : public Function<dim>
{
public:
   ExactSolution () : Function<dim>() {}
   
   double value (const Point<dim>   &p,
                 const unsigned int  component = 0) const override;
   Tensor<1,dim> gradient (const Point<dim>   &p,
                           const unsigned int  component = 0) const override;
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
   
   double value (const Point<dim>   &p,
                 const unsigned int  component = 0) const override;
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
   
   double value (const Point<dim>   &p,
                 const unsigned int  component = 0) const override;
};

template <int dim>
double BoundaryValues<dim>::value (const Point<dim>& p,
                                   const unsigned int /*component*/) const
{
   ExactSolution<dim> exact_solution;
   return exact_solution.value(p);
}

//------------------------------------------------------------------------------
template <int dim>
class LaplaceProblem
{
public:
   LaplaceProblem (int degree);
   void run (bool refine, unsigned int &ncell, unsigned int &ndofs, 
             double &L2_error, double &H1_error);
   
private:
   using PTriangulation = parallel::distributed::Triangulation<dim>;
   using PMatrix = PETScWrappers::MPI::SparseMatrix;
   using PVector = PETScWrappers::MPI::Vector;
   using Constraints = AffineConstraints<double>;

   void make_grid();
   void make_dofs ();
   void assemble_system ();
   void solve ();
   void output_results () const;
   void compute_error (double &L2_error, double &H1_error) const;

   MPI_Comm               mpi_comm;
   const unsigned int     mpi_rank;
   ConditionalOStream     pcout;

   PTriangulation         triangulation;
   FE_Q<dim>              fe;
   DoFHandler<dim>        dof_handler;
   Constraints            constraints;
   
   PMatrix                system_matrix;
   PVector                solution;
   PVector                system_rhs;
};


//------------------------------------------------------------------------------
template <int dim>
LaplaceProblem<dim>::LaplaceProblem(int degree)
    : mpi_comm(MPI_COMM_WORLD),
      mpi_rank(Utilities::MPI::this_mpi_process(mpi_comm)),
      pcout(std::cout, mpi_rank==0),
      triangulation(mpi_comm),
      fe(degree),
      dof_handler(triangulation)
{}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::make_grid ()
{
   GridGenerator::hyper_cube (triangulation, 0, 1);
   triangulation.refine_global(5);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::make_dofs ()
{
   dof_handler.distribute_dofs (fe);
   
   pcout << "   Number of active cells: "
         << triangulation.n_global_active_cells()
         << std::endl;
   pcout << "   Number of degrees of freedom: "
         << dof_handler.n_dofs()
         << std::endl;

   const auto &locally_owned_dofs = dof_handler.locally_owned_dofs();
   IndexSet locally_relevant_dofs = 
      DoFTools::extract_locally_relevant_dofs(dof_handler);

   constraints.clear();
   constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
   VectorTools::interpolate_boundary_values(dof_handler,
                                            0,
                                            BoundaryValues<dim>(),
                                            constraints);
   constraints.close();

   DynamicSparsityPattern sparsity_pattern(locally_relevant_dofs);
   DoFTools::make_sparsity_pattern(dof_handler, 
                                   sparsity_pattern,
                                   constraints, 
                                   false);
   SparsityTools::distribute_sparsity_pattern(sparsity_pattern,
                                              locally_owned_dofs,
                                              mpi_comm,
                                              locally_relevant_dofs);

   system_matrix.reinit(locally_owned_dofs,
                        locally_owned_dofs,
                        sparsity_pattern,
                        mpi_comm);
   solution.reinit(locally_owned_dofs, 
                   locally_relevant_dofs, 
                   mpi_comm);
   system_rhs.reinit(locally_owned_dofs, 
                     mpi_comm);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
   system_matrix = 0;
   system_rhs = 0;

   QGauss<dim>  quadrature_formula(fe.degree + 1);
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
   if(cell->is_locally_owned())
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
      constraints.distribute_local_to_global(cell_matrix,
                                             cell_rhs,
                                             local_dof_indices,
                                             system_matrix,
                                             system_rhs);
   }

   system_matrix.compress(VectorOperation::add);
   system_rhs.compress(VectorOperation::add);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::solve ()
{
   SolverControl solver_control(solution.size(), 1e-8 * system_rhs.l2_norm());
   PETScWrappers::SolverCG cg(solver_control);
   PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);
   PVector distributed_solution(dof_handler.locally_owned_dofs(),
                                mpi_comm);
   cg.solve(system_matrix, distributed_solution, system_rhs, preconditioner);
   constraints.distribute(distributed_solution);
   solution = distributed_solution;

   pcout
   << "   " << solver_control.last_step()
   << " CG iterations needed to obtain convergence."
   << std::endl;
}

//------------------------------------------------------------------------------
// TODO: Fix this to run in parallel
template <int dim>
void LaplaceProblem<dim>::output_results () const
{
   static int step = 0;
   static std::vector<std::pair<double, std::string>> pvtu_files;
   const unsigned int n_digits_for_counter = 2;
   PVector solution_error(dof_handler.locally_owned_dofs(), mpi_comm);
   VectorTools::interpolate(dof_handler,
                            ExactSolution<dim>(),
                            solution_error);
   solution_error -= solution;

   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.add_data_vector (solution_error, "error");
   data_out.build_patches (fe.degree);

   data_out.write_vtu_with_pvtu_record("./",
                                      "sol",
                                      step,
                                      mpi_comm,
                                      n_digits_for_counter);

   std::string fname = "sol_" + 
                       Utilities::int_to_string(step,n_digits_for_counter) +
                       ".pvtu";
   pvtu_files.emplace_back(step, fname);
   std::ofstream pvd_file("sol.pvd");
   DataOutBase::write_pvd_record(pvd_file, pvtu_files);

   pcout << "   Wrote to file " << fname << std::endl;
   ++step;
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::compute_error (double &L2_error, double &H1_error) const
{
   // compute L2 error in solution
   Vector<double> difference_per_cell (triangulation.n_active_cells());
   VectorTools::integrate_difference (dof_handler,
                                      solution,
                                      ExactSolution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree+3),
                                      VectorTools::L2_norm);
   L2_error = VectorTools::compute_global_error(triangulation,
                                                difference_per_cell,
                                                VectorTools::L2_norm);

   // compute L2 error in gradient
   VectorTools::integrate_difference (dof_handler,
                                      solution,
                                      ExactSolution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree+3),
                                      VectorTools::H1_seminorm);
   H1_error = VectorTools::compute_global_error(triangulation,
                                                difference_per_cell,
                                                VectorTools::L2_norm);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::run (bool          refine,
                               unsigned int &ncell, 
                               unsigned int &ndofs, 
                               double       &L2_error, 
                               double       &H1_error)
{
   if(refine == false)
   {
      make_grid();
   }
   else
   {
      triangulation.refine_global();
   }
   make_dofs();
   assemble_system ();
   solve ();
   output_results ();
   compute_error (L2_error, H1_error);

   ncell = triangulation.n_active_cells ();
   ndofs = dof_handler.n_dofs ();
}

//------------------------------------------------------------------------------
int main(int argc, char **argv)
{
   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
   const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
   deallog.depth_console(0);
   int degree = 1;
   ConvergenceTable  convergence_table;   
   LaplaceProblem<2> problem (degree);
   for(unsigned int n=0; n<5; ++n)
   {
      unsigned int ncell, ndofs;
      double L2_error, H1_error;
      bool refine = (n == 0) ? false : true;
      problem.run (refine, ncell, ndofs, L2_error, H1_error);

      convergence_table.add_value("cells", ncell);
      convergence_table.add_value("dofs",  ndofs);
      convergence_table.add_value("L2",    L2_error);
      convergence_table.add_value("H1",    H1_error);
      if(rank == 0)
         std::cout << "--------------------------------------------------------\n";
   }

   if(rank != 0) return 0;

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
