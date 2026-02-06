/* Write
     -Laplace(p) = f
   as
      J       = grad(p) 
      -div(J) = f
   Weak formulation
     (J, Jt) + (div(Jt),p) = 0
               (div(J),pt) = -(f,pt)
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_raviart_thomas.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_postprocessor.h>

#define DIRICHLET 1
#define NEUMANN   2

#include <fstream>
#include <iostream>

using namespace dealii;

#if PROBLEM == DIRICHLET
#include "dirichlet.h"
const std::string problem = "Dirichlet";
#elif PROBLEM == NEUMANN
#include "neumann.h"
const std::string problem = "Neumann";
#else
#error "Unknown PROBLEM"
#endif


//------------------------------------------------------------------------------
template <int dim>
class Postprocessor : public DataPostprocessor<dim>
{
public:
   Postprocessor();

   void
   evaluate_vector_field(const DataPostprocessorInputs::Vector<dim>& input_data,
                         std::vector<Vector<double>>& computed_quantities) const override
   {
      const std::vector<std::vector<Tensor<1, dim>>>& uh = input_data.solution_gradients;
      for(unsigned int i = 0; i < uh.size(); ++i)
      {
         computed_quantities[i](0) = uh[i][0][0] + uh[i][1][1];
         if(dim == 3) computed_quantities[i](0) += uh[i][2][2];
      }
   }

   std::vector<std::string>
   get_names() const override
   {
      return {"div_j"};
   }

   UpdateFlags
   get_needed_update_flags() const override
   {
      return update_gradients;
   }
};

//------------------------------------------------------------------------------
template <int dim>
Postprocessor<dim>::Postprocessor()
   : DataPostprocessor<dim>()
{}

//------------------------------------------------------------------------------
class ParameterReader : public EnableObserverPointer
{
public:
  ParameterReader(ParameterHandler&);
  void read_parameters(const std::string&);

private:
  void declare_parameters();
  ParameterHandler& prm;
};

//------------------------------------------------------------------------------
ParameterReader::ParameterReader(ParameterHandler& paramhandler)
  : prm(paramhandler)
{}

void
ParameterReader::declare_parameters()
{
  prm.declare_entry("degree",
                    "1",
                    Patterns::Integer(0),
                    "Degree of polynomial");
  prm.declare_entry("nrefine",
                    "4",
                    Patterns::Integer(0),
                    "number of refinements");
  prm.declare_entry("initial_refine",
                    "4",
                    Patterns::Integer(0),
                    "Initial refinement of mesh");
  prm.declare_entry("solver",
                    "schur",
                    Patterns::Selection("schur|umfpack"),
                    "Linear solver");
}

void
ParameterReader::read_parameters(const std::string& parameter_file)
{
  declare_parameters();
  prm.parse_input(parameter_file);
}

//------------------------------------------------------------------------------
template <int dim>
class MixedLaplaceProblem
{
public:
   MixedLaplaceProblem(const unsigned int nrefine,
                       const unsigned int degree, 
                       const unsigned int initial_refine,
                       const std::string  solver);
   void run(std::vector<int>&    ncell,  
            std::vector<double>& h_array, 
            std::vector<int>&    phi_dofs,
            std::vector<int>&    j_dofs, 
            std::vector<double>& phi_error, 
            std::vector<double>& j_error,
            std::vector<double>& d_error,
            std::vector<int>&    phi_iterations, 
            std::vector<int>&    j_iterations,
            std::vector<double>& phi_time,
            std::vector<double>& j_time);

private:
   void make_grid(const unsigned int refine);
   void setup_system();
   void assemble_system();
   void solve_schur(int&    phi_iteration, 
                    int&    j_iteration, 
                    double& phi_time,
                    double& j_time);
   void solve_umfpack(int&    phi_iteration, 
                      int&    j_iteration,
                      double& phi_time, 
                      double& j_time);
   void solve(int&    phi_iteration, 
              int&    j_iteration,
              double& phi_time, 
              double& j_time);
   void compute_errors(double& phi_err,
                       double& j_err,
                       double& d_err);
   void output_results(unsigned int n);
   void refine_grid(unsigned int refine);

   const unsigned int degree;
   const unsigned int nrefine;
   unsigned int initial_refine;
   std::string linear_solver;

   Timer               timer;
   double              h_max;
   Triangulation<dim>  triangulation;
   const FESystem<dim> fe;
   DoFHandler<dim>     dof_handler;

   AffineConstraints<double> constraints;

   BlockSparsityPattern      sparsity_pattern;
   BlockSparseMatrix<double> system_matrix;

   BlockVector<double> solution;
   BlockVector<double> system_rhs;
};

//------------------------------------------------------------------------------
template <int dim>
MixedLaplaceProblem<dim>::MixedLaplaceProblem(
      const unsigned int nrefine,
      const unsigned int degree,
      const unsigned int initial_refine,
      const std::string  solver)
   :
   degree(degree),
   nrefine(nrefine),
   initial_refine(initial_refine),
   linear_solver(solver),
   fe(FE_RaviartThomas<dim>(degree), FE_DGQ<dim>(degree)),
   dof_handler(triangulation)
{}

//------------------------------------------------------------------------------
template <int dim>
void
MixedLaplaceProblem<dim>::make_grid(const unsigned int refine)
{
   triangulation.clear();
   GridGenerator::hyper_cube(triangulation, 0.0, 1.0, false);
   triangulation.refine_global(refine);
}

//------------------------------------------------------------------------------
template<int dim>
void
MixedLaplaceProblem<dim>::setup_system()
{
   timer.start();
   dof_handler.distribute_dofs(fe);

   DoFRenumbering::component_wise(dof_handler);

   const std::vector<types::global_dof_index> dofs_per_component =
      DoFTools::count_dofs_per_fe_component(dof_handler);
   const unsigned int n_c = dofs_per_component[0],
                      n_p = dofs_per_component[dim];

   std::cout << "Number of active cells: " << triangulation.n_active_cells()
             << std::endl
             << "Total number of cells: " << triangulation.n_cells()
             << std::endl
             << "Number of degrees of freedom: " << dof_handler.n_dofs()
             << " (" << n_c << '+' << n_p << ')' << std::endl;

   constraints.clear();

   #if PROBLEM == DIRICHLET
   // p = 0 on all boundary
   // Nothing to do, weakly imposed, just dont add any boundary integral.
   #else // Neumann problem
   // J.n = 0 on all boundaries
   VectorTools::project_boundary_values_div_conforming(dof_handler,
         0,
         Functions::ZeroFunction<dim>(dim),
         types::boundary_id(0),
         constraints);

   // Add mean value constraint on phi
   Vector<double> integral_vector;
   integral_vector.reinit(dof_handler.n_dofs());

   const FEValuesExtractors::Scalar scalar_extractor(dim);
   const ComponentMask scalar_mask = fe.component_mask(scalar_extractor);
   std::vector<bool> bool_boundary_dofs;
   DoFTools::extract_dofs_with_support_on_boundary(dof_handler,
                                                   scalar_mask,
                                                   bool_boundary_dofs,
                                                   {0});

   const IndexSet all_dofs = DoFTools::extract_dofs(dof_handler, scalar_mask);
   types::global_dof_index first_dof = all_dofs.nth_index_in_set(0);

   const QGauss < dim - 1 > quadrature_formula(degree + 2);
   FEFaceValues<dim> face_fe_values(fe, quadrature_formula,
                                    update_values | update_gradients |
                                    update_quadrature_points |
                                    update_JxW_values);

   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
   const unsigned int face_n_q_points = quadrature_formula.size();
   Vector<double> local_rhs(dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

   for(const auto& cell : dof_handler.active_cell_iterators())
   {
      local_rhs = 0.0;
      for(unsigned int f = 0; f < cell->n_faces(); ++f)
      {
         if(cell->face(f)->at_boundary())
         {
            face_fe_values.reinit(cell, f);
            for(unsigned int q_point = 0; q_point < face_n_q_points; ++q_point)
            {
               for(unsigned int i = 0; i < dofs_per_cell; ++i)
               {
                  local_rhs(i) += face_fe_values[scalar_extractor].value(i, q_point) *
                                  face_fe_values.JxW(q_point);
               }
            }
         }
      }
      cell->get_dof_indices(local_dof_indices);
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
      {
         if(local_rhs(i) != 0)
         {
            integral_vector(local_dof_indices[i]) += local_rhs(i);
         }
      }
   }

   std::vector<std::pair<types::global_dof_index, double>> rhs;
   for(const types::global_dof_index i : all_dofs)
      if(i != first_dof && bool_boundary_dofs[i])
         rhs.emplace_back(i, -integral_vector(i) / integral_vector(first_dof));
   constraints.add_constraint(first_dof, rhs);
   #endif

   constraints.close();

   Table<2, DoFTools::Coupling> coupling;
   coupling.reinit(dim + 1, dim + 1);
   coupling.fill(DoFTools::always);
   coupling(dim, dim) = DoFTools::none;

   const std::vector<types::global_dof_index> block_sizes = {n_c, n_p};
   BlockDynamicSparsityPattern                dsp(block_sizes, block_sizes);
   DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp, constraints);
   sparsity_pattern.copy_from(dsp);

   system_matrix.reinit(sparsity_pattern);
   solution.reinit(block_sizes);
   system_rhs.reinit(block_sizes);

   timer.stop();
   // double time_elapsed = timer.last_wall_time();
}

//------------------------------------------------------------------------------
template <int dim>
void
MixedLaplaceProblem<dim>::assemble_system()
{
   timer.start();
   const QGauss<dim> quadrature_formula(degree + 2);
   FEValues<dim>  fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

   const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
   const unsigned int n_q_points      = quadrature_formula.size();

   FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
   Vector<double>     local_rhs(dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
   const PrescribedSolution::RHSFunction<dim> rhs_function;
   const FEValuesExtractors::Vector current(0);
   const FEValuesExtractors::Scalar potential(dim);

   h_max = 0.0;
   double aspect = 1.0;
   for(const auto& cell : dof_handler.active_cell_iterators())
   {
      h_max = std::max(h_max, cell->diameter());
      const double dx = cell->face(2)->measure();
      const double dy = cell->face(0)->measure();
      const double emin = std::min(dx, dy);
      const double emax = std::max(dx, dy);
      aspect = std::max(aspect, emax / emin);

      fe_values.reinit(cell);

      local_matrix = 0.0;
      local_rhs    = 0.0;

      for(unsigned int q = 0; q < n_q_points; ++q)
      {
         const auto rhs = rhs_function.value(fe_values.quadrature_point(q));
         for(unsigned int i = 0; i < dofs_per_cell; ++i)
         {
            const Tensor<1, dim> phi_i_c = fe_values[current].value(i, q);
            const double div_phi_i_c = fe_values[current].divergence(i, q);
            const double phi_i_p     = fe_values[potential].value(i, q);

            for(unsigned int j = 0; j < dofs_per_cell; ++j)
            {
               const Tensor<1, dim> phi_j_c =
                  fe_values[current].value(j, q);
               const double div_phi_j_c =
                  fe_values[current].divergence(j, q);
               const double phi_j_p = fe_values[potential].value(j, q);

               local_matrix(i, j) +=
                  (phi_i_c * phi_j_c
                   + phi_i_p * div_phi_j_c
                   + div_phi_i_c * phi_j_p)
                  * fe_values.JxW(q);
            }

            local_rhs(i) += - phi_i_p * rhs * fe_values.JxW(q);
         }
      }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(local_matrix,
                                             local_rhs,
                                             local_dof_indices,
                                             system_matrix,
                                             system_rhs);
   }
   std::cout << "Largest cell aspect ratio = " << aspect << std::endl;
   timer.stop();
   double time_elapsed = timer.last_wall_time();
   std::cout << "Time for assembly = " << time_elapsed << " sec\n";
}

//------------------------------------------------------------------------------
template <int dim>
void
MixedLaplaceProblem<dim>::solve_schur(int&    phi_iteration, 
                                      int&    j_iteration,
                                      double& phi_time,  
                                      double& j_time)
{
   const auto& M = system_matrix.block(0, 0);
   const auto& B = system_matrix.block(0, 1);

   const auto& F = system_rhs.block(0);
   const auto& G = system_rhs.block(1);

   auto& U = solution.block(0);
   auto& P = solution.block(1);

   const auto op_M = linear_operator(M);
   const auto op_B = linear_operator(B);

   ReductionControl         reduction_control_M(2000, 1.0e-18, 1.0e-10);
   SolverCG<Vector<double>> solver_M(reduction_control_M);
   PreconditionJacobi<SparseMatrix<double>> preconditioner_M;

   preconditioner_M.initialize(M);

   const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M);

   const auto op_S = transpose_operator(op_B) * op_M_inv * op_B;
   const auto op_aS =
      transpose_operator(op_B) * linear_operator(preconditioner_M) * op_B;

   IterationNumberControl   iteration_number_control_aS(30, 1.e-18);
   SolverCG<Vector<double>> solver_aS(iteration_number_control_aS);

   const auto preconditioner_S =
      inverse_operator(op_aS, solver_aS, PreconditionIdentity());

   const auto schur_rhs = transpose_operator(op_B) * op_M_inv * F - G;

   SolverControl            solver_control_S(6000, 1.e-8);
   SolverCG<Vector<double>> solver_S(solver_control_S);

   const auto op_S_inv = inverse_operator(op_S, solver_S, preconditioner_S);

   timer.start();
   P = op_S_inv * schur_rhs;
   timer.stop();
   phi_time = timer.last_wall_time();
   phi_iteration = solver_control_S.last_step();

   timer.start();
   U = op_M_inv * (F - op_B * P);
   constraints.distribute(solution);
   timer.stop();
   j_iteration = reduction_control_M.last_step();
   j_time = timer.last_wall_time();
}

//------------------------------------------------------------------------------
template <int dim>
void
MixedLaplaceProblem<dim>::solve_umfpack(int&    phi_iteration, 
                                        int&    j_iteration,
                                        double& phi_time,  
                                        double& j_time)
{
   timer.start();
   SparseDirectUMFPACK solver;
   solver.initialize(system_matrix);
   solver.vmult(solution, system_rhs);
   constraints.distribute(solution);
   timer.stop();

   phi_iteration = 0;
   j_iteration = 0;
   phi_time = timer.last_wall_time();
   j_time = 0.0;
}

//------------------------------------------------------------------------------
template <int dim>
void
MixedLaplaceProblem<dim>::solve(int&    phi_iteration, 
                                int&    j_iteration,
                                double& phi_time, 
                                double& j_time)
{
   if(linear_solver == "schur")
      solve_schur(phi_iteration, j_iteration, phi_time, j_time);
   else
      solve_umfpack(phi_iteration, j_iteration, phi_time, j_time);

   std::cout << "Time to solve (phi,j,total) = " << phi_time << ", " << j_time
             << ", " << phi_time + j_time << std::endl;
}

//------------------------------------------------------------------------------
template <int dim>
void
MixedLaplaceProblem<dim>::compute_errors(double& phi_err,
                                         double& j_err,
                                         double& d_err)
{
   PrescribedSolution::ExactSolution<dim> exact_solution;
   const QGauss<dim> quadrature(degree + 2);

   {
      const ComponentSelectFunction<dim> scalar_mask(dim, dim + 1);
      Vector<double> cellwise_errors(triangulation.n_active_cells());
      VectorTools::integrate_difference(dof_handler,
                                        solution,
                                        exact_solution,
                                        cellwise_errors,
                                        quadrature,
                                        VectorTools::L2_norm,
                                        &scalar_mask);
      phi_err = VectorTools::compute_global_error(triangulation,
                                                  cellwise_errors,
                                                  VectorTools::L2_norm);
   }

   {
      const ComponentSelectFunction<dim> vector_mask(std::make_pair(0, dim),
                                                     dim + 1);
      Vector<double> cellwise_errors(triangulation.n_active_cells());
      VectorTools::integrate_difference(dof_handler,
                                        solution,
                                        exact_solution,
                                        cellwise_errors,
                                        quadrature,
                                        VectorTools::L2_norm,
                                        &vector_mask);
      j_err = VectorTools::compute_global_error(triangulation,
                                                cellwise_errors,
                                                VectorTools::L2_norm);
   }

   {
      const ComponentSelectFunction<dim> vector_mask(std::make_pair(0, dim),
                                                     dim + 1);
      Vector<double> cellwise_errors(triangulation.n_active_cells());
      VectorTools::integrate_difference(dof_handler,
                                        solution,
                                        exact_solution,
                                        cellwise_errors,
                                        quadrature,
                                        VectorTools::Hdiv_seminorm,
                                        &vector_mask);
      d_err = VectorTools::compute_global_error(triangulation,
                                                cellwise_errors,
                                                VectorTools::L2_norm);
   }
}

//------------------------------------------------------------------------------
template <int dim>
void
MixedLaplaceProblem<dim>::output_results(unsigned int n)
{
   timer.start();
   std::vector<std::string> solution_names(dim, "vector");
   solution_names.emplace_back("scalar");
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
   interpretation(dim,
                  DataComponentInterpretation::component_is_part_of_vector);
   interpretation.push_back(DataComponentInterpretation::component_is_scalar);

   DataOut<dim> data_out;
   data_out.add_data_vector(dof_handler,
                            solution,
                            solution_names,
                            interpretation);

   Postprocessor<dim> div_j;
   data_out.add_data_vector(dof_handler, solution, div_j);

   data_out.build_patches(degree + 1);
   std::ofstream output("solution_" + Utilities::int_to_string(n) + ".vtu");
   data_out.write_vtu(output);

   timer.stop();
   // double time_elapsed = timer.last_wall_time();
}

//------------------------------------------------------------------------------
template <int dim>
void
MixedLaplaceProblem<dim>::run(std::vector<int>&    ncell,
                              std::vector<double>& h_array,
                              std::vector<int>&    phi_dofs,
                              std::vector<int>&    j_dofs,
                              std::vector<double>& phi_error,
                              std::vector<double>& j_error,
                              std::vector<double>& d_error,
                              std::vector<int>&    phi_iterations,
                              std::vector<int>&    j_iterations,
                              std::vector<double>& phi_time,
                              std::vector<double>& j_time)
{
   for(unsigned int i = 0; i < nrefine; ++i)
   {
      timer.reset();

      make_grid(i + initial_refine);
      std::cout << "---------------Grid level = " << i << "-----------------\n";
      setup_system();

      // Solve J,phi
      assemble_system();
      solve(phi_iterations[i], j_iterations[i], phi_time[i], j_time[i]);

      compute_errors(phi_error[i], j_error[i], d_error[i]);
      output_results(i);

      ncell[i] = triangulation.n_active_cells();
      const std::vector<types::global_dof_index> dofs_per_component =
         DoFTools::count_dofs_per_fe_component(dof_handler);
      const unsigned int n_c = dofs_per_component[0],
                         n_p = dofs_per_component[dim];
      phi_dofs[i] = n_p;
      j_dofs[i] = n_c;
      h_array[i] = h_max;
   }
}

//------------------------------------------------------------------------------
int
main()
{
   try
   {
    std::cout << "Solving " << problem << " problem\n";

    ParameterHandler prm;
    ParameterReader param(prm);
    param.read_parameters("input.prm");

    unsigned int degree = prm.get_integer("degree");
    unsigned int nrefine = prm.get_integer("nrefine");
    unsigned int initial_refine = prm.get_integer("initial_refine");
    std::cout << "Using FE = RT(" << degree << "), DG(" << degree << ")\n";

    std::vector<int> ncell(nrefine),  p_dofs(nrefine), j_dofs(nrefine),  
                     p_iterations(nrefine), j_iterations(nrefine);
    std::vector<double> p_error(nrefine),  j_error(nrefine), d_error(nrefine),
                        p_time(nrefine), j_time(nrefine), h_array(nrefine);

    MixedLaplaceProblem<2> mixed_laplace_problem(nrefine, 
                                                 degree, 
                                                 initial_refine,
                                                 prm.get("solver"));
    mixed_laplace_problem.run(ncell, 
                              h_array,
                              p_dofs, 
                              j_dofs, 
                              p_error, 
                              j_error, 
                              d_error,
                              p_iterations, 
                              j_iterations, 
                              p_time, 
                              j_time);

      ConvergenceTable  convergence_table;
      for(unsigned int n = 0; n < nrefine; ++n)
      {
         convergence_table.add_value("cells", ncell[n]);
         convergence_table.add_value("dofs(p)",  p_dofs[n]);
         convergence_table.add_value("dofs(j)",  j_dofs[n]);
         convergence_table.add_value("hmax", h_array[n]);
         convergence_table.add_value("time(p)", p_time[n]);
         convergence_table.add_value("iter(p)", p_iterations[n]);
         convergence_table.add_value("error(p)", p_error[n]);
         convergence_table.add_value("time(j)", j_time[n]);
         convergence_table.add_value("iter(j)", j_iterations[n]);
         convergence_table.add_value("error(j)", j_error[n]);
         convergence_table.add_value("error(d)", d_error[n]);
      }

      convergence_table.set_precision("error(j)", 3);
      convergence_table.set_scientific("error(j)", true);
      convergence_table.set_precision("error(p)", 3);
      convergence_table.set_scientific("error(p)", true);
      convergence_table.set_precision("error(d)", 3);
      convergence_table.set_scientific("error(d)", true);

      convergence_table.set_tex_caption("cells", "cells");
      convergence_table.set_tex_caption("dofs(p)", "dofs(p)");
      convergence_table.set_tex_caption("dofs(j)", "dofs(j)");
      convergence_table.set_tex_caption("hmax", "hmax");
      convergence_table.set_tex_caption("time(p)", "time(p)");
      convergence_table.set_tex_caption("iter(p)", "iter(p)");
      convergence_table.set_tex_caption("error(p)", "error(p)");
      convergence_table.set_tex_caption("time(j)", "time(j)");
      convergence_table.set_tex_caption("iter(j)", "iter(j)");
      convergence_table.set_tex_caption("error(j)", "error(j)");
      convergence_table.set_tex_caption("error(d)", "error(div)");

      convergence_table.set_tex_format("cells", "r");
      convergence_table.set_tex_format("dofs(p)", "r");
      convergence_table.set_tex_format("dofs(j)", "r");
      convergence_table.set_tex_format("hmax",  "r");
      convergence_table.set_tex_format("time(p)", "r");
      convergence_table.set_tex_format("iter(p)",  "r");
      convergence_table.set_tex_format("time(j)", "r");
      convergence_table.set_tex_format("iter(j)",  "r");

      convergence_table.evaluate_convergence_rates
      ("error(j)", ConvergenceTable::reduction_rate_log2);
      convergence_table.evaluate_convergence_rates
      ("error(p)", ConvergenceTable::reduction_rate_log2);
      convergence_table.evaluate_convergence_rates
      ("error(d)", ConvergenceTable::reduction_rate_log2);

      std::cout << std::endl;
      convergence_table.write_text(std::cout);

      std::ofstream tex_file(problem + ".tex");
      convergence_table.write_tex(tex_file);
      std::ofstream text_file(problem + ".txt");
      convergence_table.write_text(text_file);
   }

   catch(std::exception& exc)
   {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
   }
   catch(...)
   {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
   }

   return 0;
}
