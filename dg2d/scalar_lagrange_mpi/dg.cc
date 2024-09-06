/*
 * Solves 2d linear advection eqn
 *    u_t + div(a(x)u) = 0
 * on general quad grids, using QGauss Lagrange basis functions
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/distributed/tria.h>

#include <iostream>
#include <fstream>
#include <cmath>

// #define FORCE_USE_OF_TRILINOS

// Use petsc by default, to use trilinos uncomment above line
namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
    !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
   using namespace dealii::LinearAlgebraPETSc;
#elif defined(DEAL_II_WITH_TRILINOS)
   using namespace dealii::LinearAlgebraTrilinos;
#else
#error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
}

const double a_rk[] = {0.0, 3.0/4.0, 1.0/3.0};
const double b_rk[] = {1.0, 1.0/4.0, 2.0/3.0};

enum TestCase {expo, square, circ};

using namespace dealii;

//------------------------------------------------------------------------------
// Speed of advection
//------------------------------------------------------------------------------
template <int dim>
void advection_speed(const Point<dim>& p, Point<dim>& v)
{
   v(0) = -p(1);
   v(1) =  p(0);
}
//------------------------------------------------------------------------------
// Upwind flux
//------------------------------------------------------------------------------
double numerical_flux(double vel, double ul, double ur)
{
   if(vel > 0)
      return vel * ul;
   else
      return vel * ur;
}
//------------------------------------------------------------------------------
// Initial condition function class
//------------------------------------------------------------------------------
template <int dim>
class InitialCondition: public Function<dim>
{
public:
   InitialCondition (TestCase test_case)
   :
   test_case (test_case)
   {};
   void value_list (const std::vector<Point<dim>> &points,
                    std::vector<double> &values,
                    const unsigned int component=0) const override;
private:
   TestCase test_case;
};

// Computes boundary condition value at a list of boundary points
template <int dim>
void InitialCondition<dim>::value_list(const std::vector<Point<dim>> &points,
                                       std::vector<double> &values,
                                       const unsigned int) const
{
   Assert(values.size()==points.size(),
          ExcDimensionMismatch(values.size(),points.size()));

   if(test_case == expo)
      for (unsigned int i=0; i<values.size(); ++i)
      {
         double r2 = std::pow(points[i](0)-0.5, 2.0) + std::pow(points[i](1), 2.0);
         values[i] = 1.0 + exp(-50.0*r2);
      }
   else if(test_case == circ)
      for (unsigned int i=0; i<values.size(); ++i)
      {
         double r2 = std::pow(points[i](0)-0.5, 2.0) + std::pow(points[i](1), 2.0);
         if(r2 < 0.25*0.25)
            values[i] = 2.0;
         else
            values[i] = 1.0;
      }
   else if(test_case == square)
      for (unsigned int i=0; i<values.size(); ++i)
      {
         const Point<dim>& p = points[i];
         if(std::fabs(p(0)-0.5) < 0.20 && std::fabs(p(1)) < 0.20)
            values[i] = 2.0;
         else
            values[i] = 1.0;
      }
   else
   {
      for (unsigned int i=0; i<values.size(); ++i)
         values[i] = 1.0;
      AssertThrow(false, ExcMessage("Unknown test case"));
   }
}

//------------------------------------------------------------------------------
// Boundary condition function class
//------------------------------------------------------------------------------
template <int dim>
class BoundaryValues: public Function<dim>
{
  public:
    BoundaryValues () {};
    void value_list (const std::vector<Point<dim>> &points,
			            std::vector<double> &values,
			            const unsigned int component=0) const override;
};

// Computes boundary condition value at a list of boundary points
template <int dim>
void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
				                         std::vector<double> &values,
				                         const unsigned int) const
{
   Assert(values.size()==points.size(),
          ExcDimensionMismatch(values.size(),points.size()));

   for (unsigned int i=0; i<values.size(); ++i)
   {
      values[i]=1.0;
   }
}

//------------------------------------------------------------------------------
// Class for integrating rhs using MeshWorker
//------------------------------------------------------------------------------
template <int dim>
class RHSIntegrator
{
   public:
      RHSIntegrator (const DoFHandler<dim>& dof_handler)
         : dof_info (dof_handler) {};

      MeshWorker::IntegrationInfoBox<dim> info_box;
      MeshWorker::DoFInfo<dim> dof_info;
      MeshWorker::Assembler::ResidualSimple<LA::MPI::Vector>//< Vector<double> >
         assembler;
};

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim>
class Step12
{
public:
   Step12 (unsigned int degree,
           TestCase     test_case);
   void run ();

private:
   void setup_system ();
   void assemble_mass_matrix ();
   void set_initial_condition ();
   void setup_mesh_worker (RHSIntegrator<dim>&);
   void assemble_rhs (RHSIntegrator<dim>&);
   void compute_dt ();
   void solve ();
   void compute_min_max ();
   void output_results (double time);

   MPI_Comm 					  mpi_communicator;

   parallel::distributed::Triangulation<dim>   triangulation;

   const MappingQ<dim>  mapping;

   unsigned int         degree;
   FE_DGQArbitraryNodes<dim>  fe;
   DoFHandler<dim>      dof_handler;
   DoFHandler<dim>      dof_handler_cell;

   IndexSet 		      locally_owned_dofs;
   IndexSet 		      locally_relevant_dofs;

   LA::MPI::Vector      mass_matrix;
   LA::MPI::Vector      solution;
   LA::MPI::Vector      solution_old;
   LA::MPI::Vector      right_hand_side;
   double               dt;
   double               cfl;
   TestCase             test_case;
   double               sol_min, sol_max;
   double               h_min, h_max;

   typedef MeshWorker::DoFInfo<dim> DoFInfo;
   typedef MeshWorker::IntegrationInfo<dim> CellInfo;

   static void integrate_cell_term (DoFInfo& dinfo, CellInfo& info);
   static void integrate_boundary_term (DoFInfo& dinfo, CellInfo& info);
   static void integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
                                    CellInfo& info1, CellInfo& info2);
   ConditionalOStream 	pcout;
   TimerOutput          computing_timer;

};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
Step12<dim>::Step12 (unsigned int degree,
                     TestCase     test_case)
      :

      mpi_communicator (MPI_COMM_WORLD),
      triangulation (mpi_communicator),
      mapping (degree),
      degree (degree),
      fe (QGauss<1>(degree+1)),
      dof_handler (triangulation),
      dof_handler_cell (triangulation),
      test_case (test_case),
      pcout (std::cout,
             (Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),

      computing_timer (pcout,
                       TimerOutput::summary,
                       TimerOutput::wall_times)
{
   cfl = 0.9/(2.0*degree + 1.0);

      pcout << "Number of threads = "
            << MultithreadInfo::n_threads () << std::endl;
      if(MultithreadInfo::is_running_single_threaded())
         pcout << "Running with single thread\n";
      else
         pcout << "Running with multiple threads\n";
}

//------------------------------------------------------------------------------
// Make dofs and allocate memory
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::setup_system ()
{
   TimerOutput::Scope t(computing_timer, "setup");

   pcout << "Allocating memory ...\n";

   dof_handler.distribute_dofs (fe);

   locally_owned_dofs = dof_handler.locally_owned_dofs ();
   DoFTools::extract_locally_relevant_dofs (dof_handler,
                                            locally_relevant_dofs);

   solution.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
   solution_old.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
   right_hand_side.reinit (locally_owned_dofs, mpi_communicator);

   mass_matrix.reinit (locally_owned_dofs, mpi_communicator);

   assemble_mass_matrix ();

   pcout << "Number of active cells:       "
         << triangulation.n_global_active_cells()
         << std::endl;

   pcout << "Number of degrees of freedom: "
         << dof_handler.n_dofs()
         << std::endl;

}

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// Basis functions are based on QGauss and we use same rule to compute mass
// matrix, which makes it diagonal.
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::assemble_mass_matrix ()
{
   TimerOutput::Scope t(computing_timer, "mass matrix");

   pcout << "Constructing mass matrix ...\n";

   QGauss<dim>  quadrature_formula(fe.degree+1);

   FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                            update_values | update_JxW_values);

   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();

   Vector<double>   cell_matrix (dofs_per_cell);
   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   mass_matrix = 0;

   // Cell iterator
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (unsigned int c = 0; cell!=endc; ++cell, ++c)
   if(cell->is_locally_owned())
   {
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);

      cell_matrix = 0.0;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
               cell_matrix(i) += fe_values.shape_value (i, q_point) *
                                 fe_values.shape_value (i, q_point) *
                                 fe_values.JxW (q_point);

      mass_matrix.add(local_dof_indices,
                      cell_matrix);

   }
   mass_matrix.compress(VectorOperation::add);
}

//------------------------------------------------------------------------------
// Project initial condition by L2 projection
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::set_initial_condition ()
{
   TimerOutput::Scope t(computing_timer, "initial condition");

   pcout << "Setting initial condition ...";

   QGauss<dim> quadrature_formula (fe.degree+1);
   const unsigned int n_q_points = quadrature_formula.size();

   FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                            update_values |
                            update_quadrature_points |
                            update_JxW_values);

   // Multiply by inverse mass matrix
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   std::vector<double> initial_values (n_q_points);
   std::vector<double> cell_vector(dofs_per_cell);
   InitialCondition<dim> initial_condition(test_case);

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      cell->get_dof_indices (local_dof_indices);
      fe_values.reinit (cell);
      initial_condition.value_list (fe_values.get_quadrature_points(),
                                    initial_values);
      for(unsigned int i=0; i<dofs_per_cell; ++i)
      {
         cell_vector[i] = 0;
         for(unsigned int q=0; q<n_q_points; ++q)
            cell_vector[i] += initial_values[q] *
                              fe_values.shape_value(i,q) *
                              fe_values.JxW(q);
         cell_vector[i] /= mass_matrix(local_dof_indices[i]);
      }

      right_hand_side.set(local_dof_indices, cell_vector);
   }

   right_hand_side.compress (VectorOperation::insert);
   solution = right_hand_side;
   pcout << " Done\n";
}

//------------------------------------------------------------------------------
// Create mesh worker for integration
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::setup_mesh_worker (RHSIntegrator<dim>& rhs_integrator)
{
   pcout << "Setting up mesh worker ...\n";

   MeshWorker::IntegrationInfoBox<dim>& info_box = rhs_integrator.info_box;
   //MeshWorker::DoFInfo<dim>& dof_info = rhs_integrator.dof_info;
   MeshWorker::Assembler::ResidualSimple<LA::MPI::Vector>&
      assembler = rhs_integrator.assembler;

   const unsigned int n_gauss_points = fe.degree+1;
   info_box.initialize_gauss_quadrature(n_gauss_points,
                                        n_gauss_points,
                                        n_gauss_points);

   // Add solution vector to info_box
   AnyData solution_data;
   solution_data.add<LA::MPI::Vector*> (&solution, "solution");
   info_box.cell_selector.add     ("solution", true, false, false);
   info_box.boundary_selector.add ("solution", true, false, false);
   info_box.face_selector.add     ("solution", true, false, false);

   info_box.initialize_update_flags ();
   info_box.add_update_flags_all      (update_quadrature_points);
   info_box.add_update_flags_cell     (update_gradients);
   info_box.add_update_flags_boundary (update_values);
   info_box.add_update_flags_face     (update_values);

   info_box.initialize (fe, mapping, solution_data, LA::MPI::Vector());

   // Attach rhs vector to assembler
   AnyData rhs;
   rhs.add<LA::MPI::Vector*> (&right_hand_side, "RHS");

   assembler.initialize (rhs);
}

//------------------------------------------------------------------------------
// Compute time-step
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::compute_dt ()
{
   TimerOutput::Scope t(computing_timer, "dt");

   pcout << "Computing local time-step ...\n";

   dt = 1.0e20;

   // Cell iterator
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      double h = cell->diameter () / std::sqrt(2.0);
      const Point<dim>& cell_center = cell->center();
      Point<dim> beta;
      advection_speed(cell_center, beta);
      double dt_cell = 1.0 / (std::fabs(beta(0))/h + std::fabs(beta(1))/h);
      dt = std::min (dt, dt_cell);
   }

   dt *= cfl;
   dt  = -dt;
   dt = -Utilities::MPI::max (dt, mpi_communicator);
   pcout << "      dt = " << dt << std::endl;
}

//------------------------------------------------------------------------------
// Assemble rhs of the problem
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::assemble_rhs (RHSIntegrator<dim>& rhs_integrator)
{
   TimerOutput::Scope t(computing_timer, "assemble");

   right_hand_side = 0.0;

   typedef
   FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
   CellFilter;

   MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>,
                    MeshWorker::IntegrationInfoBox<dim>>
      (CellFilter (IteratorFilters::LocallyOwnedCell(),dof_handler.begin_active()),
       CellFilter (IteratorFilters::LocallyOwnedCell(),dof_handler.end()),
       rhs_integrator.dof_info,
       rhs_integrator.info_box,
       &Step12<dim>::integrate_cell_term,
       &Step12<dim>::integrate_boundary_term,
       &Step12<dim>::integrate_face_term,
       rhs_integrator.assembler);

   right_hand_side.compress (VectorOperation::add);

   // Multiply by inverse mass matrix
   const unsigned int dofs_per_cell = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
         right_hand_side(local_dof_indices[i]) /= mass_matrix(local_dof_indices[i]);
   }

   right_hand_side.compress (VectorOperation::insert);

}

//------------------------------------------------------------------------------
// Compute cell integral
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_cell_term (DoFInfo& dinfo, CellInfo& info)
{
   const FEValuesBase<dim>& fe_v  = info.fe_values();
   const std::vector<double>& sol = info.values[0][0];
   Vector<double>& local_vector   = dinfo.vector(0).block(0);
   const std::vector<double>& JxW = fe_v.get_JxW_values ();

   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      advection_speed(fe_v.quadrature_point(point), beta);

      for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
         local_vector(i) += beta *
                            fe_v.shape_grad(i,point) *
                            sol[point] *
                            JxW[point];
   }
}

//------------------------------------------------------------------------------
// Compute boundary integral
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_boundary_term (DoFInfo& dinfo, CellInfo& info)
{
   const FEValuesBase<dim>& fe_v  = info.fe_values();
   const std::vector<double>& sol = info.values[0][0];

   Vector<double>& local_vector = dinfo.vector(0).block(0);

   const std::vector<double>& JxW = fe_v.get_JxW_values ();
   const std::vector<Tensor<1,dim>>& normals = fe_v.get_normal_vectors ();

   std::vector<double> g(fe_v.n_quadrature_points);

   static BoundaryValues<dim> boundary_function;
   boundary_function.value_list (fe_v.get_quadrature_points(), g);

   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      advection_speed(fe_v.quadrature_point(point), beta);

      const double beta_n = beta * normals[point];
      const double flux = numerical_flux(beta_n, sol[point], g[point]);
      for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
         local_vector(i) -= flux *
                            fe_v.shape_value(i,point) *
                            JxW[point];
   }
}

//------------------------------------------------------------------------------
// Compute integral over internal faces
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::integrate_face_term (DoFInfo& dinfo1, DoFInfo& dinfo2,
				                           CellInfo& info1, CellInfo& info2)
{
   const FEValuesBase<dim>& fe_v          = info1.fe_values();
   const FEValuesBase<dim>& fe_v_neighbor = info2.fe_values();

   const std::vector<double>& sol1 = info1.values[0][0];
   const std::vector<double>& sol2 = info2.values[0][0];

   Vector<double>& local_vector1 = dinfo1.vector(0).block(0);
   Vector<double>& local_vector2 = dinfo2.vector(0).block(0);

   const std::vector<double>& JxW = fe_v.get_JxW_values ();
   const std::vector<Tensor<1,dim>>& normals = fe_v.get_normal_vectors ();

   for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
   {
      Point<dim> beta;
      advection_speed(fe_v.quadrature_point(point), beta);

      const double beta_n = beta * normals[point];
      const double flux = numerical_flux(beta_n, sol1[point], sol2[point]);
         for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
            local_vector1(i) -= flux *
                                fe_v.shape_value(i,point) *
                                JxW[point];

         for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
            local_vector2(k) += flux *
                                fe_v_neighbor.shape_value(k,point) *
                                JxW[point];
   }
}

//------------------------------------------------------------------------------
// Compute min and max of cell average values
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::compute_min_max ()
{
   std::vector<unsigned int> dof_indices (fe.dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

   sol_min =  1.0e20; sol_max = -1.0e20;
   h_min   =  1.0e20; h_max   = -1.0e20;

   for(; cell != endc; ++cell)
   if(cell->is_locally_owned())
   {
      cell->get_dof_indices (dof_indices);
      double avg = solution(dof_indices[0]);
      sol_min = std::min( sol_min, avg );
      sol_max = std::max( sol_max, avg );
      h_min = std::min(h_min, cell->diameter()/std::sqrt(2));
      h_max = std::max(h_max, cell->diameter()/std::sqrt(2));
   }

   const double local_min[2] = {-sol_min, -h_min};
   const double local_max[2] = { sol_max,  h_max};

   double global_min[2], global_max[2];
   Utilities::MPI::max (local_min, mpi_communicator, global_min);
   Utilities::MPI::max (local_max, mpi_communicator, global_max);

   sol_min = -global_min[0]; h_min = -global_min[1];
   sol_max =  global_max[0]; h_max =  global_max[1];

}

//------------------------------------------------------------------------------
// Solve the problem to convergence by RK time integration
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::solve ()
{
   RHSIntegrator<dim> rhs_integrator (dof_handler);
   setup_mesh_worker (rhs_integrator);
   compute_dt ();

   pcout << "Solving by RK ...\n";

   double final_time = 2.0*M_PI;
   unsigned int iter = 0;
   double time = 0;
   while (time < final_time)
   {
      // We want to reach final_time exactly
      if(time + dt > final_time) dt = final_time - time;

      solution_old = solution;

      // 3-stage RK scheme
      for(unsigned int r=0; r<3; ++r)
      {
         assemble_rhs (rhs_integrator);
         {
            TimerOutput::Scope t(computing_timer, "update");

            // right_hand_side = dt*right_hand_side + solution
            right_hand_side.sadd (dt, solution);
            // right_hand_side = b_rk*right_hand_side + a_rk*solution_old
            right_hand_side.sadd (b_rk[r], a_rk[r], solution_old);
            solution = right_hand_side;
         }
      }
      compute_min_max();

      ++iter; time += dt;

      pcout << "It=" << iter
            << ", t= " << time
            << ", dt= " << dt
            << ", min,max u= " << sol_min << "  " << sol_max
            << std::endl;
      if(std::fmod(iter,100)==0 || std::fabs(time-final_time) < 1.0e-14)
         output_results(time);
   }

}

//------------------------------------------------------------------------------
// Save results to file
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::output_results (double time)
{
   TimerOutput::Scope t(computing_timer, "output");

   static unsigned int cycle = 0;
   static std::vector<std::vector<std::string>> all_files;

   {
      // Output of the polynomial solution
      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "u");

      Vector<float> subdomain (triangulation.n_active_cells());
      for (unsigned int i=0; i<subdomain.size();++i)
         subdomain(i) = triangulation.locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");

      data_out.build_patches(mapping, fe.degree);

      DataOutBase::VtkFlags flags;
      flags.time = time;
      flags.cycle = cycle;
      flags.write_higher_order_cells = false;
      data_out.set_flags(flags);
      data_out.write_vtu_with_pvtu_record("./", "sol", cycle,
                                          MPI_COMM_WORLD, 4);

      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
         std::vector<std::string> filenames;
         for (unsigned int i=0;
              i<Utilities::MPI::n_mpi_processes(mpi_communicator);
              ++i)
            filenames.push_back ("sol_" +
                                 Utilities::int_to_string (cycle, 4) +
                                 "." +
                                 Utilities::int_to_string (i) +
                                 ".vtu");
         all_files.push_back (filenames);
         std::ofstream visit_output ("solution.visit");
         DataOutBase::write_visit_record(visit_output, all_files);
      }

      static std::vector<std::pair<double, std::string>> pvtu_files;
      std::string pvtu_filename = "sol_" +
                                  Utilities::int_to_string(cycle, 4) +
                                  ".pvtu";
      pvtu_files.emplace_back(time, pvtu_filename);
      std::ofstream pvd_file("solution.pvd");
      DataOutBase::write_pvd_record(pvd_file, pvtu_files);
   }

   pcout << "Wrote solution at time,cycle = " << time << " " << cycle << std::endl;

   ++cycle;
}

//------------------------------------------------------------------------------
// Actual computation starts from here
//------------------------------------------------------------------------------
template <int dim>
void Step12<dim>::run ()
{
   int grid = 2;

   if(grid == 0)
   {
      GridGenerator::subdivided_hyper_cube (triangulation,100,-1.0,+1.0);
   }
   else if(grid == 1)
   {
      GridGenerator::hyper_shell(triangulation,Point<dim>(0.0,0.0), 0.1, 1.0);
      triangulation.refine_global(5);
   }
   else if(grid == 2)
   {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation (triangulation);
      std::ifstream gfile ("annulus.msh");
      AssertThrow(gfile.is_open(), ExcMessage("Grid file not found"));
      grid_in.read_msh(gfile);

      const Point<dim> center(0.0, 0.0);
      const SphericalManifold<dim> manifold(center);
      triangulation.set_all_manifold_ids(0);
      triangulation.set_all_manifold_ids_on_boundary(1);
      triangulation.set_manifold(1, manifold);
   }

   setup_system ();
   set_initial_condition ();
   output_results(0);
   solve ();

   computing_timer.print_summary ();
   computing_timer.reset ();

   pcout << std::endl;
}

//------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------
int main (int argc, char *argv[])
{
   try
   {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      //dealllog.depth_console(0);

      unsigned int degree = 1;
      TestCase test_case = expo;
      Step12<2> dgmethod(degree, test_case);
      dgmethod.run ();
   }
   catch (std::exception &exc)
   {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
		          << std::endl;
      std::cerr << "Exception on processing: " << std::endl
		          << exc.what() << std::endl
		          << "Aborting!" << std::endl
		          << "----------------------------------------------------"
		          << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
   };

   return 0;
}
