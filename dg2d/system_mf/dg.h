//------------------------------------------------------------------------------
// Solves system of PDE of the form
//    u_t + div(f(u,x,t)) = s(u,x,t)
//------------------------------------------------------------------------------
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_interface_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/table.h>
#include <deal.II/base/timer.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/distributed/tria.h>


#include <fstream>
#include <iostream>

#include "pde.h"
#include "models/problem_base.h"

#define sign(a)   (((a) > 0.0) ? 1 : -1)

using namespace dealii;

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const unsigned int n_rk_stages = 3;
const double a_rk[3] = {0.0, 3.0 / 4.0, 1.0 / 3.0};
const double b_rk[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};

// Numerical flux functions
enum class LimiterType {none, tvd};

//------------------------------------------------------------------------------
// Scheme parameters
//------------------------------------------------------------------------------
struct Parameter
{
   int          degree;
   std::string  basis;
   int          mapping_degree;
   std::string  mapping;
   double       cfl;
   double       final_time;
   std::string  grid;
   unsigned int n_cells_x, n_cells_y, n_cells_z;
   unsigned int n_refine;
   unsigned int output_step;
   unsigned int output_number;
   double       output_interval;
   LimiterType  limiter_type;
   double       Mlim;
   FluxType     flux_type;
};

//------------------------------------------------------------------------------
// minmod of three numbers
//------------------------------------------------------------------------------
double
minmod(const double a, const double b, const double c, const double Mh2 = 0.0)
{
   double aa = std::fabs(a);
   if(aa < Mh2) return a;

   int sa = sign(a);
   int sb = sign(b);
   int sc = sign(c);

   double result;

   if(sa != sb || sb != sc)
   {
      result = 0.0;
   }
   else
   {
      result  = sa * std::min(aa, std::min(std::fabs(b), std::fabs(c)));
   }

   return result;
}

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim, int degree>
class DGSystem
{
public:
   DGSystem(Parameter&        param,
            ProblemBase<dim>& problem,
            Quadrature<1>&    quadrature_1d);
   void run();

private:
   typedef parallel::distributed::Triangulation<dim> PTriangulation;
   typedef LinearAlgebra::distributed::Vector<double> PVector;

   void make_grid_and_dofs();
   const Mapping<dim, dim> &mapping() const;
   void initialize();
   void assemble_mass_matrix();
   void assemble_rhs();
   void compute_averages();
   void compute_dt();
   void apply_limiter();
   void apply_TVD_limiter();
   void update(const unsigned int rk_stage);
   bool call_output();
   void output_results(const double time) const;

   const MPI_Comm                       mpi_comm;
   Parameter*                           param;
   double                               time, stage_time, dt, next_output_time;
   unsigned int                         time_step;
   ProblemBase<dim>*                    problem;
   ConditionalOStream                   pcout;
   mutable TimerOutput                  computing_timer;
   PTriangulation                       triangulation;
   const FESystem<dim>                  fe;
   DoFHandler<dim>                      dof_handler;
   AffineConstraints<double>            constraints;
   const Quadrature<dim>                cell_quadrature;
   const Quadrature<dim-1>              face_quadrature;
   PVector                              solution;
   PVector                              solution_old;
   PVector                              rhs;
   PVector                              imm;
   std::vector<Tensor<1, nvar, double>> average;
   MatrixFree<dim, double>              matrix_free;
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim, int degree>
DGSystem<dim, degree>::DGSystem(Parameter &param,
                                ProblemBase<dim> &problem,
                                Quadrature<1> &quadrature_1d)
   :
   mpi_comm(MPI_COMM_WORLD),
   param(&param),
   problem(&problem),
   pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm) == 0)),
   computing_timer(mpi_comm,
                     pcout,
                     TimerOutput::summary,
                     TimerOutput::wall_times),
   triangulation(mpi_comm),
   fe(FE_DGQArbitraryNodes<dim>(quadrature_1d), nvar),
   dof_handler(triangulation),
   cell_quadrature(quadrature_1d),
   face_quadrature(quadrature_1d)
{
   time = 0.0;
   time_step = 0;
   next_output_time = param.output_interval;
}

//------------------------------------------------------------------------------
// Make grid and allocate memory for solution variables
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::make_grid_and_dofs()
{
   TimerOutput::Scope t(computing_timer, "make grid and dofs");
   pcout << "Making initial grid ...\n";
   if(param->grid == "user")
   {
      pcout << "   User specified code for grid generation ...\n";
      problem->make_grid(triangulation);
   }
   else if(dim == 2 && param->grid == "box")
   {
      pcout << "   Making grid using subdivided_hyper_rectangle ...\n";
      pcout << "      Grid size = " << param->n_cells_x << " x "
            << param->n_cells_y << "\n";
      const Point<dim> p1(problem->get_xmin(), problem->get_ymin());
      const Point<dim> p2(problem->get_xmax(), problem->get_ymax());
      std::vector<unsigned int> ncells({param->n_cells_x,param->n_cells_y});
      GridGenerator::subdivided_hyper_rectangle(triangulation, ncells,
                                                p1, p2, true);
   }
   else if(dim == 3 && param->grid == "box")
   {
      pcout << "   Making grid using subdivided_hyper_rectangle ...\n";
      pcout << "      Grid size = " << param->n_cells_x << " x "
            << param->n_cells_y << " x "
            << param->n_cells_z << "\n";
      const Point<dim> p1(problem->get_xmin(),
                          problem->get_ymin(),
                          problem->get_zmin());
      const Point<dim> p2(problem->get_xmax(),
                          problem->get_ymax(),
                          problem->get_zmax());
      std::vector<unsigned int> ncells({param->n_cells_x,
                                        param->n_cells_y,
                                        param->n_cells_z});
      GridGenerator::subdivided_hyper_rectangle(triangulation, ncells,
                                                p1, p2, true);
   }
   else
   {
      pcout << "Reading gmsh grid from file " << param->grid << std::endl;
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream gfile(param->grid);
      AssertThrow(gfile.is_open(), ExcMessage("Grid file not found"));
      grid_in.read_msh(gfile);
   }

   if(problem->get_periodic())
   {
      typedef typename PTriangulation::cell_iterator Iter;
      std::vector<GridTools::PeriodicFacePair<Iter>> periodicity_vector;
      if(problem->get_periodic_x())
      {
         pcout << "   Applying periodic in x\n";
         GridTools::collect_periodic_faces(triangulation,
                                          0,
                                          1,
                                          0,
                                          periodicity_vector);
      }
      if(problem->get_periodic_y())
      {
         pcout << "   Applying periodic in y\n";
         GridTools::collect_periodic_faces(triangulation,
                                          2,
                                          3,
                                          1,
                                          periodicity_vector);
      }
      if(dim == 3 && problem->get_periodic_z())
      {
         pcout << "   Applying periodic in z\n";
         GridTools::collect_periodic_faces(triangulation,
                                          4,
                                          5,
                                          2,
                                          periodicity_vector);
      }
      triangulation.add_periodicity(periodicity_vector);
   }

   // Attach any manifold
   pcout << "   Setting manifolds\n";
   problem->set_manifolds(triangulation);

   // User specified transformation
   pcout << "   Transforming grid\n";
   problem->transform_grid(triangulation);

   if(param->n_refine > 0)
   {
      pcout << "   Refining initial grid\n";
      triangulation.refine_global(param->n_refine);
   }

   // Set an index for each cell, used to access cell averages
   unsigned int counter = 0;
   for(auto & cell : triangulation.active_cell_iterators())
      if(cell->is_locally_owned() || cell->is_ghost())
         cell->set_user_index(counter++);

   pcout << "   Number of active cells: "
         << triangulation.n_global_active_cells()
         << std::endl
         << "   Total number of cells: "
         << triangulation.n_cells()
         << std::endl;

    dof_handler.distribute_dofs(fe);

    pcout << "   Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;

    // We dont have any constraints in DG.
    constraints.clear();
    constraints.close();
    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_quadrature_points |
                                           update_JxW_values |
                                           update_values;
    additional_data.mapping_update_flags_inner_faces =
                                                     update_quadrature_points |
                                                     update_JxW_values |
                                                     update_values |
                                                     update_normal_vectors;
    additional_data.mapping_update_flags_boundary_faces =
                                                     update_quadrature_points |
                                                     update_JxW_values |
                                                     update_values |
                                                     update_normal_vectors;

    matrix_free.reinit(mapping(), dof_handler, constraints, cell_quadrature,
                       additional_data);
    average.resize(counter);

    matrix_free.initialize_dof_vector(solution);
    solution_old.reinit(solution);
    rhs.reinit(solution);
    imm.reinit(solution);

    pcout << "Mapping type   = " << param->mapping << std::endl;
    pcout << "Mapping degree = " << param->mapping_degree << std::endl;

    // TODO: This check is probably not needed in matrixfree code.
    // check support point order. We assume that the order of cell_quadrature
    // points is same as the order of lagrange basis points. This allows us to
    // directly get solution at quadrature points without using
    // get_function_values.
    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    {
      auto ind_i   = fe.system_to_component_index(i).second;
      auto q_point = cell_quadrature.point(ind_i);
      auto value   = fe.shape_value(i, q_point);
      AssertThrow(fabs(value - 1.0) < 1.0e-13,
                  ExcMessage("Support point order assumption wrong"));
    }
}

//------------------------------------------------------------------------------
// Return mapping type based on selected type
//------------------------------------------------------------------------------
template <int dim, int degree>
const Mapping<dim, dim>& DGSystem<dim,degree>::mapping() const
{
   if (param->mapping == "q")
   {
      static MappingQ<dim> m(param->mapping_degree);
      return m;
   }
   else if (param->mapping == "cartesian")
   {
      static MappingCartesian<dim> m;
      return m;
   }
   else
   {
      AssertThrow(false, ExcMessage("Requested mapping type is unknown"));
      static MappingQ1<dim> m;
      return m;
   }
}

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// With Lagrange basis at collocation points, mass matrix is diagonal,
// with entries equal to JxW at the corresponding quadrature point.
// Invert it and store.
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::assemble_mass_matrix()
{
   TimerOutput::Scope t(computing_timer, "mass matrix");
   pcout << "Constructing mass matrix ...\n";

   FEEvaluation<dim, degree, degree+1, nvar, double> phi(matrix_free);
   imm = 0.0;

   for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
   {
      phi.reinit(cell);
      const auto n_lanes = matrix_free.n_active_entries_per_cell_batch(cell);
      for (unsigned int lane = 0; lane < n_lanes; ++lane)
      {
         auto cell_it = matrix_free.get_cell_iterator(cell, lane);
         std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
         cell_it->get_dof_indices(dof_indices);
         for (unsigned int q = 0; q < phi.n_q_points; ++q)
         {
            const auto JxW = phi.JxW(q)[lane];
            for (unsigned int c = 0; c < nvar; ++c)
            {
               const auto idx = fe.component_to_system_index(c, q);
               imm(dof_indices[idx]) += JxW;
            }
         }
      }
   }

   imm.compress(VectorOperation::add);

   for (unsigned int i = 0; i < imm.locally_owned_size(); ++i)
     imm.local_element(i) = 1.0 / imm.local_element(i);
}

//------------------------------------------------------------------------------
// Set initial conditions
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::initialize()
{
   TimerOutput::Scope t(computing_timer, "initialize");
   pcout << "Interpolating initial condition ...\n";

   FEEvaluation<dim, degree, degree+1, nvar, double> phi(matrix_free);
   for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
   {
      phi.reinit(cell);
      const auto n_lanes = matrix_free.n_active_entries_per_cell_batch(cell);
      for (unsigned int lane = 0; lane < n_lanes; ++lane)
      {
         auto cell_it = matrix_free.get_cell_iterator(cell, lane);
         std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
         cell_it->get_dof_indices(dof_indices);

         for (unsigned int q = 0; q < phi.n_q_points; ++q)
         {
           Point<dim> p;
           for (unsigned int d = 0; d < dim; ++d)
             p[d] = phi.quadrature_point(q)[d][lane];
           Vector<double> initial_value(nvar);
           problem->initial_value(p, initial_value);
           for (unsigned int i = 0; i < nvar; ++i)
           {
             const auto idx = fe.component_to_system_index(i, q);
             solution(dof_indices[idx]) = initial_value[i];
           }
         }
      }
   }

   solution.compress(VectorOperation::insert);
}

//------------------------------------------------------------------------------
// DG operator for MatrixFree assembly
//------------------------------------------------------------------------------
template <int dim, int degree>
class DGOperator
{
public:
   DGOperator(const Parameter *param,
              const std::vector<Tensor<1, nvar, double>> *average,
              const double time)
         : param(param), average(average), time(time)
   {}

   //---------------------------------------------------------------------------
   // Assembly over cells
   //---------------------------------------------------------------------------
   void cell(const MatrixFree<dim, double> &mf,
             LinearAlgebra::distributed::Vector<double> &dst,
             const LinearAlgebra::distributed::Vector<double> &src,
             const std::pair<unsigned int, unsigned int> &range) const
   {
      FEEvaluation<dim, degree, degree+1, nvar, double> phi(mf);
      for (auto cell = range.first; cell < range.second; ++cell)
      {
         phi.reinit(cell);
         phi.gather_evaluate(src, EvaluationFlags::values);
         for (unsigned int q = 0; q < phi.n_q_points; ++q)
         {
            auto u = phi.get_value(q);
            FluxData<dim, VectorizedArray<double>> flux_data;
            flux_data.p = phi.quadrature_point(q);
            typename PDE::FluxMatrix<dim, VectorizedArray<double>> f;
            PDE::physical_flux<dim>(u, flux_data, f);
            phi.submit_gradient(f, q);
            #if defined(ADD_SOURCE)
            Tensor<1, nvar, VectorizedArray<double>> src_val;
            PDE::source(phi.quadrature_point(q), time, u, src_val);
            phi.submit_value(src_val, q);
            #endif
         }
         phi.integrate_scatter(EvaluationFlags::gradients
                              #if defined(ADD_SOURCE)
                              | EvaluationFlags::values
                              #endif
                              , dst);
      }
   }

   //---------------------------------------------------------------------------
   // Assembly over interior faces
   //---------------------------------------------------------------------------
   void face(const MatrixFree<dim, double> &mf,
             LinearAlgebra::distributed::Vector<double> &dst,
             const LinearAlgebra::distributed::Vector<double> &src,
             const std::pair<unsigned int, unsigned int> &range) const
   {
      FEFaceEvaluation<dim, degree, degree+1, nvar, double>
         phi_m(mf, true), phi_p(mf, false);
      for (auto face = range.first; face < range.second; ++face)
      {
         phi_m.reinit(face);
         phi_p.reinit(face);
         phi_m.gather_evaluate(src, EvaluationFlags::values);
         phi_p.gather_evaluate(src, EvaluationFlags::values);

         FluxData<dim, VectorizedArray<double>> flux_data;
         const auto n_lanes = mf.n_active_entries_per_face_batch(face);
         for (unsigned int lane = 0; lane < n_lanes; ++lane)
         {
            const auto c_l = mf.get_face_iterator(face, lane, true).first->user_index();
            const auto c_r = mf.get_face_iterator(face, lane, false).first->user_index();
            for (unsigned int c = 0; c < nvar; ++c)
            {
               flux_data.ul[c][lane] = (*average)[c_l][c];
               flux_data.ur[c][lane] = (*average)[c_r][c];
            }
         }

         for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
         {
            const auto ul = phi_m.get_value(q);
            const auto ur = phi_p.get_value(q);
            const auto n = phi_m.normal_vector(q);

            flux_data.p = phi_m.quadrature_point(q);

            typename PDE::NormalFlux<dim, VectorizedArray<double>> num_flux;
            PDE::numerical_flux(param->flux_type, ul, ur, n, flux_data, num_flux);

            phi_m.submit_value(-num_flux, q);
            phi_p.submit_value(num_flux, q);
         }
         phi_m.integrate_scatter(EvaluationFlags::values, dst);
         phi_p.integrate_scatter(EvaluationFlags::values, dst);
      }
   }

   //---------------------------------------------------------------------------
   // Assembly over boundary faces
   //---------------------------------------------------------------------------
   void boundary(const MatrixFree<dim, double> &mf,
                 LinearAlgebra::distributed::Vector<double> &dst,
                 const LinearAlgebra::distributed::Vector<double> &src,
                 const std::pair<unsigned int, unsigned int> &range) const
   {
      FEFaceEvaluation<dim, degree, degree+1, nvar, double> phi_m(mf, true);
      for (auto face = range.first; face < range.second; ++face)
      {
         phi_m.reinit(face);
         phi_m.gather_evaluate(src, EvaluationFlags::values);

         FluxData<dim, VectorizedArray<double>> flux_data;
         const auto n_lanes = mf.n_active_entries_per_face_batch(face);
         Tensor<1, nvar, VectorizedArray<double>> avg_l;
         for (unsigned int lane = 0; lane < n_lanes; ++lane)
         {
            const auto cell_user_index =
               mf.get_face_iterator(face, lane, true).first->user_index();
            for (unsigned int c = 0; c < nvar; ++c)
            {
               flux_data.ur[c][lane] = (*average)[cell_user_index][c];
               flux_data.ul[c][lane] = (*average)[cell_user_index][c];
            }
         }

         for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
         {
            const auto ul = phi_m.get_value(q);
            const auto ur = ul;
            const auto n  = phi_m.normal_vector(q);

            flux_data.p = phi_m.quadrature_point(q);

            typename PDE::NormalFlux<dim, VectorizedArray<double>> num_flux;
            PDE::boundary_flux(ul, ur, n, flux_data, num_flux);

            phi_m.submit_value(-num_flux, q);
         }
         phi_m.integrate_scatter(EvaluationFlags::values, dst);
      }
   }

private:
   const Parameter *param;
   const std::vector<Tensor<1, nvar, double>> *average;
   const double time;
}; // class DGOperator

//------------------------------------------------------------------------------
// Assemble system rhs
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::assemble_rhs()
{
   TimerOutput::Scope t(computing_timer, "assemble rhs");
   solution.update_ghost_values();
   rhs = 0.0;
   DGOperator<dim, degree> op(param, &average, stage_time);

   matrix_free.loop(&DGOperator<dim, degree>::cell,
                    &DGOperator<dim, degree>::face,
                    &DGOperator<dim, degree>::boundary,
                    &op,
                    rhs,
                    solution,
                    false,
                    MatrixFree<dim>::DataAccessOnFaces::values);

   rhs.compress(VectorOperation::add);
   rhs.scale(imm);
}

//------------------------------------------------------------------------------
// Compute cell average values
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::compute_averages()
{
   TimerOutput::Scope t(computing_timer, "compute averages");
   FEEvaluation<dim, degree, degree + 1, nvar, double> phi(matrix_free);

   const auto total =
       matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
   for (unsigned int cell = 0; cell < total; ++cell)
   {
      phi.reinit(cell);
      phi.gather_evaluate(solution, EvaluationFlags::values);

      Tensor<1, nvar, VectorizedArray<double>> cell_integral;
      VectorizedArray<double> cell_volume = 0.0;
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
         const auto val = phi.get_value(q);
         const auto JxW = phi.JxW(q);
         cell_volume += JxW;
         for (unsigned int i = 0; i < nvar; ++i)
            cell_integral[i] += val[i] * JxW;
      }

      for (unsigned int i = 0; i < nvar; ++i)
         cell_integral[i] /= cell_volume;

      const auto n_lanes = matrix_free.n_active_entries_per_cell_batch(cell);
      for (unsigned int lane = 0; lane < n_lanes; ++lane)
      {

         const auto cell_it = matrix_free.get_cell_iterator(cell, lane);
         const auto idx = cell_it->user_index();
         for (unsigned int i = 0; i < nvar; ++i)
            average[idx][i] = cell_integral[i][lane];
      }
  }
}

//------------------------------------------------------------------------------
// Apply TVD limiter
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::apply_TVD_limiter()
{
   if (degree == 0) return;
   AssertThrow(false, ExcNotImplemented());
}

//------------------------------------------------------------------------------
// Apply limiter
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::apply_limiter()
{
   TimerOutput::Scope t(computing_timer, "limiter");
   if (degree == 0) return;
   if (param->limiter_type == LimiterType::none) return;
   apply_TVD_limiter();
}

//------------------------------------------------------------------------------
// Compute time step from cfl condition
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::compute_dt()
{
   TimerOutput::Scope t(computing_timer, "compute dt");
   dt = 1.0e20;

   for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
   {
      unsigned int n_lanes = matrix_free.n_active_entries_per_cell_batch(cell);

      Tensor<1, nvar, VectorizedArray<double>> avg;
      Point<dim, VectorizedArray<double>> center;
      VectorizedArray<double> h;

      for (unsigned int lane = 0; lane < n_lanes; ++lane)
      {
         const auto cell_it = matrix_free.get_cell_iterator(cell, lane);
         const auto c = cell_it->user_index();
         const auto cell_center = cell_it->center();
         for (unsigned int i = 0; i < nvar; ++i)
            avg[i][lane] = average[c][i];
         for (unsigned int d = 0; d < dim; ++d)
            center[d][lane] = cell_center[d];
         h[lane] = cell_it->minimum_vertex_distance();
      }

      Tensor<1, dim, VectorizedArray<double>> jac;
      PDE::max_speed(avg, center, jac);

      VectorizedArray<double> dtcell = h / (jac.norm() + 1.0e-20);

      for (unsigned int lane = 0; lane < n_lanes; ++lane)
         dt = std::min(dt, dtcell[lane]);
   }

    dt *= param->cfl;
    dt  = Utilities::MPI::min(dt, mpi_comm);

   if (time + dt > param->final_time)
   {
      dt = param->final_time - time;
   }
   else if (param->output_interval > 0)
   {
      if (time + dt > next_output_time)
         dt = next_output_time - time;
   }
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::update(const unsigned int rk_stage)
{
   TimerOutput::Scope t(computing_timer, "update");

   // solution = solution + dt * rhs
   solution.add(dt, rhs);

   // solution = b_rk * solution + a_rk * solution_old
   solution.sadd(b_rk[rk_stage], a_rk[rk_stage], solution_old);

   stage_time = a_rk[rk_stage] * time + b_rk[rk_stage] * (stage_time + dt);
}

//-----------------------------------------------------------------------------
// Decide if solution needs to be saved
//-----------------------------------------------------------------------------
template <int dim, int degree>
bool
DGSystem<dim,degree>::call_output()
{
   // Save initial condition
   if (time_step == 0)
      return true;

   // Save final solution
   if (fabs(time - param->final_time) < 1.0e-13)
      return true;

   if (param->output_step > 0)
      if (time_step % param->output_step == 0)
         return true;

   if (param->output_interval > 0)
      if (fabs(time - next_output_time) < 1.0e-13)
      {
         next_output_time += param->output_interval;
         next_output_time = std::min(next_output_time, param->final_time);
         return true;
      }

   return false;
}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::output_results(const double time) const
{
   TimerOutput::Scope t(computing_timer, "output results");
   static unsigned int counter = 0;
   static std::vector<XDMFEntry> xdmf_entries;
   std::string mesh_filename = "mesh.h5";
   std::string solution_filename = ("vars-" +
                                   Utilities::int_to_string(counter, 4) +
                                   ".h5");
   bool write_mesh_file = (counter == 0) ? true : false;

   DataOut<dim> data_out;
   PDE::Postprocessor<dim> postprocessor;
   data_out.add_data_vector(dof_handler, solution, postprocessor);
   data_out.build_patches(mapping(), degree,
                          DataOut<dim>::curved_inner_cells);

   DataOutBase::DataOutFilter data_filter(DataOutBase::DataOutFilterFlags(true, true));
  // Filter the data and store it in data_filter
  data_out.write_filtered_data(data_filter);
  // Write the filtered data to HDF5
  data_out.write_hdf5_parallel(data_filter,
                               write_mesh_file,
                               mesh_filename,
                               solution_filename,
                               mpi_comm);
  // Create an XDMF entry detailing the HDF5 file
  XDMFEntry new_xdmf_entry = data_out.create_xdmf_entry(data_filter,
                                                        mesh_filename,
                                                        solution_filename,
                                                        time,
                                                        mpi_comm);
  // Add the XDMF entry to the list
  xdmf_entries.push_back(new_xdmf_entry);
  // Create an XDMF file from all stored entries
  data_out.write_xdmf_file(xdmf_entries, "solution.xdmf", mpi_comm);

  pcout << "Wrote " << solution_filename << " at t = " << time << "\n";
   ++counter;
}

//------------------------------------------------------------------------------
// Start solving the problem
//------------------------------------------------------------------------------
template <int dim, int degree>
void
DGSystem<dim,degree>::run()
{
  pcout << "Solving " << PDE::name << " for " << problem->get_name() << "\n";
  pcout << "Number of threads   = " << MultithreadInfo::n_threads() << "\n";
  pcout << "Number of ranks     = " << Utilities::MPI::n_mpi_processes(mpi_comm)
                                    << "\n";
  pcout << "Vectorization width = " << VectorizedArray<double>::size() << "\n";

  if (Utilities::MPI::this_mpi_process(mpi_comm) == 0)
    PDE::print_info();
  make_grid_and_dofs();
  assemble_mass_matrix();
  initialize();
  solution.update_ghost_values();
  compute_averages();
  output_results(0.0);

  while (time < param->final_time)
  {
    solution_old = solution;
    stage_time = time;
    compute_dt();

    for (unsigned int rk = 0; rk < n_rk_stages; ++rk)
    {
      assemble_rhs();
      update(rk);
      solution.update_ghost_values();
      compute_averages();
      apply_limiter();
    }

    time += dt, ++time_step;
    pcout << "Iter = " << time_step
          << " dt = " << dt
          << " time = " << time
          << std::endl;
    if (call_output())
      output_results(time);
  }

  computing_timer.print_summary();
}

//------------------------------------------------------------------------------
// Declare input parameters
//------------------------------------------------------------------------------
void
declare_parameters(ParameterHandler& prm)
{
   prm.declare_entry("degree", "0", Patterns::Integer(0),
                     "Polynomial degree");
   prm.declare_entry("basis", "gl", Patterns::Selection("gl|gll"),
                     "Specify basis: gl or gll");
   prm.declare_entry("mapping", "q,1", Patterns::Anything(),
                     "Specify mapping: cartesian or q or q,1 or q,2 etc.");
   prm.declare_entry("grid", "0", Patterns::Anything(),
                     "Specify grid: 100,100 or user or foo.msh");
   prm.declare_entry("initial refine", "0", Patterns::Integer(0),
                     "Number of grid refinements");
   prm.declare_entry("output step", "0", Patterns::Integer(0),
                     "Iteration frequency to save solution");
   prm.declare_entry("output number", "0", Patterns::Integer(0),
                     "How many time to save solution");
   prm.declare_entry("output interval", "0.0", Patterns::Double(0),
                     "Time frequency to save solution");
   prm.declare_entry("cfl", "0.0", Patterns::Double(),
                     "CFL number");
   prm.declare_entry("final time", "0.0", Patterns::Double(0),
                     "Final time");
   prm.declare_entry("limiter", "none",
                     Patterns::Selection("none|tvd"),
                     "Limiter");
   prm.declare_entry("numflux", "central",
                     Patterns::Anything(),
                     "Numerical flux");
   prm.declare_entry("tvb parameter", "0.0", Patterns::Double(0),
                     "TVB parameter");
}

//------------------------------------------------------------------------------
void
parse_parameters(const ParameterHandler& ph, Parameter& param)
{
   param.degree = ph.get_integer("degree");
   param.basis = ph.get("basis");

   {
      auto value = ph.get("mapping");
      if(value == "cartesian")
      {
         param.mapping = "cartesian";
         param.mapping_degree = 0; // not needed
      }
      else if(value == "q")
      {
         param.mapping = "q";
         param.mapping_degree = param.degree;
      }
      else
      {
         auto values = Utilities::split_string_list(value, ",");
         if(values[0] == "q")
         {
            param.mapping = "q";
            param.mapping_degree = Utilities::string_to_int(values[1]);
            AssertThrow(param.mapping_degree >= 1,
                        ExcMessage("Need mapping degree >= 1"));
            AssertThrow(param.mapping_degree <= param.degree,
                        ExcMessage("Need mapping degree <= degree"));
         }
         else
         {
            AssertThrow(false, ExcMessage("Unknown mapping"));
         }
      }
   }

   auto grid = ph.get("grid");
   AssertThrow(grid != "0", ExcMessage("Grid is not specified."));
   auto grid_size = Utilities::split_string_list(grid, ",");
   if(grid_size.size() == 2)
   {
      param.grid = "box";
      param.n_cells_x = Utilities::string_to_int(grid_size[0]);
      param.n_cells_y = Utilities::string_to_int(grid_size[1]);
   }
   else if(grid == "user")
   {
      param.grid = "user";
   }
   else
   {
      param.grid = grid; // gmsh grid file name
   }

   double final_time = ph.get_double("final time");
   if(final_time > 0.0)
      param.final_time = final_time;

   param.n_refine = ph.get_integer("initial refine");

   param.output_step = ph.get_integer("output step");
   param.output_number = ph.get_integer("output number");
   param.output_interval = ph.get_double("output interval");
   if (param.output_step > 0)
   {
      param.output_number = 0;
      param.output_interval = 0.0;
   }
   else if (param.output_number > 0)
   {
      param.output_step = 0;
      param.output_interval = param.final_time / (param.output_number - 1);
   }
   else if (param.output_interval > 0.0)
   {
      param.output_step = 0;
      param.output_number = 0;
   }
   else
   {
      AssertThrow(false, ExcMessage("No output settings given"));
   }

   param.cfl = ph.get_double("cfl");
   if(param.cfl == 0.0) param.cfl = 0.95 / (2 * param.degree + 1);
   if(param.cfl < 0.0) param.cfl = abs(param.cfl) / (2 * param.degree + 1);

   {
      std::string value = ph.get("numflux");
      auto search = FluxTypeList.find(value);
      if (search != FluxTypeList.end())
         param.flux_type = search->second;
      else
      {
         std::cout << "Available num fluxes\n";
         for (const auto &v : FluxTypeList)
            std::cout << "   * " << v.first << std::endl;
         AssertThrow(false, ExcMessage("Unknown flux type"));
      }
   }

   {
      std::string value = ph.get("limiter");
      if (value == "none") param.limiter_type = LimiterType::none;
      else if (value == "tvd") param.limiter_type = LimiterType::tvd;
      else AssertThrow(false, ExcMessage("Unknown limiter"));
   }

   param.Mlim = ph.get_double("tvb parameter");
}

