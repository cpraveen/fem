#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_interface_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <fstream>
#include <iostream>

#include "pde.h"
#include "problem.h"

#define dsign(a)   (((a) > 0.0) ? 1 : -1)

using namespace dealii;

// Coefficients for 3-stage SSP RK scheme of Shu-Osher
const double a_rk[3] = {0.0, 3.0 / 4.0, 1.0 / 3.0};
const double b_rk[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};

// Numerical flux functions
enum class LimiterType {none, tvd};

//------------------------------------------------------------------------------
// Scheme parameters
//------------------------------------------------------------------------------
struct Parameter
{
   int           degree;
   double        cfl;
   double        final_time;
   unsigned int  n_cells;
   unsigned int  output_step;
   LimiterType   limiter_type;
   double        Mlim;
   PDE::FluxType flux_type;
};

//------------------------------------------------------------------------------
// minmod of three numbers
//------------------------------------------------------------------------------
double
minmod(const double a, const double b, const double c, const double Mh2 = 0.0)
{
   double aa = std::fabs(a);
   if(aa < Mh2) return a;

   int sa = dsign(a);
   int sb = dsign(b);
   int sc = dsign(c);

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
template <int dim>
class ScalarProblem
{
public:
   ScalarProblem(Parameter&           param,
                 Quadrature<1>&       cell_quadrature);
   void run();

private:
   void make_grid_and_dofs();
   void initialize();
   void assemble_mass_matrix();
   void assemble_rhs();
   void compute_averages();
   void compute_dt();
   void apply_limiter();
   void apply_TVD_limiter();
   void update(const unsigned int rk_stage);
   void output_results(const double time) const;
   void process_solution(unsigned int step);

   Parameter*           param;
   double               time, stage_time, dt;
   double               dx;
   unsigned int         n_rk_stages;

   const Quadrature<dim>       cell_quadrature;
   Triangulation<dim>          triangulation;
   FESystem<dim>               fe;
   DoFHandler<dim>             dof_handler;
   Vector<double>              solution;
   Vector<double>              solution_old;
   Vector<double>              rhs;
   Vector<double>              imm;
   std::vector<Vector<double>> average;
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
ScalarProblem<dim>::ScalarProblem(Parameter&           param,
                                  Quadrature<1>&       cell_quadrature)
   :
   param(&param),
   cell_quadrature(cell_quadrature),
   fe(FE_DGP<dim>(param.degree),nvar),
   dof_handler(triangulation)
{
   AssertThrow(dim == 1, ExcIndexRange(dim, 0, 1));

   n_rk_stages = 3;
}

//------------------------------------------------------------------------------
// Make grid and allocate memory for solution variables
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::make_grid_and_dofs()
{
   std::cout << "Making grid ...\n";
   GridGenerator::subdivided_hyper_cube(triangulation, param->n_cells, 
                                        Problem::xmin, Problem::xmax);
   if(Problem::periodic)
   {
      std::cout << "Applying periodicity to grid\n";
      typedef typename Triangulation<dim>::cell_iterator Iter;
      std::vector<GridTools::PeriodicFacePair<Iter>> periodicity_vector;
      GridTools::collect_periodic_faces(triangulation,
                                          0,
                                          1,
                                          0,
                                          periodicity_vector);
      triangulation.add_periodicity(periodicity_vector);
   }
   dx = (Problem::xmax - Problem::xmin) / triangulation.n_active_cells();

   unsigned int counter = 0;
   for(auto & cell : triangulation.active_cell_iterators())
      cell->set_user_index(counter++);

   std::cout << "   Number of active cells: "
             << triangulation.n_active_cells()
             << std::endl
             << "   Total number of cells: "
             << triangulation.n_cells()
             << std::endl;

   dof_handler.distribute_dofs(fe);

   std::cout << "   Number of degrees of freedom: "
             << dof_handler.n_dofs()
             << std::endl;

   // Solution variables
   solution.reinit(dof_handler.n_dofs());
   solution_old.reinit(dof_handler.n_dofs());
   rhs.reinit(dof_handler.n_dofs());
   imm.reinit(dof_handler.n_dofs());

   average.resize(triangulation.n_active_cells(), Vector<double>(nvar));
}

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// With Legendre basis, mass matrix is diagonal, we only store diagonal part.
// Invert it and store
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::assemble_mass_matrix()
{
   std::cout << "Constructing mass matrix ...\n";

   FEValues<dim> fe_values(fe, cell_quadrature,
                           update_values | update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = cell_quadrature.size();
   Vector<double>   cell_matrix(dofs_per_cell);
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   imm = 0.0;

   // Cell iterator
   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      cell_matrix = 0.0;

      for(unsigned int q = 0; q < n_q_points; ++q)
         for(unsigned int i = 0; i < dofs_per_cell; ++i)
         {
            cell_matrix(i) += fe_values.shape_value(i, q) *
                              fe_values.shape_value(i, q) *
                              fe_values.JxW(q);
         }

      cell->get_dof_indices(dof_indices);
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
         imm[dof_indices[i]] = 1.0 / cell_matrix(i);
   }
}

//------------------------------------------------------------------------------
// Set initial conditions
// L2 projection of initial condition onto dofs
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::initialize()
{
   std::cout << "Projecting initial condition ...\n";

   QGauss<dim>  quadrature_formula(2 * param->degree + 1);
   FEValues<dim> fe_values(fe, quadrature_formula,
                           update_values   |
                           update_quadrature_points |
                           update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   Vector<double>       cell_rhs(dofs_per_cell);
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      cell_rhs  = 0.0;

      // integral over cell
      for(unsigned int q = 0; q < n_q_points; ++q)
      {
         // Get primitive variable at quadrature point
         Vector<double> initial_value(nvar);
         Problem::initial_value(fe_values.quadrature_point(q),
                                initial_value);
         for(unsigned int i = 0; i < dofs_per_cell; ++i)
         {
            const auto c = fe.system_to_component_index(i).first;
            cell_rhs(i) += fe_values.shape_value(i, q) *
                           initial_value[c] *
                           fe_values.JxW(q);
         }
      }

      // Multiply by inverse mass matrix and add to rhs
      cell->get_dof_indices(dof_indices);
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
      {
         auto ig = dof_indices[i];
         solution(ig) = imm(ig) * cell_rhs(i);
      }
   }

}

//------------------------------------------------------------------------------
// Assemble system rhs
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::assemble_rhs()
{
   FEValues<dim> fe_values(fe, cell_quadrature,
                           update_values   | 
                           update_gradients |
                           update_quadrature_points |
                           update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = cell_quadrature.size();
   std::vector<Vector<double>>  solution_values(n_q_points,Vector<double>(nvar));

   // for getting neighbor cell solution using trapezoidal rule
   Vector<double>       cell_rhs(dofs_per_cell);
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
   std::vector<types::global_dof_index> dof_indices_nbr(dofs_per_cell);

   FEFaceValues<dim> fe_face_values0(fe,
                                     Quadrature<dim-1> (Point<dim>(0.0)),
                                     update_values);
   FEFaceValues<dim> fe_face_values1(fe,
                                     Quadrature<dim-1> (Point<dim>(0.0)),
                                     update_values);
   std::vector<Vector<double>>  left_state(1,Vector<double>(nvar));
   std::vector<Vector<double>>  right_state(1,Vector<double>(nvar));

   rhs = 0.0;

   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);

      // Compute conserved variables at quadrature points
      fe_values.get_function_values(solution,  solution_values);

      // Flux integral over cell
      cell_rhs  = 0.0;
      for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
         Vector<double> flux(nvar);
         PDE::physical_flux(solution_values[q_point],
                            fe_values.quadrature_point(q_point),
                            flux);
         for(unsigned int i = 0; i < dofs_per_cell; ++i)
         {
            auto c = fe.system_to_component_index(i).first;
            cell_rhs(i) += (fe_values.shape_grad(i, q_point)[0] *
                            flux[c] *
                            fe_values.JxW(q_point));
         }
      }

      // Add cell residual to rhs
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
         rhs(dof_indices[i]) += cell_rhs(i);

      // Add face residual to rhs
      // assemble flux at right face only since we have periodic bc
      auto ncell = cell->neighbor_or_periodic_neighbor(1);
      fe_face_values0.reinit(cell, 1);
      fe_face_values1.reinit(ncell, cell->neighbor_of_neighbor(1));
      fe_face_values0.get_function_values(solution, left_state);
      fe_face_values1.get_function_values(solution, right_state);
      Vector<double> num_flux(nvar);
      PDE::numerical_flux(param->flux_type, left_state[0], right_state[0],
                          cell->face(1)->center(), num_flux);
      // Add to left cell
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
      {
         auto c = fe.system_to_component_index(i).first;
         rhs(dof_indices[i]) -= num_flux[c] * 
                                fe_face_values0.shape_value(i, 0);
      }
      // Add to right cell
      ncell->get_dof_indices(dof_indices_nbr);
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
      {
         auto c = fe.system_to_component_index(i).first;
         rhs(dof_indices_nbr[i]) += num_flux[c] * 
                                    fe_face_values1.shape_value(i, 0);
      }
   }

   // Multiply by inverse mass matrix
   rhs.scale(imm);
}

//------------------------------------------------------------------------------
// Compute cell average values
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::compute_averages()
{
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
   const unsigned int dofs_per_component = param->degree + 1;

   for(auto & cell : dof_handler.active_cell_iterators())
   {
      cell->get_dof_indices(dof_indices);
      const auto c = cell->user_index();
      unsigned int ii = 0;
      for (unsigned int i = 0; i < nvar; ++i, ii += dofs_per_component)
         average[c][i] = solution(dof_indices[ii]);
   }
}

//------------------------------------------------------------------------------
// Apply TVD limiter
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::apply_TVD_limiter()
{
   if(param->degree == 0) return;

   const double Mh2 = param->Mlim * dx * dx;
   const double sqrt_3 = std::sqrt(3.0);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
   Vector<double> solution_new(dofs_per_cell);
   const unsigned int dofs_per_component = param->degree + 1;
   Vector<double> db(nvar), df(nvar), Dx(nvar);
   Vector<double> db1(nvar), df1(nvar), Dx1(nvar);
   FullMatrix<double> R(nvar,nvar), L(nvar,nvar);
   Vector<double> Dx1_new(nvar), Dx_new(nvar);

   for(auto & cell : dof_handler.active_cell_iterators())
   {
      auto c  = cell->user_index();
      auto cl = cell->neighbor_or_periodic_neighbor(0)->user_index();
      auto cr = cell->neighbor_or_periodic_neighbor(1)->user_index();
      cell->get_dof_indices(dof_indices);

      unsigned int idx = 1;
      for(unsigned int comp=0; comp<nvar; ++comp, idx+=dofs_per_component)
      {
         db[comp] = average[c][comp]  - average[cl][comp];
         df[comp] = average[cr][comp] - average[c][comp];
         Dx[comp] = solution(dof_indices[idx]);
      }

      PDE::char_mat(average[c], cell->center(), R, L);
      L.vmult(db1, db);
      L.vmult(df1, df);
      L.vmult(Dx1, Dx);

      bool tolimit = false;
      for(unsigned int comp=0; comp<nvar; ++comp)
      {
         Dx1_new[comp] = minmod(sqrt_3 * Dx1[comp], db1[comp], df1[comp], Mh2) / sqrt_3;
         if(fabs(Dx1[comp] - Dx1_new[comp]) > 1.0e-6 * fabs(Dx1[comp]))
            tolimit = true;
      }

      if(tolimit)
      {
         R.vmult(Dx_new, Dx1_new);
         solution_new = 0.0;
         idx = 0;
         for(unsigned int comp=0; comp<nvar; ++comp, idx+=dofs_per_component)
         {
            solution(dof_indices[idx]) = average[c][comp];
            solution(dof_indices[idx+1]) = Dx_new[comp];
         }
      }
   }

}

//------------------------------------------------------------------------------
// Apply TVD limiter
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::apply_limiter()
{
   if(param->degree == 0 || param->limiter_type == LimiterType::none) return;
   apply_TVD_limiter();
}

//------------------------------------------------------------------------------
// Compute time step from cfl condition
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::compute_dt()
{
   dt = 1.0e20;

   for(auto &cell : dof_handler.active_cell_iterators())
   {
      auto c = cell->user_index();
      double dtcell = cell->measure() 
                      / (PDE::max_speed(average[c], cell->center()) + 1.0e-20);
      dt = std::min(dt, dtcell);
   }

   dt *= param->cfl;
}

//------------------------------------------------------------------------------
// Update solution by one stage of RK
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::update(const unsigned int rk_stage)
{
   // Update conserved variables
   for(unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
   {
      solution(i)  = a_rk[rk_stage] * solution_old(i) +
                     b_rk[rk_stage] * (solution(i) + dt* rhs(i));
   }

   stage_time = a_rk[rk_stage] * time + b_rk[rk_stage] * (stage_time + dt);
}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::output_results(const double time) const
{
   static unsigned int counter = 0;

   DataOut<dim> data_out;
   data_out.attach_dof_handler(dof_handler);
   data_out.add_data_vector(solution, "solution");

   if(param->degree <= 1)
      data_out.build_patches(1);
   else
      data_out.build_patches(2 * param->degree);

   std::string filename = "sol_" + Utilities::int_to_string(counter) + ".gpl";
   std::cout << "Output at t = " << time << "  " << filename << std::endl;

   std::ofstream output(filename);
   data_out.write_gnuplot(output);

   std::ofstream fo;
   filename = "avg_" + Utilities::int_to_string(counter) + ".gpl";
   fo.open(filename);

   for(auto & cell : triangulation.active_cell_iterators())
   {
      Point<dim> x = cell->center();
      auto c = cell->user_index();
      fo << x(0) << " ";
      for(unsigned int i=0; i<nvar; ++i)
         fo << average[c][i] << " ";
      fo << std::endl;
   }

   fo.close();

   ++counter;
}

//------------------------------------------------------------------------------
// Start solving the problem
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::run()
{
   std::cout << "Solving 1-D scalar problem ...\n";

   make_grid_and_dofs();
   assemble_mass_matrix();
   initialize();
   compute_averages();
   output_results(0.0);

   time = 0.0;
   unsigned int iter = 0;

   while(time < param->final_time)
   {
      solution_old  = solution;
      stage_time = time;

      compute_dt();
      if(time + dt > param->final_time) dt = param->final_time - time;

      for(unsigned int rk = 0; rk < n_rk_stages; ++rk)
      {
         assemble_rhs();
         update(rk);
         compute_averages();
         apply_limiter();
      }

      time += dt;
      ++iter;
      if(iter % param->output_step == 0) output_results(time);
   }
   output_results(time);
   std::cout << "Iter = " << iter << " time = " << time
             << std::endl;
}

//------------------------------------------------------------------------------
// Declare input parameters
//------------------------------------------------------------------------------
void
declare_parameters(ParameterHandler& prm)
{
   prm.declare_entry("degree", "0", Patterns::Integer(0),
                     "Polynomial degree");
   prm.declare_entry("ncells", "100", Patterns::Integer(2),
                     "Number of elements");
   prm.declare_entry("output step", "10", Patterns::Integer(0),
                     "Frequency to save solution");
   prm.declare_entry("cfl", "0.0", Patterns::Double(0, 1.0),
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
   param.n_cells = ph.get_integer("ncells");
   param.output_step = ph.get_integer("output step");
   param.cfl = ph.get_double("cfl");
   if(param.cfl == 0.0) param.cfl = 0.98 / (2 * param.degree + 1);
   if(param.cfl < 0.0) param.cfl = param.cfl / (2 * param.degree + 1);

   param.final_time = ph.get_double("final time");
   param.Mlim = ph.get_double("tvb parameter");

   {
      std::string value = ph.get("numflux");
      auto search = PDE::FluxTypeList.find(value);
      if(search != PDE::FluxTypeList.end())
         param.flux_type = search->second;
      else
      {
         std::cout << "Available num fluxes\n";
         for(const auto& v : PDE::FluxTypeList)
            std::cout << v.first << std::endl;
         AssertThrow(false, ExcMessage("Unknown flux type"));
      }
   }

   {
      std::string value = ph.get("limiter");
      if (value == "none") param.limiter_type = LimiterType::none;
      else if (value == "tvd") param.limiter_type = LimiterType::tvd;
      else AssertThrow(false, ExcMessage("Unknown limiter"));
   }
}
