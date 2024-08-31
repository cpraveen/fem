#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_interface_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <fstream>
#include <iostream>

#include "pde.h"

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
   double       xmin, xmax;
   std::string  basis;
   int          degree;
   double       cfl;
   double       final_time;
   unsigned int n_cells;
   unsigned int n_refine;
   unsigned int output_step;
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
                 Quadrature<1>&       quadrature_1d,
                 const Function<dim>& initial_condition,
                 const Function<dim>& exact_solution);
   void run();

private:
   void solve();
   void make_grid_and_dofs(unsigned int step);
   void initialize();
   void assemble_mass_matrix();
   void assemble_rhs();
   void compute_averages();
   void compute_dt();
   void apply_limiter();
   void apply_TVD_limiter();
   void mark_troubled_cells();
   void update(const unsigned int rk_stage);
   void output_results(const double time) const;
   void process_solution(unsigned int step);

   Parameter*           param;
   double               time, stage_time, dt;
   double               xmin, xmax, dx;
   unsigned int         n_rk_stages;

   const Quadrature<dim>       cell_quadrature;
   const Function<dim>*        initial_condition;
   const Function<dim>*        exact_solution;

   Triangulation<dim>          triangulation;
   FE_DGQArbitraryNodes<dim>   fe;
   DoFHandler<dim>             dof_handler;
   Vector<double>              solution;
   Vector<double>              solution_old;
   Vector<double>              rhs;
   Vector<double>              imm;
   Vector<double>              average;
   std::vector<bool>           troubled_cell;
   ConvergenceTable            convergence_table;
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
ScalarProblem<dim>::ScalarProblem(Parameter&           param,
                                  Quadrature<1>&       quadrature_1d,
                                  const Function<dim>& initial_condition,
                                  const Function<dim>& exact_solution)
   :
   param(&param),
   cell_quadrature(quadrature_1d),
   initial_condition(&initial_condition),
   exact_solution(&exact_solution),
   fe(quadrature_1d),
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
ScalarProblem<dim>::make_grid_and_dofs(unsigned int step)
{
   if(step == 0)
   {
      std::cout << "Making initial grid ...\n";
      GridGenerator::subdivided_hyper_cube(triangulation, param->n_cells, xmin, xmax);
      typedef typename Triangulation<dim>::cell_iterator Iter;
      std::vector<GridTools::PeriodicFacePair<Iter>> periodicity_vector;
      GridTools::collect_periodic_faces(triangulation,
                                        0,
                                        1,
                                        0,
                                        periodicity_vector);
      triangulation.add_periodicity(periodicity_vector);
   }
   else
   {
      triangulation.refine_global(1);
   }
   dx = (xmax - xmin) / triangulation.n_active_cells();

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

   average.reinit(triangulation.n_active_cells());
   troubled_cell.resize(triangulation.n_active_cells());
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
   std::cout << "  Quadrature using " << fe.degree + 1 << " points\n";

   FEValues<dim> fe_values(fe, cell_quadrature,
                           update_JxW_values);
   std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);

   // Cell iterator
   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);
      for(unsigned int i = 0; i < fe.dofs_per_cell; ++i)
         imm[dof_indices[i]] = 1.0 / fe_values.JxW(i);
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
   std::cout << "Interpolate initial condition ...\n";

   FEValues<dim> fe_values(fe, cell_quadrature,
                           update_quadrature_points);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = cell_quadrature.size();
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);
      for(unsigned int q = 0; q < n_q_points; ++q)
      {
         double initial_value = initial_condition->value(fe_values.quadrature_point(q));
         solution(dof_indices[q]) = initial_value;
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
                           update_gradients |
                           update_quadrature_points |
                           update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = cell_quadrature.size();
   std::vector<double>  solution_values(n_q_points);

   // for getting neighbor cell solution using trapezoidal rule
   Vector<double>       cell_rhs(dofs_per_cell);
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   FEInterfaceValues<dim> fe_face_values(fe,
                                         Quadrature<dim-1> (Point<dim>(0.0)),
                                         update_values);

   rhs = 0.0;

   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);

      // Compute conserved variables at quadrature points
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
         solution_values[i] = solution(dof_indices[i]);

      // Flux integral over cell
      cell_rhs  = 0.0;
      for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
         double flux = physical_flux(solution_values[q_point],
                                     fe_values.quadrature_point(q_point));
         for(unsigned int i = 0; i < dofs_per_cell; ++i)
         {
            cell_rhs(i) += (fe_values.shape_grad(i, q_point)[0] *
                            flux *
                            fe_values.JxW(q_point));
         }
      }

      // Add cell residual to rhs
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
      {
         auto ig = dof_indices[i];
         rhs(ig) += cell_rhs(i);
      }

      // Add face residual to rhs
      // assemble flux at right face only since we have periodic bc
      auto ncell = cell->neighbor_or_periodic_neighbor(1);
      fe_face_values.reinit(cell, 1, 0,
                            ncell, cell->neighbor_of_neighbor(1), 0);
      const unsigned int n_face_dofs = fe_face_values.n_current_interface_dofs();
      const auto& face_dof_indices = fe_face_values.get_interface_dof_indices();
      std::vector<double> left_state(1), right_state(1);
      fe_face_values.get_fe_face_values(0).get_function_values(solution, left_state);
      fe_face_values.get_fe_face_values(1).get_function_values(solution, right_state);
      double num_flux;
      numerical_flux(param->flux_type, left_state[0], right_state[0], 
                     cell->face(1)->center(), num_flux);
      for(unsigned int i = 0; i < n_face_dofs; ++i)
      {
         auto ig = face_dof_indices[i];
         rhs(ig) += -num_flux *
                    fe_face_values.jump_in_shape_values(i, 0);
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
   FEValues<dim> fe_values(fe, cell_quadrature,
                           update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);
      double integral = 0.0;
      double measure = 0.0;
      for(unsigned int q = 0; q < dofs_per_cell; ++q)
      {
         integral += solution(dof_indices[q]) * fe_values.JxW(q);
         measure += fe_values.JxW(q);
      }
      average[cell->user_index()] = integral / measure;
   }
}

//------------------------------------------------------------------------------
// Use TVB limiter to identify troubled cells
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::mark_troubled_cells()
{
   AssertThrow(false, ExcNotImplemented());
   const double EPS = 1.0e-14;
   QTrapezoid<dim>  quadrature_formula;
   FEValues<dim> fe_values(fe, quadrature_formula, update_values);
   std::vector<double> face_values(2), limited_face_values(2);

   unsigned int n_troubled_cells = 0;
   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      fe_values.get_function_values(solution, face_values);

      auto c = cell->user_index();
      auto cl = cell->neighbor_or_periodic_neighbor(0)->user_index();
      auto cr = cell->neighbor_or_periodic_neighbor(1)->user_index();

      double db = average[c]  - average[cl];
      double df = average[cr] - average[c];

      double DF = face_values[1] - average[c];
      double DB = average[c] - face_values[0];

      double dl = minmod(DB, db, df);
      double dr = minmod(DF, db, df);

      limited_face_values[0] = average[c] - dl;
      limited_face_values[1] = average[c] + dr;
      bool c0 = std::fabs(limited_face_values[0] - face_values[0])
                > EPS * std::fabs(face_values[0]);
      bool c1 = std::fabs(limited_face_values[1] - face_values[1])
                > EPS * std::fabs(face_values[1]);

      troubled_cell[c] = false;
      if(c0 || c1)
      {
         troubled_cell[c] = true;
         ++n_troubled_cells;
      }
   }
   //std::cout << "No of troubled cells = " << n_troubled_cells << std::endl;
}

//------------------------------------------------------------------------------
// Apply TVD limiter
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::apply_TVD_limiter()
{
   if(fe.degree == 0) return;

   const double EPS = 1.0e-6;
   QTrapezoid<dim> quadrature_formula;
   FEValues<dim> fe_values(fe, quadrature_formula, update_values);
   std::vector<double> face_values(2), limited_face_values(2);

   const double Mh2 = param->Mlim * dx * dx;
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
   const auto xi = fe.get_unit_support_points();

   for(auto & cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit(cell);
      fe_values.get_function_values(solution, face_values);

      auto c  = cell->user_index();
      auto cl = cell->neighbor_or_periodic_neighbor(0)->user_index();
      auto cr = cell->neighbor_or_periodic_neighbor(1)->user_index();

      double db = average[c] - average[cl];
      double df = average[cr] - average[c];

      double DF = face_values[1] - average[c];
      double DB = average[c] - face_values[0];

      double dl = minmod(DB, db, df, Mh2);
      double dr = minmod(DF, db, df, Mh2);

      limited_face_values[0] = average[c] - dl;
      limited_face_values[1] = average[c] + dr;
      bool c0 = std::fabs(limited_face_values[0] - face_values[0]) 
                          > EPS * std::fabs(face_values[0]);
      bool c1 = std::fabs(limited_face_values[1] - face_values[1]) 
                          > EPS * std::fabs(face_values[1]);

      if(c0 || c1)
      {
         double Dx = 0.5 * ( dl + dr); // limited average slope
         cell->get_dof_indices(dof_indices);
         for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
            solution(dof_indices[i]) = average[c] + (xi[i][0] - 0.5) * Dx;
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
   if(fe.degree == 0 || param->limiter_type == LimiterType::none) return;
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
                      / (max_speed(average[c], cell->center()) + 1.0e-20);
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

   if(fe.degree <= 1)
      data_out.build_patches(1);
   else
      data_out.build_patches(2 * fe.degree);

   std::string filename = "sol_" + Utilities::int_to_string(counter) + ".gpl";
   std::cout << "Outout at t = " << time << "  " << filename << std::endl;

   std::ofstream output(filename);
   data_out.write_gnuplot(output);

   std::ofstream fo;
   filename = "avg_" + Utilities::int_to_string(counter) + ".gpl";
   fo.open(filename);

   for(auto & cell : triangulation.active_cell_iterators())
   {
      Point<dim> x = cell->center();
      auto c = cell->user_index();
      int tc = (troubled_cell[c] ? 1 : 0);
      fo << x(0) << " " << average[c] << "  " << tc << std::endl;
   }

   fo.close();

   ++counter;
}

//------------------------------------------------------------------------------
// Start solving the problem
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::solve()
{
   std::cout << "Solving 1-D scalar problem ...\n";

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
// Compute error norms
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::process_solution(unsigned int step)
{
   // compute error in solution
   Vector<double> difference_per_cell(triangulation.n_active_cells());
   VectorTools::integrate_difference(dof_handler,
                                     solution,
                                     *exact_solution,
                                     difference_per_cell,
                                     QGauss<dim>(2 * fe.degree + 1),
                                     VectorTools::L2_norm);
   const double L2_error = difference_per_cell.l2_norm();

   // compute error in gradient
   VectorTools::integrate_difference(dof_handler,
                                     solution,
                                     *exact_solution,
                                     difference_per_cell,
                                     QGauss<dim>(2 * fe.degree + 1),
                                     VectorTools::H1_norm);
   const double H1_error = difference_per_cell.l2_norm();

   const unsigned int n_active_cells = triangulation.n_active_cells();
   const unsigned int n_dofs = dof_handler.n_dofs();
   convergence_table.add_value("step", step);
   convergence_table.add_value("cells", n_active_cells);
   convergence_table.add_value("dofs", n_dofs);
   convergence_table.add_value("L2", L2_error);
   convergence_table.add_value("H1", H1_error);

   convergence_table.set_scientific("L2", true);
   convergence_table.set_scientific("H1", true);

   std::cout << std::endl;
}

//------------------------------------------------------------------------------
// Run with several refinements
// Compute error and convergence rate
//------------------------------------------------------------------------------
template <int dim>
void
ScalarProblem<dim>::run()
{
   for(unsigned int step = 0; step < param->n_refine; ++step)
   {
      make_grid_and_dofs(step);
      solve();
      process_solution(step);
   }

   convergence_table.set_precision("L2", 3);
   convergence_table.set_scientific("L2", true);
   convergence_table.evaluate_convergence_rates("L2",
                    ConvergenceTable::reduction_rate_log2);

   convergence_table.set_precision("H1", 3);
   convergence_table.set_scientific("H1", true);
   convergence_table.evaluate_convergence_rates("H1",
                    ConvergenceTable::reduction_rate_log2);

   convergence_table.set_tex_caption("cells", "\\# cells");
   convergence_table.set_tex_caption("dofs", "\\# dofs");
   convergence_table.set_tex_caption("L2", "$L^2$-error");
   convergence_table.set_tex_caption("H1", "$H^1$-error");

   convergence_table.set_tex_format("cells", "r");
   convergence_table.set_tex_format("dofs", "r");

   std::cout << std::endl;
   convergence_table.write_text(std::cout);

   std::ofstream error_table_file("error.tex");
   convergence_table.write_tex(error_table_file);
}

//------------------------------------------------------------------------------
// Declare input parameters
//------------------------------------------------------------------------------
void
declare_parameters(ParameterHandler& prm)
{
   prm.declare_entry("basis", "gl",
                     Patterns::Selection("gl|gll"),
                     "GL vs GLL points");
   prm.declare_entry("degree", "0", Patterns::Integer(0),
                     "Polynomial degree");
   prm.declare_entry("ncells", "100", Patterns::Integer(2),
                     "Number of elements");
   prm.declare_entry("nrefine", "1", Patterns::Integer(1, 10),
                     "Number of grid refinements");
   prm.declare_entry("output step", "10", Patterns::Integer(0),
                     "Frequency to save solution");
   prm.declare_entry("test case", "sine",
                     Patterns::Anything(),
                     "Test case");
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
   param.basis = ph.get("basis");
   param.degree = ph.get_integer("degree");
   param.n_cells = ph.get_integer("ncells");
   param.n_refine = ph.get_integer("nrefine");
   param.output_step = ph.get_integer("output step");
   param.cfl = ph.get_double("cfl");
   if(param.cfl == 0.0) param.cfl = 0.98 / (2 * param.degree + 1);
   if(param.cfl < 0.0) param.cfl = param.cfl / (2 * param.degree + 1);

   param.final_time = ph.get_double("final time");
   param.Mlim = ph.get_double("tvb parameter");

   {
      std::string value = ph.get("numflux");
      auto search = FluxTypeList.find(value);
      if(search != FluxTypeList.end())
         param.flux_type = search->second;
      else
      {
         std::cout << "Available num fluxes\n";
         for(const auto& v : FluxTypeList)
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
