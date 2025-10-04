//------------------------------------------------------------------------------
// Solves system of PDE of the form
//    u_t + div(f(u,x)) = 0
//------------------------------------------------------------------------------
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_interface_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/distributed/tria.h>


#include <fstream>
#include <iostream>

#include "pde.h"
#include "../models/problem_base.h"

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
   double       cfl;
   double       final_time;
   std::string  grid;
   unsigned int n_cells_x, n_cells_y;
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
// Find cell size dx, dy for cartesian grid
//------------------------------------------------------------------------------
template <typename CIterator>
void cell_size(const CIterator &cell, double &dx, double &dy)
{
   const auto dr1 = cell->face(1)->center() - cell->face(0)->center();
   const auto dr2 = cell->face(3)->center() - cell->face(2)->center();
   dx = std::max(fabs(dr1[0]), fabs(dr2[0]));
   dy = std::max(fabs(dr1[1]), fabs(dr2[1]));
}

//------------------------------------------------------------------------------
template <int dim>
struct ScratchData
{
   ScratchData(const Mapping<dim> &mapping,
               const FiniteElement<dim> &fe,
               const Quadrature<dim> &cell_quadrature,
               const Quadrature<dim-1> &face_quadrature,
               const UpdateFlags update_flags = update_values |
                                                update_gradients |
                                                update_quadrature_points |
                                                update_JxW_values,
               const UpdateFlags interface_update_flags = update_values | 
                                                          update_quadrature_points |
                                                          update_JxW_values | 
                                                          update_normal_vectors)
       : 
       fe_values(mapping, fe, cell_quadrature, update_flags), 
       fe_interface_values(mapping,
                           fe,
                           face_quadrature,
                           interface_update_flags),
      solution_values(cell_quadrature.size(), Vector<double>(nvar)),
      left_state(face_quadrature.size(), Vector<double>(nvar)),
      right_state(face_quadrature.size(), Vector<double>(nvar))
   {
   }

   ScratchData(const ScratchData<dim> &scratch_data)
       : fe_values(scratch_data.fe_values.get_mapping(),
                   scratch_data.fe_values.get_fe(),
                   scratch_data.fe_values.get_quadrature(),
                   scratch_data.fe_values.get_update_flags()),
         fe_interface_values(scratch_data.fe_interface_values.get_mapping(),
                             scratch_data.fe_interface_values.get_fe(),
                             scratch_data.fe_interface_values.get_quadrature(),
                             scratch_data.fe_interface_values.get_update_flags()),
         solution_values(scratch_data.fe_values.get_quadrature().size(),
                         Vector<double>(nvar)),
         left_state(scratch_data.fe_interface_values.get_quadrature().size(),
                    Vector<double>(nvar)),
         right_state(scratch_data.fe_interface_values.get_quadrature().size(),
                     Vector<double>(nvar))
   {
   }

   FEValues<dim> fe_values;
   FEInterfaceValues<dim> fe_interface_values;
   std::vector<Vector<double>> solution_values;
   std::vector<Vector<double>> left_state;
   std::vector<Vector<double>> right_state;
};

//------------------------------------------------------------------------------
struct CopyDataFace
{
   std::vector<types::global_dof_index> joint_dof_indices;
   Vector<double> cell_rhs;
};

//------------------------------------------------------------------------------
struct CopyData
{
   Vector<double> cell_rhs;
   std::vector<types::global_dof_index> local_dof_indices;
   std::vector<CopyDataFace> face_data;

   template <class Iterator>
   void reinit(const Iterator &cell, unsigned int dofs_per_cell)
   {
      cell_rhs.reinit(dofs_per_cell);

      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
   }
};

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
template <int dim>
class DGSystem
{
public:
   DGSystem(Parameter&        param,
            ProblemBase<dim>& problem);
   void run();

private:
   typedef parallel::distributed::Triangulation<dim> PTriangulation;
   typedef LinearAlgebra::distributed::Vector<double> PVector;

   void make_grid_and_dofs();
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

   template <class Iterator>
   void cell_worker(const Iterator &cell,
                    ScratchData<dim> &scratch_data,
                    CopyData &copy_data);

   template <class Iterator>
   void boundary_worker(const Iterator &cell,
                        const unsigned int &f,
                        ScratchData<dim> &scratch_data,
                        CopyData &copy_data);

   template <class Iterator>
   void face_worker(const Iterator &cell,
                    const unsigned int &f,
                    const unsigned int &sf,
                    const Iterator &ncell,
                    const unsigned int &nf,
                    const unsigned int &nsf,
                    ScratchData<dim> &scratch_data,
                    CopyData &copy_data);

   const MPI_Comm              mpi_comm;
   Parameter*                  param;
   double                      time, stage_time, dt, next_output_time;
   unsigned int                time_step;
   ProblemBase<dim>*           problem;
   ConditionalOStream          pcout;
   PTriangulation              triangulation;
   FESystem<dim>               fe;
   DoFHandler<dim>             dof_handler;
   MappingCartesian<dim>       mapping;
   AffineConstraints<double>   constraints;
   PVector                     solution;
   PVector                     solution_old;
   PVector                     rhs;
   PVector                     imm;
   std::vector<Vector<double>> average;
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <int dim>
DGSystem<dim>::DGSystem(Parameter&        param,
                        ProblemBase<dim>& problem)
   :
   mpi_comm(MPI_COMM_WORLD),
   param(&param),
   problem(&problem),
   pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm) == 0)),
   triangulation(mpi_comm),
   fe(FE_DGP<dim>(param.degree),nvar),
   dof_handler(triangulation)
{
   AssertThrow(dim == 2, ExcIndexRange(dim, 0, 2));

   time = 0.0;
   time_step = 0;
   next_output_time = param.output_interval;
}

//------------------------------------------------------------------------------
// Make grid and allocate memory for solution variables
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::make_grid_and_dofs()
{
   pcout << "Making initial grid ...\n";
   if(param->grid == "user")
   {
      pcout << "   User specified code for grid generation ...\n";
      problem->make_grid(triangulation);
   }
   else if(param->grid == "box")
   {
      pcout << "   Making grid using subdivided_hyper_rectangle ...\n";
      pcout << "      Grid size = " << param->n_cells_x << " x " 
            << param->n_cells_y << "\n";
      const Point<dim> p1(problem->get_xmin(), problem->get_ymin());
      const Point<dim> p2(problem->get_xmax(), problem->get_ymax());
      std::vector<unsigned int> ncells2d({param->n_cells_x,param->n_cells_y});
      GridGenerator::subdivided_hyper_rectangle(triangulation, ncells2d,
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
      triangulation.add_periodicity(periodicity_vector);
   }

   // User specified transformation. NOTE: Cells must remain rectangles.
   pcout << "   Transforming grid\n";
   problem->transform_grid(triangulation);

   if(param->n_refine > 0)
   {
      pcout << "   Refining initial grid\n";
      triangulation.refine_global(param->n_refine);
   }

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

   const auto& locally_owned_dofs = dof_handler.locally_owned_dofs();
   IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

   // Solution and rhs variables
   solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
   solution_old.reinit(locally_owned_dofs, mpi_comm);
   rhs.reinit(solution);
   imm.reinit(solution_old);
   average.resize(counter, Vector<double>(nvar));

   // We dont have any constraints in DG.
   constraints.clear();
   constraints.close();
}

//------------------------------------------------------------------------------
// Assemble mass matrix for each cell
// With Legendre basis, mass matrix is diagonal, we only store diagonal part.
// Invert it and store
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::assemble_mass_matrix()
{
   pcout << "Constructing mass matrix ...\n";
   pcout << "  Quadrature using " << param->degree + 1 << " points\n";

   QGauss<dim>  quadrature_formula(param->degree + 1);
   FEValues<dim> fe_values(mapping, fe, quadrature_formula,
                           update_values | update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   Vector<double>   cell_matrix(dofs_per_cell);
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   imm = 0.0;

   // Cell iterator
   for(auto & cell : dof_handler.active_cell_iterators())
   if(cell->is_locally_owned())
   {
      fe_values.reinit(cell);
      cell_matrix = 0.0;

      for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
         for(unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_matrix(i) += fe_values.shape_value(i, q_point) *
                              fe_values.shape_value(i, q_point) *
                              fe_values.JxW(q_point);

      cell->get_dof_indices(dof_indices);
      imm.add(dof_indices, cell_matrix);
   }

   imm.compress(VectorOperation::add);

   // Invert mass matrix
   for (unsigned int i = 0; i < imm.locally_owned_size(); ++i)
   {
      imm.local_element(i) = 1.0 / imm.local_element(i);
   }
}

//------------------------------------------------------------------------------
// Set initial conditions
// L2 projection of initial condition onto dofs
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::initialize()
{
   pcout << "Projecting initial condition ...\n";

   QGauss<dim>  quadrature_formula(2 * param->degree + 1);
   FEValues<dim> fe_values(mapping, fe, quadrature_formula,
                           update_values   |
                           update_quadrature_points |
                           update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
   Vector<double>       cell_rhs(dofs_per_cell);
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   solution = 0.0;

   for(auto & cell : dof_handler.active_cell_iterators())
   if(cell->is_locally_owned())
   {
      fe_values.reinit(cell);
      cell_rhs  = 0.0;

      // integral over cell
      for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
         // Get primitive variable at quadrature point
         Vector<double> initial_value(nvar);
         problem->initial_value(fe_values.quadrature_point(q_point),
                                initial_value);
         for(unsigned int i = 0; i < dofs_per_cell; ++i)
         {
            auto c = fe.system_to_component_index(i).first;
            cell_rhs(i) += fe_values.shape_value(i, q_point) *
                           initial_value[c] *
                           fe_values.JxW(q_point);
         }
      }

      // Multiply by inverse mass matrix and add to rhs
      cell->get_dof_indices(dof_indices);
      solution.add(dof_indices, cell_rhs);
   }

   solution.compress(VectorOperation::add);
   solution.scale(imm);
}

//------------------------------------------------------------------------------
template <int dim>
template <class Iterator>
void DGSystem<dim>::cell_worker(const Iterator &cell,
                                ScratchData<dim> &scratch_data,
                                CopyData &copy_data)
{
   FEValues<dim> &fe_values = scratch_data.fe_values;
   fe_values.reinit(cell);

   const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell();
   const unsigned int n_q_points = fe_values.get_quadrature().size();

   copy_data.reinit(cell, dofs_per_cell);

   auto &cell_rhs = copy_data.cell_rhs;
   auto &solution_values = scratch_data.solution_values;
   fe_values.get_function_values(solution,  solution_values);

   for (unsigned int q = 0; q < n_q_points; ++q)
   {
      FluxData<dim> data;
      data.p = fe_values.quadrature_point(q);
      data.t = stage_time;
      ndarray<double,nvar,dim> flux;
      PDE::physical_flux(solution_values[q], data, flux);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
         const auto c = fe_values.get_fe().system_to_component_index(i).first;
         const auto& shape_grad = fe_values.shape_grad_component(i,q,c);
         double tmp = 0.0;
         for(unsigned int d=0; d<dim; ++d) tmp += shape_grad[d] * flux[c][d];
         cell_rhs(i) += tmp * fe_values.JxW(q);
      }

   }
}

//------------------------------------------------------------------------------
template <int dim>
template <class Iterator>
void DGSystem<dim>::face_worker(const Iterator &cell,
                                const unsigned int &f,
                                const unsigned int &sf,
                                const Iterator &ncell,
                                const unsigned int &nf,
                                const unsigned int &nsf,
                                ScratchData<dim> &scratch_data,
                                CopyData &copy_data)
{
   FEInterfaceValues<dim> &fe_face_values = scratch_data.fe_interface_values;
   fe_face_values.reinit(cell, f, sf, ncell, nf, nsf);

   const unsigned int n_cell_dofs = fe_face_values.get_fe().n_dofs_per_cell();
   const unsigned int n_face_dofs = fe_face_values.n_current_interface_dofs();
   const unsigned int n_q_points = fe_face_values.get_quadrature().size();
   const auto &q_points = fe_face_values.get_quadrature_points();

   auto &left_state = scratch_data.left_state;
   auto &right_state = scratch_data.right_state;
   fe_face_values.get_fe_face_values(0).get_function_values(solution, left_state);
   fe_face_values.get_fe_face_values(1).get_function_values(solution, right_state);

   copy_data.face_data.emplace_back();
   CopyDataFace &copy_data_face = copy_data.face_data.back();
   copy_data_face.joint_dof_indices = fe_face_values.get_interface_dof_indices();
   copy_data_face.cell_rhs.reinit(n_face_dofs);
   auto &cell_rhs = copy_data_face.cell_rhs;

   for(unsigned int q=0; q<n_q_points; ++q)
   {
      FluxData<dim> data;
      data.p = q_points[q];
      data.t = stage_time;
      data.ul = &average[cell->user_index()];
      data.ur = &average[ncell->user_index()];
      Vector<double> num_flux(nvar);
      PDE::numerical_flux(param->flux_type, 
                          left_state[q], 
                          right_state[q], 
                          fe_face_values.normal(q),
                          data,
                          num_flux);
      for (unsigned int i = 0; i < n_face_dofs; ++i)
      {
         unsigned int ii = (i < n_cell_dofs) ? i : i - n_cell_dofs;
         const auto c = fe_face_values.get_fe().system_to_component_index(ii).first;
         cell_rhs(i) -= num_flux[c] *
                        fe_face_values.jump_in_shape_values(i, q, c) *
                        fe_face_values.JxW(q);
      }
   }
}

//------------------------------------------------------------------------------
template <int dim>
template <class Iterator>
void DGSystem<dim>::boundary_worker(const Iterator &cell,
                                    const unsigned int &f,
                                    ScratchData<dim> &scratch_data,
                                    CopyData &copy_data)
{
   scratch_data.fe_interface_values.reinit(cell, f);
   const auto &fe_face_values 
      = scratch_data.fe_interface_values.get_fe_face_values(0);

   const unsigned int n_face_dofs = fe_face_values.get_fe().n_dofs_per_cell();
   const unsigned int n_q_points = fe_face_values.get_quadrature().size();
   const auto &q_points = fe_face_values.get_quadrature_points();

   auto &left_state = scratch_data.left_state;
   auto &right_state = scratch_data.right_state;
   fe_face_values.get_function_values(solution, left_state);
   auto &cell_rhs = copy_data.cell_rhs;

   for (unsigned int q = 0; q < n_q_points; ++q)
   {
      problem->boundary_value(cell->face(f)->boundary_id(),
                              q_points[q],
                              stage_time,
                              fe_face_values.normal_vector(q),
                              left_state[q],
                              right_state[q]);
      FluxData<dim> data;
      data.p = q_points[q];
      data.t = stage_time;
      data.ul = &average[cell->user_index()];
      data.ur = &average[cell->user_index()];
      Vector<double> num_flux(nvar);
      PDE::boundary_flux(left_state[q], //todo
                         right_state[q],
                         fe_face_values.normal_vector(q),
                         data,
                         num_flux);
      for (unsigned int i = 0; i < n_face_dofs; ++i)
      {
         const auto c = fe_face_values.get_fe().system_to_component_index(i).first;
         cell_rhs(i) -= num_flux[c] *
                        fe_face_values.shape_value_component(i, q, c) *
                        fe_face_values.JxW(q);
      }
   }
}

//------------------------------------------------------------------------------
// Assemble system rhs
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::assemble_rhs()
{
   using Iterator = typename DoFHandler<dim>::active_cell_iterator;

   auto cell_worker =
       [&](const Iterator &cell,
           ScratchData<dim> &scratch_data,
           CopyData &copy_data)
   {
      this->cell_worker(cell, scratch_data, copy_data);
   };

   auto face_worker =
       [&](const Iterator &cell,
           const unsigned int f,
           const unsigned int sf,
           const Iterator &ncell,
           const unsigned int nf,
           const unsigned int nsf,
           ScratchData<dim> &scratch_data,
           CopyData &copy_data)
   {
      this->face_worker(cell, f, sf, ncell, nf, nsf, scratch_data, copy_data);
   };

   auto boundary_worker =
       [&](const Iterator &cell,
           const unsigned int f,
           ScratchData<dim> &scratch_data,
           CopyData &copy_data)
   {
      this->boundary_worker(cell, f, scratch_data, copy_data);
   };

   auto copier = [&](const CopyData &cd)
   {
      this->constraints.distribute_local_to_global(cd.cell_rhs,
                                                   cd.local_dof_indices,
                                                   this->rhs);
      for (auto &cdf : cd.face_data)
      {
         this->constraints.distribute_local_to_global(cdf.cell_rhs,
                                                      cdf.joint_dof_indices,
                                                      this->rhs);
      }
   };

   const unsigned int n_gauss_points = param->degree + 1;
   const QGauss<dim> cell_quadrature(n_gauss_points);
   const QGauss<dim-1> face_quadrature(n_gauss_points);

   ScratchData<dim> scratch_data(mapping,
                                 fe,
                                 cell_quadrature,
                                 face_quadrature);

   const auto iterator_range =
        filter_iterators(dof_handler.active_cell_iterators(),
                         IteratorFilters::LocallyOwnedCell());

   rhs = 0.0;
   MeshWorker::mesh_loop(iterator_range,
                         cell_worker,
                         copier,
                         scratch_data,
                         CopyData(),
                         MeshWorker::assemble_own_cells |
                         MeshWorker::assemble_boundary_faces |
                         MeshWorker::assemble_own_interior_faces_once |
                         MeshWorker::assemble_ghost_faces_once,
                         boundary_worker,
                         face_worker);

   // Reduce over all MPI ranks
   rhs.compress(VectorOperation::add);

   // Multiply by inverse mass matrix
   rhs.scale(imm);
}

//------------------------------------------------------------------------------
// Compute cell average values
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::compute_averages()
{
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
   const unsigned int dofs_per_comp = ((param->degree + 1)* (param->degree + 2)) / 2;

   for(auto & cell : dof_handler.active_cell_iterators())
   if(cell->is_locally_owned() || cell->is_ghost())
   {
      cell->get_dof_indices(dof_indices);
      const auto c = cell->user_index();
      unsigned int j = 0;
      for(unsigned int i = 0; i < nvar; ++i, j+=dofs_per_comp)
         average[c][i] = solution(dof_indices[j]);
   }
}

//------------------------------------------------------------------------------
// Apply TVD limiter: 2d case only
// TODO: Make it work on locally refined grids
//------------------------------------------------------------------------------
template <>
void
DGSystem<2>::apply_TVD_limiter()
{
   if(param->degree == 0) return;

   const double sqrt_3 = std::sqrt(3.0);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
   const unsigned int degree = param->degree;
   const unsigned int dofs_per_comp = ((degree+1)*(degree+2))/2;
   Vector<double> dbx(nvar), dfx(nvar), Dx(nvar), Dx_new(nvar);
   Vector<double> dby(nvar), dfy(nvar), Dy(nvar), Dy_new(nvar);
   Vector<double> dbx1(nvar), dfx1(nvar), Dx1(nvar), Dx1_new(nvar);
   Vector<double> dby1(nvar), dfy1(nvar), Dy1(nvar), Dy1_new(nvar);
   FullMatrix<double> Rx(nvar,nvar), Lx(nvar,nvar), Ry(nvar,nvar), Ly(nvar,nvar);

   for(auto & cell : dof_handler.active_cell_iterators())
   if(cell->is_locally_owned())
   {
      double dx, dy;
      cell_size(cell, dx, dy);
      const double h = std::max(dx, dy);
      const double Mh2 = param->Mlim * h * h;
      const auto c  = cell->user_index();

      unsigned int cl, cr, cb, ct;

      // left cell
      if (cell->face(0)->at_boundary() && cell->has_periodic_neighbor(0) == false)
      {
         cl = c; // TODO: assuming neumann-like bc
      }
      else
      {
         cl = cell->neighbor_or_periodic_neighbor(0)->user_index();
      }

      // right cell
      if (cell->face(1)->at_boundary() && cell->has_periodic_neighbor(1) == false)
      {
         cr = c; // TODO: assuming neumann-like bc
      }
      else
      {
         cr = cell->neighbor_or_periodic_neighbor(1)->user_index();
      }

      // bottom cell
      if (cell->face(2)->at_boundary() && cell->has_periodic_neighbor(2) == false)
      {
         cb = c; // TODO: assuming neumann-like bc
      }
      else
      {
         cb = cell->neighbor_or_periodic_neighbor(2)->user_index();
      }

      // top cell
      if (cell->face(3)->at_boundary() && cell->has_periodic_neighbor(3) == false)
      {
         ct = c; // TODO: assuming neumann-like bc
      }
      else
      {
         ct = cell->neighbor_or_periodic_neighbor(3)->user_index();
      }

      cell->get_dof_indices(dof_indices);

      for(unsigned int i=0, j=0; i<nvar; ++i, j+=dofs_per_comp)
      {
         dbx[i] = average[c][i]  - average[cl][i];
         dfx[i] = average[cr][i] - average[c][i];
         Dx[i] = solution(dof_indices[j+1]);

         dby[i] = average[c][i]  - average[cb][i];
         dfy[i] = average[ct][i] - average[c][i];
         Dy[i] = solution(dof_indices[j+degree+1]);
      }

      // TODO: Transform to characteristic
      const auto drx = cell->face(1)->center() - cell->face(0)->center();
      const auto ex = drx / drx.norm();
      const auto dry = cell->face(3)->center() - cell->face(2)->center();
      const auto ey = dry / dry.norm();
      PDE::char_mat(average[c], cell->center(), ex, ey, Rx, Lx, Ry, Ly);
      Lx.vmult(dbx1, dbx);
      Lx.vmult(dfx1, dfx);
      Lx.vmult(Dx1,  Dx);
      Ly.vmult(dby1, dby);
      Ly.vmult(dfy1, dfy);
      Ly.vmult(Dy1,  Dy);

      bool tolimit = false;
      for(unsigned int i=0; i<nvar; ++i)
      {
         Dx1_new[i] = minmod(sqrt_3 * Dx1[i], dbx1[i], dfx1[i], Mh2) / sqrt_3;
         Dy1_new[i] = minmod(sqrt_3 * Dy1[i], dby1[i], dfy1[i], Mh2) / sqrt_3;
         if(fabs(Dx1[i] - Dx1_new[i]) > 1.0e-6 * fabs(Dx1[i]) || 
            fabs(Dy1[i] - Dy1_new[i]) > 1.0e-6 * fabs(Dy1[i]))
            tolimit = true;
      }

      if(tolimit)
      {
         Rx.vmult(Dx_new, Dx1_new);
         Ry.vmult(Dy_new, Dy1_new);
         for(unsigned int i = 0; i < dofs_per_cell; ++i)
            solution(dof_indices[i]) = 0;
         for(unsigned int i=0, j=0; i<nvar; ++i, j+=dofs_per_comp)
         {
            solution(dof_indices[j]) = average[c][i];
            solution(dof_indices[j+1]) = Dx_new[i];
            solution(dof_indices[j+degree+1]) = Dy_new[i];
         }
      }
   }
}

//------------------------------------------------------------------------------
// Apply TVD limiter
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::apply_limiter()
{
   if(param->degree == 0 || param->limiter_type == LimiterType::none) return;
   apply_TVD_limiter();
}

//------------------------------------------------------------------------------
// Compute time step from cfl condition
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::compute_dt()
{
   dt = 1.0e20;

   for(auto &cell : dof_handler.active_cell_iterators())
   if(cell->is_locally_owned())
   {
      auto c = cell->user_index();
      double dx, dy;
      cell_size(cell, dx, dy);
      Tensor<1,dim> jac;
      PDE::max_speed(average[c], cell->center(), jac);
      double dtcell = 1.0 / (fabs(jac[0])/dx + fabs(jac[1])/dy + 1.0e-20);
      dt = std::min(dt, dtcell);
   }

   dt *= param->cfl;
   dt = Utilities::MPI::min(dt, mpi_comm);

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
template <int dim>
void
DGSystem<dim>::update(const unsigned int rk_stage)
{
   // solution = solution + dt * rhs
   solution.add(dt, rhs);

   // solution = b_rk * solution + a_rk * solution_old
   solution.sadd(b_rk[rk_stage], a_rk[rk_stage], solution_old);

   stage_time = a_rk[rk_stage] * time + b_rk[rk_stage] * (stage_time + dt);
}

//-----------------------------------------------------------------------------
// Decide if solution needs to be saved
//-----------------------------------------------------------------------------
template <int dim>
bool DGSystem<dim>::call_output()
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
template <int dim>
void
DGSystem<dim>::output_results(const double time) const
{
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
   data_out.build_patches(mapping, param->degree);

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
template <int dim>
void
DGSystem<dim>::run()
{
   pcout << "Solving " << PDE::name << " for " << problem->get_name() << "\n";
   pcout << "Number of threads = " << MultithreadInfo::n_threads() << "\n";

   if (Utilities::MPI::this_mpi_process(mpi_comm) == 0)
      PDE::print_info();
   make_grid_and_dofs();
   assemble_mass_matrix();
   initialize();
   solution.update_ghost_values();
   compute_averages();
   output_results(0.0);

   while(time < param->final_time)
   {
      solution_old  = solution;
      stage_time = time;
      compute_dt();

      for(unsigned int rk = 0; rk < n_rk_stages; ++rk)
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
            << " time = " << time << std::endl;
      if(call_output()) output_results(time);
   }
}

//------------------------------------------------------------------------------
// Declare input parameters
//------------------------------------------------------------------------------
void
declare_parameters(ParameterHandler& prm)
{
   prm.declare_entry("degree", "0", Patterns::Integer(0),
                     "Polynomial degree");
   prm.declare_entry("basis", "legendre", Patterns::Anything(),
                     "Specify basis: NOT USED, always Legendre");
   prm.declare_entry("mapping", "cartesian", Patterns::Anything(),
                     "Specify mapping: NOT USED, always cartesian");
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
