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
   std::string  basis;
   int          mapping_degree;
   std::string  mapping;
   double       cfl;
   double       final_time;
   std::string  grid;
   unsigned int n_cells_x, n_cells_y;
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
            ProblemBase<dim>& problem,
            Quadrature<1>&    quadrature_1d);
   void run();

private:
   typedef parallel::distributed::Triangulation<dim> PTriangulation;
   typedef LinearAlgebra::distributed::Vector<double> PVector;

   void make_grid_and_dofs();
   const Mapping<dim, dim>& mapping() const;
   void initialize();
   void assemble_mass_matrix();
   void assemble_rhs();
   void compute_averages();
   void compute_dt();
   void apply_limiter();
   void apply_TVD_limiter();
   void update(const unsigned int rk_stage);
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
   double                      time, stage_time, dt;
   ProblemBase<dim>*           problem;
   ConditionalOStream          pcout;
   PTriangulation              triangulation;
   FESystem<dim>               fe;
   DoFHandler<dim>             dof_handler;
   AffineConstraints<double>   constraints;
   const Quadrature<dim>       cell_quadrature;
   const Quadrature<dim-1>     face_quadrature;
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
                        ProblemBase<dim>& problem,
                        Quadrature<1>&    quadrature_1d)
   :
   mpi_comm(MPI_COMM_WORLD),
   param(&param),
   problem(&problem),
   pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm) == 0)),
   triangulation(mpi_comm),
   fe(FE_DGQArbitraryNodes<dim>(quadrature_1d),nvar),
   dof_handler(triangulation),
   cell_quadrature(quadrature_1d),
   face_quadrature(quadrature_1d)
{
   AssertThrow(dim == 2, ExcIndexRange(dim, 0, 2));
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
   IndexSet locally_relevant_dofs;
   DoFTools::extract_locally_relevant_dofs(dof_handler,
                                           locally_relevant_dofs);

   // Solution and rhs variables
   solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
   solution_old.reinit(locally_owned_dofs, mpi_comm);
   rhs.reinit(solution);
   imm.reinit(solution_old);
   average.resize(counter, Vector<double>(nvar));

   // We dont have any constraints in DG.
   constraints.clear();
   constraints.close();

   pcout << "Mapping type   = " << param->mapping << std::endl;
   pcout << "Mapping degree = " << param->mapping_degree << std::endl;

   // check support point order. We assume that the order of cell_quadrature
   // points is same as the order of lagrange basis points. This allows us to 
   // directly get solution at quadrature points without using get_function_values.
   for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
   {
      auto ind_i = fe.system_to_component_index(i).second;
      auto q_point = cell_quadrature.point(ind_i);
      auto value = fe.shape_value(i, q_point);
      AssertThrow(fabs(value-1.0) < 1.0e-13, 
                  ExcMessage("Support point order assumption wrong"));
   }
}

//------------------------------------------------------------------------------
// Return mapping type based on selected type
//------------------------------------------------------------------------------
template <int dim>
const Mapping<dim, dim>& DGSystem<dim>::mapping() const
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
// With Legendre basis, mass matrix is diagonal, we only store diagonal part.
// Invert it and store
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::assemble_mass_matrix()
{
   pcout << "Constructing mass matrix ...\n";

   FEValues<dim> fe_values(mapping(), fe, cell_quadrature,
                           update_values | update_JxW_values);
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = cell_quadrature.size();
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
//------------------------------------------------------------------------------
template <int dim>
void
DGSystem<dim>::initialize()
{
   pcout << "Interpolating initial condition ...\n";

   FEValues<dim> fe_values(mapping(), fe, cell_quadrature,
                           update_quadrature_points);
   const unsigned int   n_q_points    = cell_quadrature.size();
   std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);

   for(auto & cell : dof_handler.active_cell_iterators())
   if(cell->is_locally_owned())
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);

      for(unsigned int q = 0; q < n_q_points; ++q)
      {
         Vector<double> initial_value(nvar);
         problem->initial_value(fe_values.quadrature_point(q),
                                initial_value);
         for(unsigned int i = 0; i < nvar; ++i)
         {
            auto idx = fe.component_to_system_index(i, q);
            solution(dof_indices[idx]) = initial_value[i];
         }
      }
   }

   solution.compress(VectorOperation::insert);
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
   auto &dof_indices = copy_data.local_dof_indices;

   for(unsigned int i=0; i<dofs_per_cell; ++i)
   {
      auto comp_i = fe_values.get_fe().system_to_component_index(i).first;
      auto indx_i = fe_values.get_fe().system_to_component_index(i).second;
      solution_values[indx_i][comp_i] = solution(dof_indices[i]);
   }

   for (unsigned int q = 0; q < n_q_points; ++q)
   {
      ndarray<double,nvar,dim> flux;
      PDE::physical_flux(solution_values[q],
                         fe_values.quadrature_point(q),
                         flux);
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
      Vector<double> num_flux(nvar);
      PDE::numerical_flux(param->flux_type, 
                          left_state[q], 
                          right_state[q], 
                          q_points[q], 
                          fe_face_values.normal(q),
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
      Vector<double> num_flux(nvar);
      PDE::boundary_flux(left_state[q], //todo
                         right_state[q],
                         q_points[q],
                         fe_face_values.normal_vector(q),
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

   ScratchData<dim> scratch_data(mapping(),
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
   FEValues<dim> fe_values(mapping(), fe, cell_quadrature,
                           update_JxW_values);
   std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
   const unsigned int n_q_points = cell_quadrature.size();

   for(auto & cell : dof_handler.active_cell_iterators())
   if(cell->is_locally_owned() || cell->is_ghost())
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);
      const auto c = cell->user_index();
      average[c] = 0.0;
      double cell_measure = 0.0;

      for(unsigned int q = 0; q < n_q_points; ++q)
      {
         cell_measure += fe_values.JxW(q);
         for(unsigned int i = 0; i < nvar; ++i)
         {
            auto idx = fe.component_to_system_index(i,q);
            average[c][i] += solution(dof_indices[idx]) * fe_values.JxW(q);
         }
      }

      average[c] /= cell_measure;
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
   AssertThrow(false, ExcNotImplemented());
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
      Tensor<1,dim> jac;
      PDE::max_speed(average[c], cell->center(), jac);
      double dtcell = cell->minimum_vertex_distance() / (jac.norm() + 1.0e-20);
      dt = std::min(dt, dtcell);
   }

   dt *= param->cfl;
   dt = Utilities::MPI::min(dt, mpi_comm);
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
   data_out.build_patches(mapping(), param->degree,
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
         solution.update_ghost_values();
         compute_averages();
         apply_limiter();
      }

      time += dt;
      ++iter;
      if(iter % param->output_step == 0) output_results(time);
      pcout << "Iter = " << iter 
                << " dt = " << dt
                << " time = " << time << std::endl;
   }
   output_results(time);
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
   prm.declare_entry("output step", "10", Patterns::Integer(0),
                     "Frequency to save solution");
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

   param.n_refine = ph.get_integer("initial refine");
   param.output_step = ph.get_integer("output step");

   param.cfl = ph.get_double("cfl");
   if(param.cfl == 0.0) param.cfl = 0.95 / (2 * param.degree + 1);
   if(param.cfl < 0.0) param.cfl = abs(param.cfl) / (2 * param.degree + 1);

   double final_time = ph.get_double("final time");
   if(final_time > 0.0)
      param.final_time = final_time;
   param.Mlim = ph.get_double("tvb parameter");

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
}
