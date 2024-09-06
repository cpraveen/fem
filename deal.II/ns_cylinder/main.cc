/*
 Solution of incompressible Navier-Stokes equations using Taylor-Hood
 
 BDF1 in first time step; BDF2 from second time step onwards
 Convective term is linearized by extrapolation
 Second order accurate in time
 
 Author: Praveen. C, 
         TIFR-CAM, Bangalore
         http://praveen.tifrbng.res.in
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parameter_handler.h>


#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>

using namespace dealii;

enum runMode { norun, steady, unsteady };
runMode run_mode;
bool restart;
std::string parameter_file;

//------------------------------------------------------------------------------------
// Get command line arguments
//------------------------------------------------------------------------------------
void
parse_command_line (const int     argc,
                    char *const  *argv)
{
   if (argc < 3)
   {
      std::cerr << "Not enough parameters given\n";
      exit (1);
   }
   
   std::list<std::string> args;
   for (int i=1; i<argc; ++i)
      args.push_back (argv[i]);
   
   run_mode = norun;
   restart = false;
   parameter_file = "null";
   
   while (args.size())
   {
      if (args.front() == std::string("-steady"))
      {
         run_mode = steady;
         args.pop_front ();
      }
      else if (args.front() == std::string("-unsteady"))
      {
         run_mode = unsteady;
         args.pop_front ();
      }
      else if (args.front() == std::string("-restart"))
      {
         restart = true;
         args.pop_front ();
      }
      else if (args.front() == std::string("-p"))
      {
         if (args.size() == 1)
         {
            std::cerr << "Error: flag '-p' must be followed by the\n"
                      << "       name of a parameter file."
                      << std::endl;
            exit (1);
         }
         args.pop_front ();
         parameter_file = args.front ();
         args.pop_front ();
      }
      else
      {
         std::cerr << "Unknown parameter: " << args.front() << std::endl;
         exit(1);
      }
   }
   
   if(run_mode == norun)
   {
      std::cerr << "Run mode must be specified\n";
      std::cerr << "    -steady or -unsteady\n";
      exit(1);
   }

   if(parameter_file == "null")
   {
      std::cerr << "Parameter file must be specified\n";
      std::cerr << "          -p filename\n";
      exit(1);
   }
}

//------------------------------------------------------------------------------
// Declare parameters for simulation
//------------------------------------------------------------------------------
void declare_parameters(ParameterHandler& prm)
{
   prm.declare_entry("reynolds no", "1.0", Patterns::Double(0.0),
                     "Reynolds number");
   prm.declare_entry("reference length", "1.0", Patterns::Double(0.0),
                     "Reference length");
   prm.declare_entry("reference velocity", "1.0", Patterns::Double(0.0),
                     "Reference velocity");
   
   prm.declare_entry("pressure degree", "1", Patterns::Integer(1,6),
                     "Pressure polynomial degree");
   
   prm.declare_entry("mesh file", "turek.msh",
                     Patterns::Anything(),
                     "Mesh file name");
   prm.declare_entry("linear solver", "umfpack",
                     Patterns::Selection("umfpack"),
                     "Linear solver");

   prm.declare_entry("time step", "0.1", Patterns::Double(0.0),
                     "Time step");
   prm.declare_entry("final time", "1.0e20", Patterns::Double(0.0),
                     "Final time");
   prm.declare_entry("no of iterations", "1000000", Patterns::Integer(0),
                     "Number of iterations");
}

//------------------------------------------------------------------------------------
// Initial condition. Sets parabolic velocity profile which is x direction only
// Pressure is set to zero which is not important.
//------------------------------------------------------------------------------------
template <int dim>
class InitialCondition : public Function<dim>
{
public:
   InitialCondition () : Function<dim>(dim+1) {}
   double value (const Point<dim>   &p,
                 const unsigned int  component = 0) const override;
   void vector_value (const Point<dim> &p,
                      Vector<double>   &value) const override;
};

template <int dim>
double
InitialCondition<dim>::value (const Point<dim>  &p,
                            const unsigned int component) const
{
   Assert (component < this->n_components,
           ExcIndexRange (component, 0, this->n_components));
   if (component == 0)
      return (1.5/std::pow(0.205,2))*(p[1]+0.2)*(0.21-p[1]);
   else
      return 0.0;
}

template <int dim>
void
InitialCondition<dim>::vector_value (const Point<dim> &p,
                                   Vector<double>   &values) const
{
   for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = InitialCondition<dim>::value (p, c);
}

//------------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------------
template <int dim>
class NS
{
public:
   NS (ParameterHandler &prm);
   ~NS() {};
   void run ();
   
private:
   struct AssemblyScratchData
   {
      AssemblyScratchData (const unsigned int degree,
                           const FESystem<dim> &fe,
                           MappingQ<dim>& mapping,
                           const unsigned int order);
      AssemblyScratchData (const AssemblyScratchData &scratch_data);
      FEValues<dim>     fe_values;
      double            a0, a1;
   };
   
   struct AssemblyCopyData
   {
      FullMatrix<double>                   local_matrix;
      Vector<double>                       local_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
   };
   
   void run_steady ();
   void run_unsteady ();
   void make_grid_dofs ();
   void assemble_mass_matrix ();
   void assemble_matrix (unsigned int order);
   void local_assemble_system
      (const typename DoFHandler<dim>::active_cell_iterator &cell,
       AssemblyScratchData                                  &scratch_data,
       AssemblyCopyData                                     &copy_data);
   void copy_local_to_global (const AssemblyCopyData &copy_data);
   void assemble_matrix_and_rhs (unsigned int order);
   void solve ();
   void compute_vorticity ();
   void output_results() const;
   
   ParameterHandler           *parameters;
   unsigned int               degree;
   FESystem<dim>              fe;
   FE_Q<dim>                  fe_scalar;
   Triangulation<dim>         triangulation;
   DoFHandler<dim>            dof_handler;
   DoFHandler<dim>            dof_handler_scalar;
   MappingQ<dim>              mapping;
   
   AffineConstraints<double>  constraints;
   BlockSparsityPattern       sparsity_pattern;
   BlockSparseMatrix<double>  system_matrix_constant;
   BlockSparseMatrix<double>  system_matrix;
   BlockVector<double>        solution0, solution1, solution2;
   BlockVector<double>        system_rhs;
   
   SparsityPattern            sparsity_pattern_scalar;
   SparseMatrix<double>       mass_matrix;
   Vector<double>             vorticity;
   SparseDirectUMFPACK        vorticity_solver;
   
   // Parameters
   double                     dt, Uref, Lref, Re, viscosity, final_time;
};

//------------------------------------------------------------------------------------
// Constructor of the main class.
//------------------------------------------------------------------------------------
template <int dim>
NS<dim>::NS (ParameterHandler &prm)
:
   parameters (&prm),
   degree (prm.get_integer("pressure degree")),
   fe( FE_Q<dim>(QGaussLobatto<1>(degree+2)), dim,
       FE_Q<dim>(QGaussLobatto<1>(degree+1)),   1),
   fe_scalar (FE_Q<dim>(QGaussLobatto<1>(degree+2))),
   dof_handler (triangulation),
   dof_handler_scalar (triangulation),
   mapping (degree+1)
{
   dt = parameters->get_double("time step");
   Re = parameters->get_double("reynolds no");
   Uref = parameters->get_double("reference velocity");
   Lref = parameters->get_double("reference length");
   final_time = parameters->get_double("final time");
   std::string grid_file = parameters->get("mesh file");
   
   std::cout << "Reading grid from " << grid_file << std::endl;
   GridIn<dim> grid_in;
   grid_in.attach_triangulation (triangulation);
   std::ifstream input_file (grid_file.c_str());
   grid_in.read_msh (input_file);
   
   viscosity = Uref*Lref/Re;

   // Set cylinder boundary description
   Point<dim> center (0.0, 0.0);
   static const SphericalManifold<dim> boundary_description (center);
   triangulation.set_manifold (2, boundary_description);
   
   std::string grid_output_file = "grid.eps";
   std::ofstream grid_output(grid_output_file);
   GridOut grid_out;
   grid_out.write_eps (triangulation, grid_output);
   std::cout << "Saved grid into " << grid_output_file << std::endl;
}

//------------------------------------------------------------------------------------
// Allocate memory for matrices and vectors.
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::make_grid_dofs()
{
   dof_handler.distribute_dofs (fe);
   
   DoFRenumbering::Cuthill_McKee (dof_handler);
   std::vector<unsigned int> block_component (dim+1,0);
   block_component[dim] = 1;
   DoFRenumbering::component_wise (dof_handler, block_component);

   std::vector<types::global_dof_index> dofs_per_block (2);
   dofs_per_block = DoFTools::count_dofs_per_fe_block (dof_handler, block_component);
   const unsigned int n_u = dofs_per_block[0],
                      n_p = dofs_per_block[1];
   std::cout << "   Number of active cells: "
             << triangulation.n_active_cells()
             << std::endl
             << "   Number of degrees of freedom: "
             << dof_handler.n_dofs()
             << " (" << n_u << '+' << n_p << ')'
             << std::endl;
   
   {
      BlockDynamicSparsityPattern csp (2,2);
      csp.block(0,0).reinit (n_u, n_u);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(1,1).reinit (n_p, n_p);
      csp.collect_sizes();
      DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
      sparsity_pattern.copy_from (csp);
   }
   
   system_matrix_constant.reinit (sparsity_pattern);
   system_matrix.reinit (sparsity_pattern);
   
   solution0.reinit (2);
   solution0.block(0).reinit (n_u);
   solution0.block(1).reinit (n_p);
   solution0.collect_sizes ();

   solution1.reinit (2);
   solution1.block(0).reinit (n_u);
   solution1.block(1).reinit (n_p);
   solution1.collect_sizes ();

   solution2.reinit (2);
   solution2.block(0).reinit (n_u);
   solution2.block(1).reinit (n_p);
   solution2.collect_sizes ();

   system_rhs.reinit (2);
   system_rhs.block(0).reinit (n_u);
   system_rhs.block(1).reinit (n_p);
   system_rhs.collect_sizes ();
   
   // These are needed for computing vorticity
   dof_handler_scalar.distribute_dofs (fe_scalar);
   DoFRenumbering::Cuthill_McKee (dof_handler_scalar);

   {
      DynamicSparsityPattern csp (dof_handler_scalar.n_dofs());
      DoFTools::make_sparsity_pattern (dof_handler_scalar, csp);
      sparsity_pattern_scalar.copy_from (csp);
   }
   
   mass_matrix.reinit (sparsity_pattern_scalar);
   vorticity.reinit (dof_handler_scalar.n_dofs());
   std::cout << "   Number of vorticity dofs: "
             << dof_handler_scalar.n_dofs()
             << std::endl;
}

//------------------------------------------------------------------------------------
// Assemble mass matrix which is used for vorticity projection.
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::assemble_mass_matrix ()
{
   mass_matrix = 0;
   
   QGauss<dim>   quadrature_formula(degree+2);
   FEValues<dim> fe_values (mapping, fe_scalar, quadrature_formula,
                            update_values    |
                            update_JxW_values);
   const unsigned int   dofs_per_cell   = fe_scalar.dofs_per_cell;
   const unsigned int   n_q_points      = quadrature_formula.size();
   FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler_scalar.begin_active(),
      endc = dof_handler_scalar.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      local_matrix = 0;
      
      for (unsigned int q=0; q<n_q_points; ++q)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<=i; ++j)
               local_matrix(i,j) +=   fe_values.shape_value(i,q)
                                    * fe_values.shape_value(j,q)
                                    * fe_values.JxW(q);
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=i+1; j<dofs_per_cell; ++j)
            local_matrix(i,j) = local_matrix(j,i);
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            mass_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             local_matrix(i,j));
   }
}

//------------------------------------------------------------------------------------
// Assemble part of matrix which is independent of time. We assume that time step
// is fixed. Only the convective terms are not assembled here.
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::assemble_matrix (unsigned int order)
{
   double a2;
   if(order == 0)
      a2 = 0.0;
   else if(order == 1)
      a2 = 1.0;
   else if(order == 2)
      a2 = 1.5;
   else
      Assert(false, ExcMessage("Not implemented"));
   
   system_matrix_constant = 0;
   
   QGauss<dim>   quadrature_formula(degree+2);
   FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                            update_values    |
                            update_gradients |
                            update_JxW_values);
   const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
   const unsigned int   n_q_points      = quadrature_formula.size();
   FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   
   const FEValuesExtractors::Vector velocities (0);
   const FEValuesExtractors::Scalar pressure (dim);
   
   std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
   std::vector<double>                  div_phi_u     (dofs_per_cell);
   std::vector<double>                  phi_p         (dofs_per_cell);
   std::vector<Tensor<1,dim> >          phi_u         (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      local_matrix = 0;

      for (unsigned int q=0; q<n_q_points; ++q)
      {
         for (unsigned int k=0; k<dofs_per_cell; ++k)
         {
            symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
            div_phi_u[k]     = fe_values[velocities].divergence (k, q);
            phi_p[k]         = fe_values[pressure].value (k, q);
            phi_u[k]         = fe_values[velocities].value (k, q);
         }
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            for (unsigned int j=0; j<=i; ++j)
            {
               local_matrix(i,j) += ((a2/dt) * phi_u[i] * phi_u[j]
                                     + 2.0 * viscosity * symgrad_phi_u[i] * symgrad_phi_u[j]
                                     - div_phi_u[i] * phi_p[j]
                                     - phi_p[i] * div_phi_u[j]
                                     )
                                     * fe_values.JxW(q);
            }
         }
      }
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=i+1; j<dofs_per_cell; ++j)
            local_matrix(i,j) = local_matrix(j,i);
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix_constant.add (local_dof_indices[i],
                                        local_dof_indices[j],
                                        local_matrix(i,j));
   }
}

//------------------------------------------------------------------------------------
// Data used in parallel assembly
//------------------------------------------------------------------------------------
template <int dim>
NS<dim>::AssemblyScratchData::
AssemblyScratchData (const unsigned int   degree,
                     const FESystem<dim>& fe,
                     MappingQ<dim>&       mapping,
                     const unsigned int   order)
:
fe_values (mapping, fe,
           QGauss<dim>(degree+2),
           update_values   | update_gradients |
           update_JxW_values)
{
   if(order == 0)
   {
      a0 = 0.0;
      a1 = 0.0;
   }
   else if(order == 1)
   {
      a0 =-1.0;
      a1 = 0.0;
   }
   else if(order == 2)
   {
      a0 =  0.5;
      a1 = -2.0;
   }
   else
      Assert(false, ExcMessage("Not implemented"));
}

//------------------------------------------------------------------------------------
// Data used in parallel assembly
//------------------------------------------------------------------------------------
template <int dim>
NS<dim>::AssemblyScratchData::
AssemblyScratchData (const AssemblyScratchData &scratch_data)
:
fe_values (scratch_data.fe_values.get_mapping(),
           scratch_data.fe_values.get_fe(),
           scratch_data.fe_values.get_quadrature(),
           update_values   | update_gradients |
           update_JxW_values),
a0 (scratch_data.a0),
a1 (scratch_data.a1)
{}

//------------------------------------------------------------------------------------
// Assemble matrix/rhs on one cell
//------------------------------------------------------------------------------------
template <int dim>
void
NS<dim>::
local_assemble_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                       AssemblyScratchData                                  &scratch_data,
                       AssemblyCopyData                                     &copy_data)
{
   const double a0 = scratch_data.a0;
   const double a1 = scratch_data.a1;
   
   FEValues<dim>&      fe_values         = scratch_data.fe_values;
   const unsigned int  dofs_per_cell     = fe.dofs_per_cell;
   const unsigned int  n_q_points        = scratch_data.fe_values.get_quadrature().size();
   FullMatrix<double>& local_matrix      = copy_data.local_matrix;
   Vector<double>&     local_rhs         = copy_data.local_rhs;
   std::vector<types::global_dof_index>&
                       local_dof_indices = copy_data.local_dof_indices;
   
   local_matrix.reinit (dofs_per_cell, dofs_per_cell);
   local_rhs.reinit (dofs_per_cell);
   local_dof_indices.resize (dofs_per_cell);
   
   const FEValuesExtractors::Vector velocities (0);
   const FEValuesExtractors::Scalar pressure (dim);
   
   std::vector<Tensor<2,dim> >  grad_phi_u (dofs_per_cell);
   std::vector<Tensor<1,dim> >  phi_u      (dofs_per_cell);
   std::vector<Tensor<1,dim> >  velocity   (n_q_points, Tensor<1,dim>());
   
   fe_values.reinit (cell);
   fe_values[velocities].get_function_values (solution2, velocity);

   cell->get_dof_indices (local_dof_indices);
   
   for (unsigned int q=0; q<n_q_points; ++q)
   {
      for (unsigned int k=0; k<dofs_per_cell; ++k)
      {
         grad_phi_u[k] = fe_values[velocities].gradient (k, q);
         phi_u[k]      = fe_values[velocities].value (k, q);
      }
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         for (unsigned int j=0; j<dofs_per_cell; ++j)
         {
            local_matrix(i,j) += (grad_phi_u[j] * velocity[q]) * phi_u[i]
                                 * fe_values.JxW(q);
            local_rhs(i) += (1.0/dt) *
                            (- a1 * solution1(local_dof_indices[j]) 
                             - a0 * solution0(local_dof_indices[j]))
                            * phi_u[i] * phi_u[j] * fe_values.JxW(q);
         }
   }
}

//------------------------------------------------------------------------------------
// Copy cell assembly to global position
//------------------------------------------------------------------------------------
template <int dim>
void
NS<dim>::copy_local_to_global (const AssemblyCopyData &copy_data)
{
   for (unsigned int i=0; i<copy_data.local_dof_indices.size(); ++i)
   {
      for (unsigned int j=0; j<copy_data.local_dof_indices.size(); ++j)
         system_matrix.add (copy_data.local_dof_indices[i],
                            copy_data.local_dof_indices[j],
                            copy_data.local_matrix(i,j));
      system_rhs(copy_data.local_dof_indices[i]) += copy_data.local_rhs(i);
   }
}

//------------------------------------------------------------------------------------
// Assemble matrix and rhs. Only the convective term is assembled in the matrix.
// The other terms are independent of time and have been assembled before. This is
// thread parallel using WorkStream.
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::assemble_matrix_and_rhs (unsigned int order)
{
   // use solution2 for extrapolated velocity
   if(order == 0)
      ; // nothing to do; we use solution2
   else if(order == 1)
      solution2.block(0) = solution0.block(0);
   else
      // solution2 = 2 * solution1 - solution0
      // Step 1: solution2 = 2 * solution1
      solution2.block(0).equ(2.0, solution1.block(0));
      // Step 2: solution2 += (-1) * solution0
      solution2.block(0).add(-1.0, solution0.block(0));
   
   system_matrix.copy_from(system_matrix_constant);
   system_rhs    = 0;
   
   WorkStream::run(dof_handler.begin_active(),
                   dof_handler.end(),
                   *this,
                   &NS::local_assemble_system,
                   &NS::copy_local_to_global,
                   AssemblyScratchData(degree, fe, mapping, order),
                   AssemblyCopyData());
   
   const FEValuesExtractors::Vector velocities (0);
   
   // Apply boundary conditions
   std::map<types::global_dof_index,double> boundary_values;
   VectorTools::interpolate_boundary_values (mapping,
                                             dof_handler,
                                             1,
                                             InitialCondition<dim>(),
                                             boundary_values,
                                             fe.component_mask(velocities));
   VectorTools::interpolate_boundary_values (mapping,
                                             dof_handler,
                                             2,
                                             Functions::ZeroFunction<dim>(dim+1),
                                             boundary_values,
                                             fe.component_mask(velocities));
   VectorTools::interpolate_boundary_values (mapping,
                                             dof_handler,
                                             3,
                                             Functions::ZeroFunction<dim>(dim+1),
                                             boundary_values,
                                             fe.component_mask(velocities));
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution2,
                                       system_rhs);
}

//------------------------------------------------------------------------------------
// Solve linear system. At present only direct solver is implemented.
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::solve()
{
   SparseDirectUMFPACK  solver;
   solver.initialize (system_matrix);
   solver.vmult (solution2, system_rhs);
}

//------------------------------------------------------------------------------------
// Compute vorticity by doing an L2 projection. Vorticity space has same degree as
// velocity space.
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::compute_vorticity ()
{
   static unsigned int status = 0;
   
   if(status == 0)
   {
      assemble_mass_matrix ();
      vorticity_solver.initialize (mass_matrix);
      status = 1;
   }
   
   Vector<double> vorticity_rhs (dof_handler_scalar.n_dofs());
   
   QGauss<dim>   quadrature_formula(degree+2);
   FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                            update_gradients);
   FEValues<dim> fe_values_vorticity (mapping, fe_scalar, quadrature_formula,
                            update_values | update_JxW_values);
   const unsigned int   dofs_per_cell   = fe_scalar.dofs_per_cell;
   const unsigned int   n_q_points      = quadrature_formula.size();
   Vector<double>       local_rhs (dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   
   const FEValuesExtractors::Vector velocities (0);
   std::vector<typename FEValuesViews::Vector<dim>::curl_type> vorticity_values (n_q_points);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end(),
      cell_vorticity = dof_handler_scalar.begin_active();
   for (; cell!=endc; ++cell, ++cell_vorticity)
   {
      fe_values.reinit (cell);
      fe_values_vorticity.reinit (cell_vorticity);

      fe_values[velocities].get_function_curls (solution2, vorticity_values);
      
      local_rhs    = 0;
      
      for (unsigned int q=0; q<n_q_points; ++q)
         for (unsigned int i=0; i<dofs_per_cell; ++i)
            local_rhs(i) +=   vorticity_values[q][0]
                            * fe_values_vorticity.shape_value(i,q)
                            * fe_values_vorticity.JxW(q);
      
      cell_vorticity->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
         vorticity_rhs(local_dof_indices[i]) += local_rhs(i);
   }
   
   vorticity_solver.vmult (vorticity, vorticity_rhs);
}

//------------------------------------------------------------------------------------
// Save solution to vtk file
//------------------------------------------------------------------------------------
template <int dim>
void
NS<dim>::output_results ()  const
{
   static unsigned int cycle = 0;

   std::vector<std::string> solution_names (dim, "velocity");
   solution_names.push_back ("pressure");
   
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
   data_component_interpretation
   (dim, DataComponentInterpretation::component_is_part_of_vector);
   data_component_interpretation
   .push_back (DataComponentInterpretation::component_is_scalar);
   
   DataOut<dim> data_out;
   data_out.add_data_vector (dof_handler,
                             solution2,
                             solution_names,
                             data_component_interpretation);
   data_out.add_data_vector (dof_handler_scalar,
                             vorticity,
                             "vorticity");

   data_out.build_patches (mapping, degree+1);
   
   std::string filename = "solution-" 
                          + Utilities::int_to_string (cycle, 3) + ".vtk";
   std::ofstream output (filename);
   data_out.write_vtk (output);
   std::cout << "Wrote solution into " << filename << std::endl;
   
   ++cycle;
}

//------------------------------------------------------------------------------------
// Run steady state computation. The convection term is treated linearly and a Picard
// iteration is performed on it.
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::run_steady ()
{
   unsigned int order = 0;
   assemble_matrix (order);
   
   // Set initial condition
   std::cout << "Setting initial condition ..." << std::endl;
   VectorTools::interpolate(mapping, dof_handler,
                            InitialCondition<dim>(), solution0);
   solution1 = solution0;
   solution2 = solution0;
   compute_vorticity ();
   output_results ();
   
   unsigned int iter = 0;
   
   while (iter < 10)
   {
      // Assemble matrix and rhs
      assemble_matrix_and_rhs (order);
      
      // solve
      solve ();
      
      ++iter;
      std::cout << iter << std::endl;
      
      compute_vorticity ();
      output_results ();
   }
   
   // save solution to file
   std::ofstream output_file("steady.dat");
   solution2.block_write (output_file);
   
}

//------------------------------------------------------------------------------------
// Run unsteady computations
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::run_unsteady ()
{
   unsigned int order = 1;
   assemble_matrix (order);
   
   // Set initial condition
   if(restart == false)
   {
      std::cout << "Setting initial condition ..." << std::endl;
      VectorTools::interpolate(mapping, dof_handler,
                               InitialCondition<dim>(), solution0);
   }
   else
   {
      std::cout << "Read initial condition from file steady.dat\n";
      std::ifstream restart_file("steady.dat");
      solution0.block_read (restart_file);
   }
   
   solution1 = solution0;
   solution2 = solution0;
   
   compute_vorticity ();
   output_results ();
   
   double time = 0;
   unsigned int iter = 0;
   
   while (time < final_time)
   {
      // Assemble matrix and rhs
      assemble_matrix_and_rhs (order);
      
      // solve
      solve ();
      
      time += dt;
      ++iter;
      std::cout << iter << "  " << time << "  " << std::endl;
      
      if(iter == 1)
      {
         solution1 = solution2;

         order = 2;
         assemble_matrix (order);
      }
      else
      {
         solution0 = solution1;
         solution1 = solution2;
      }
      
      if(iter%10 == 0)
      {
         compute_vorticity ();
         output_results ();
      }
   }
   
}

//------------------------------------------------------------------------------------
// Run the code in specified mode
//------------------------------------------------------------------------------------
template <int dim>
void NS<dim>::run ()
{
   make_grid_dofs ();

   if(run_mode == steady)
      run_steady ();
   else if(run_mode == unsteady)
      run_unsteady ();
   else
      AssertThrow(false, ExcMessage("Unknown run mode"));
}

//------------------------------------------------------------------------------------
// Main function
//------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   try
   {
      deallog.depth_console (0);

      parse_command_line (argc, argv);
      
      ParameterHandler prm;
      declare_parameters (prm);
      prm.parse_input (parameter_file);
      
      NS<2> ns_problem (prm);
      ns_problem.run();
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
   }
   return 0;
}
