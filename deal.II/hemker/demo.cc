#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

#include <fstream>
#include <iostream>

using namespace dealii;

const Tensor<1,2> velocity({1.0, 0.0});
double epsilon = 1.0e-3;
int supg = 0;
int degree = 1;
int nrefine = 0;

//------------------------------------------------------------------------------
template <int dim>
class LaplaceProblem
{
public:
   LaplaceProblem (int degree);
   void run ();
   
private:
   void make_grid_and_dofs ();
   void assemble_system ();
   void solve ();
   void output_results () const;
   
   Triangulation<dim>     triangulation;
   FE_Q<dim>              fe;
   DoFHandler<dim>        dof_handler;
   
   SparsityPattern        sparsity_pattern;
   SparseMatrix<double>   system_matrix;
   
   Vector<double>         solution;
   Vector<double>         system_rhs;
};

//------------------------------------------------------------------------------
template <int dim>
LaplaceProblem<dim>::LaplaceProblem (int degree) :
fe (degree),
dof_handler (triangulation)
{}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
   GridIn<dim> grid_in;
   grid_in.attach_triangulation(triangulation);
   std::ifstream gfile("hemker.msh");
   AssertThrow(gfile.is_open(), ExcMessage("Grid file not found"));
   grid_in.read_msh(gfile);

   // circle has boundary id = 2, see geo file. Set its manifold id = 1
   triangulation.set_all_manifold_ids_on_boundary(2, 1);
   triangulation.set_manifold(1, SphericalManifold<dim>());

   if(nrefine > 0) triangulation.refine_global(nrefine);
   
   std::cout
   << "Problem size:" << std::endl
   << "   Number of active cells: "
   << triangulation.n_active_cells()
   << std::endl
   << "   Total number of cells: "
   << triangulation.n_cells()
   << std::endl;
   
   dof_handler.distribute_dofs (fe);
   
   std::cout
   << "   Number of degrees of freedom: "
   << dof_handler.n_dofs()
   << std::endl;
   
   DynamicSparsityPattern dsp(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, dsp);
   sparsity_pattern.copy_from(dsp);
   
   system_matrix.reinit (sparsity_pattern);
   solution.reinit (dof_handler.n_dofs());
   system_rhs.reinit (dof_handler.n_dofs());
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
   system_matrix = 0;
   system_rhs = 0;
   QGauss<dim>    cell_quadrature_formula(2*fe.degree);
   
   FEValues<dim> fe_values (fe, cell_quadrature_formula,
                            update_values   | update_gradients |
                            update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = cell_quadrature_formula.size();
   
   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

   
   for (const auto &cell : dof_handler.active_cell_iterators())
   {
      const double h = cell->diameter();
      const double tau = supg * h / velocity.norm();

      fe_values.reinit (cell);
      cell_matrix = 0;
      
      // Cell integral
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
               cell_matrix(i,j) += 
                  (epsilon *
                   fe_values.shape_grad (i, q_point) *
                   fe_values.shape_grad (j, q_point)
                   +
                   (velocity * 
                   fe_values.shape_grad(j,q_point)) *
                   fe_values.shape_value(i,q_point)
                   +
                   tau * 
                   (velocity * fe_values.shape_grad(i,q_point)) *
                   (velocity * fe_values.shape_grad(j,q_point))) *
                   fe_values.JxW (q_point);
         }
      }
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));
      }
   }
   
   // Dirichlet boundary condition
   std::map<types::global_dof_index,double> boundary_values;
   // Left boundary
   VectorTools::interpolate_boundary_values(dof_handler,
                                            1,
                                            Functions::ZeroFunction<dim>(),
                                            boundary_values);
   // Cylinder
   VectorTools::interpolate_boundary_values(dof_handler,
                                            2,
                                            Functions::ConstantFunction<dim>(1.0),
                                            boundary_values);
   // Modify matrix and rhs
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::solve ()
{
   SparseDirectUMFPACK solver;
   solver.initialize(system_matrix);
   solver.vmult(solution, system_rhs);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::output_results () const
{
   DataOut<dim> data_out;
   
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches (fe.degree);
   std::ofstream output ("solution.vtu");
   data_out.write_vtu (output);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::run ()
{
   make_grid_and_dofs();
   assemble_system ();
   solve ();
   output_results ();
}

//------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
   for(int i=1; i<argc; i+=2)
   {
      if (std::string(argv[i]) == "-supg")
      {
         supg = Utilities::string_to_int(argv[i+1]);
         AssertThrow(supg==0 || supg==1, ExcMessage("Invalid value of -supg"));
      }
      else if (std::string(argv[i]) == "-eps")
         epsilon = Utilities::string_to_double(argv[i+1]);
      else if (std::string(argv[i]) == "-degree")
         degree = Utilities::string_to_int(argv[i+1]);
      else if (std::string(argv[i]) == "-nrefine")
         nrefine = Utilities::string_to_int(argv[i+1]);
      else
      {
         std::cout << "Unknown option given: " << argv[i] << std::endl;
         exit(0);
      }
   }

   std::cout << "degree  = " << degree << std::endl;
   std::cout << "epsilon = " << epsilon << std::endl;
   std::cout << "supg    = " << supg << std::endl;
   std::cout << "nrefine = " << nrefine << std::endl;

   if(supg == 1) AssertThrow(degree==1, ExcMessage("supg needs degree=1"));

   deallog.depth_console (0);
   LaplaceProblem<2> problem (degree);
   problem.run ();
   
   return 0;
}
