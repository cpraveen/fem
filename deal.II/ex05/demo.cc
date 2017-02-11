/* 
 Solve 2d laplace equation
 -Laplace(u) = f(x) in (0,1)x(0,1)
 Exact solution is u = sin(2*pi*x) * sin(2*pi*y)
 f is obtained from exact solution.
 Boundary condition:
   dirichlet: on x=0 and x=1
   neumann  : on y=0 and y=1
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
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

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
                 const unsigned int  component = 0) const;
   Tensor<1,dim> gradient (const Point<dim>   &p,
                           const unsigned int  component = 0) const;
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
   
   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
};

template <>
double RightHandSide<2>::value (const Point<2> &p,
                                  const unsigned int /*component*/) const
{
   return 8*M_PI*M_PI*sin(2*M_PI*p[0])*sin(2*M_PI*p[1]);
}



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
   GridGenerator::hyper_cube (triangulation, 0, 1, true);
   triangulation.refine_global (5);
   
   std::cout
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
   QGauss<dim-1>  face_quadrature_formula(2*fe.degree);
   const RightHandSide<dim> right_hand_side;
   const ExactSolution<dim> exact_solution;
   
   FEValues<dim> fe_values (fe, cell_quadrature_formula,
                            update_values   | update_gradients |
                            update_quadrature_points | update_JxW_values);
   FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
                                     update_values   | update_normal_vectors |
                                     update_quadrature_points | update_JxW_values);
   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = cell_quadrature_formula.size();
   const unsigned int   n_face_q_points = face_quadrature_formula.size();
   
   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
   Vector<double>       cell_rhs (dofs_per_cell);
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;
      
      // Cell integral
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
         const double f = right_hand_side.value(fe_values.quadrature_point(q_point));
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
               cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                    fe_values.shape_grad (j, q_point) *
                                    fe_values.JxW (q_point));
            
            cell_rhs(i) += (f *
                            fe_values.shape_value (i, q_point) *
                            fe_values.JxW (q_point));
         }
      }
      
      // Loop over cell faces and assemble only if on Neumann part
      for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
         if(cell->face(f)->at_boundary() &&
            (cell->face(f)->boundary_id()==2 || cell->face(f)->boundary_id()==3))
         {
            fe_face_values.reinit(cell, f);
            for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
            {
               const double g = (exact_solution.gradient(fe_face_values.quadrature_point(q_point)) *
                                 fe_face_values.normal_vector(q_point));
               for (unsigned int i=0; i<dofs_per_cell; ++i)
                  cell_rhs(i) += (g *
                                  fe_face_values.shape_value(i, q_point) *
                                  fe_face_values.JxW(q_point));
               
            }
         }
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
         for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));
         
         system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }
   }
   
   // Dirichlet boundary condition
   std::map<unsigned int,double> boundary_values;
   // Left boundary
   VectorTools::interpolate_boundary_values (dof_handler,
                                             0,
                                             exact_solution,
                                             boundary_values);
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
   // Right boundary
   VectorTools::interpolate_boundary_values (dof_handler,
                                             1,
                                             exact_solution,
                                             boundary_values);
   MatrixTools::apply_boundary_values (boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::solve ()
{
   SolverControl           solver_control (1000, 1e-12);
   SolverCG<>              cg (solver_control);
   cg.solve (system_matrix, solution, system_rhs,
             PreconditionIdentity());
   
   std::cout
   << "   " << solver_control.last_step()
   << " CG iterations needed to obtain convergence."
   << std::endl;
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::output_results () const
{
   DataOut<dim> data_out;
   
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches (fe.degree);
   std::ofstream output ("solution.vtk");
   data_out.write_vtk (output);
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
int main ()
{
   deallog.depth_console (0);
   int degree = 1;
   LaplaceProblem<2> problem (degree);
   problem.run ();
   
   return 0;
}
