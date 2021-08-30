/*
 Solve 2d/3d laplace equation
 -Laplace(u) = f(x,y) in (0,1)x(0,1)
 Exact solution is u = sin(2*pi*x) * sin(2*pi*y)
 f is obtained from exact solution.
 Boundary condition is dirichlet and taken from exact solution.
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
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>

#include <fstream>
#include <iostream>


using namespace dealii;

//------------------------------------------------------------------------------
template <int dim>
Tensor<1,dim> velocity(const Point<dim>& p)
{
   Tensor<1,dim> v;
   v[0] =  p[1];
   v[1] = -p[0];
   return v;
}

//------------------------------------------------------------------------------
// RHS function f in Poisson equation
template <int dim>
class RightHandSide : public Function<dim>
{
public:
   RightHandSide () : Function<dim>() {}

   virtual double value (const Point<dim>&   p,
                         const unsigned int  component = 0) const;
};

// RHS in 2-D
template <>
double RightHandSide<2>::value (const Point<2>& /*p*/,
                                const unsigned int /*component*/) const
{
   return 0.0;
}

//------------------------------------------------------------------------------
// Boundary condition
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
   BoundaryValues () : Function<dim>() {}

   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
};

template <int dim>
double BoundaryValues<dim>::value (const Point<dim>& p,
                                   const unsigned int /*component*/) const
{
   const double EPS = 1.0e-13;
   const double x = p[0];
   const double y = p[1];
   // top sides
   if(fabs(y - 1.0) < EPS)
   {
      return 0.0;
   }
   else if(fabs(x) < EPS) // left side
   {
      if(fabs(y - 0.5) < 0.25)
         return 1.0;
      else
         return 0.0;
   }
   else
   {
      AssertThrow(false, ExcMessage("Wrong bc"));
   }
}

//------------------------------------------------------------------------------
template <int dim>
class LaplaceProblem
{
public:
   LaplaceProblem (int degree, std::string method);
   void run ();

private:
   void make_grid_and_dofs ();
   void assemble_system ();
   void solve ();
   void output_results () const;

   std::string            method;
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
LaplaceProblem<dim>::LaplaceProblem (int degree, std::string method)
:
method(method),
fe (degree),
dof_handler (triangulation)
{
   std::cout << "Degree = " << degree << std::endl;
   std::cout << "Method = " << method << std::endl;
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
   GridGenerator::hyper_cube (triangulation, 0, 1, true);
   triangulation.refine_global (7);

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
   system_rhs    = 0;

   QGauss<dim>  quadrature_formula(2*fe.degree);
   const RightHandSide<dim> right_hand_side;
   FEValues<dim> fe_values (fe, quadrature_formula,
                            update_values   | update_gradients |
                            update_quadrature_points | update_JxW_values);

   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();

   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
   Vector<double>       cell_rhs (dofs_per_cell);
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   std::vector<double>  rhs_values (n_q_points);

   const int supg = (method == "supg") ? 1 : 0;

   for (const auto &cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;
      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);

      const double h = cell->diameter();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
         const Tensor<1,dim> vel = velocity<dim>(fe_values.quadrature_point(q_point));
         const double mod_vel = vel.norm();
         const double f_supg = supg * h / mod_vel;
         for (unsigned int i=0; i<dofs_per_cell; ++i)
         {
            const double v_grad_phi_i = vel * fe_values.shape_grad(i, q_point);
            // phi_i + (h/|vel|) * v . grad(phi_i)
            const double test_i = (fe_values.shape_value(i, q_point)
                                  + f_supg * v_grad_phi_i);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
               const double v_grad_phi_j = vel * fe_values.shape_grad(j, q_point);
               cell_matrix(i,j) += test_i * v_grad_phi_j * fe_values.JxW (q_point);
            }

            cell_rhs(i) += (test_i *                   // test fun
                            rhs_values[q_point] *      // rhs func
                            fe_values.JxW (q_point));  // det(J) * w
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

   // boundary condition
   std::map<types::global_dof_index,double> boundary_values;
   VectorTools::interpolate_boundary_values (dof_handler,
                                             0,  // left boundary
                                             BoundaryValues<dim>(),
                                             boundary_values);
   VectorTools::interpolate_boundary_values (dof_handler,
                                             3,  // top boundary
                                             BoundaryValues<dim>(),
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
   std::cout << "Solving matrix problem" << std::endl;
   SparseDirectUMFPACK  solver;
   solver.initialize(system_matrix);
   solver.vmult (solution, system_rhs);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::output_results () const
{
   DataOut<dim> data_out;

   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches (fe.degree);
   std::string filename = method + ".vtk";
   std::ofstream output (filename);
   data_out.write_vtk (output);
   std::cout << "Saved solution into " << filename << std::endl;
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
   {
      LaplaceProblem<2> problem (degree, "galerkin");
      problem.run ();
      std::cout << "------------------------------------\n";
   }
   {
      LaplaceProblem<2> problem (degree, "supg");
      problem.run ();
      std::cout << "------------------------------------\n";
   }

   return 0;
}
