/* Solve 1d laplace equation
   -u_xx = 4 sin(2x)  for x in (0,1)
   Exact solution is u = sin(2x)
   Boundary condition is dirichlet and taken from exact solution
*/
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <base/logstream.h>

using namespace dealii;

//------------------------------------------------------------------------------
template <int dim>
class RightHandSide : public Function<dim>
{
public:
   RightHandSide () : Function<dim>() {}
   
   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
};

template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int /*component*/) const
{
   return 4*sin(2*p[0]);
}

//------------------------------------------------------------------------------
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
   BoundaryValues () : Function<dim>() {}
   
   virtual double value (const Point<dim>   &p,
                         const unsigned int  component = 0) const;
};

template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
   return sin(2*p[0]);
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

    int                  degree;
    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};


//------------------------------------------------------------------------------
template <int dim>
LaplaceProblem<dim>::LaplaceProblem (int degree) :
    degree (degree),
    fe (degree),
    dof_handler (triangulation)
{}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
    GridGenerator::hyper_cube (triangulation, 0, 1);
    triangulation.refine_global (8);

    std::cout << "   Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "   Total number of cells: "
              << triangulation.n_cells()
              << std::endl;

    dof_handler.distribute_dofs (fe);

    std::cout << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    sparsity_pattern.copy_from(c_sparsity);

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
    QGauss<dim>  quadrature_formula(2*degree);

    const RightHandSide<dim> right_hand_side;

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
        fe_values.reinit (cell);
        cell_matrix = 0;
        cell_rhs = 0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                         fe_values.shape_grad (j, q_point) *
                                         fe_values.JxW (q_point));

                cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                                right_hand_side.value (fe_values.quadrature_point (q_point)) *
                                fe_values.JxW (q_point));
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

	// left boundary condition
    std::map<unsigned int,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
            0,
            BoundaryValues<dim>(),
            boundary_values);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
	// right boundary condition
    VectorTools::interpolate_boundary_values (dof_handler,
            1,
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
    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              cg (solver_control);
    cg.solve (system_matrix, solution, system_rhs,
              PreconditionIdentity());

    std::cout << "   " << solver_control.last_step()
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

    data_out.build_patches (degree);

    std::ofstream output ("solution.gnuplot");
    data_out.write_gnuplot (output);
}

//------------------------------------------------------------------------------
template <int dim>
void LaplaceProblem<dim>::run ()
{
    std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;

    make_grid_and_dofs();
    assemble_system ();
    solve ();
    output_results ();
}

//------------------------------------------------------------------------------
int main ()
{
    deallog.depth_console (0);
    {
        int degree = 2;
        LaplaceProblem<1> laplace_problem_1d (degree);
        laplace_problem_1d.run ();
    }

    return 0;
}

