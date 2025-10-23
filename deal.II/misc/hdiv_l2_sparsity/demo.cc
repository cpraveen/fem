#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

using namespace dealii;

const int dim = 2;

int
main()
{
   const unsigned int degree = 1;
   const FE_RaviartThomas<dim> fe1(degree);
   const FE_DGQ<dim> fe2(degree);
   const FESystem<dim> fe(fe1, fe2);

   Triangulation<dim> triangulation;
   GridGenerator::hyper_cube(triangulation);
   triangulation.refine_global(2);

   DoFHandler<dim> dof_handler(triangulation);
   dof_handler.distribute_dofs(fe);
   DoFRenumbering::component_wise(dof_handler);

   AffineConstraints<double> constraints;
   constraints.clear();
   constraints.close();

   const std::vector<types::global_dof_index> dofs_per_component =
      DoFTools::count_dofs_per_fe_component(dof_handler);
   const unsigned int n_u = dofs_per_component[0],
                      n_p = dofs_per_component[dim];

   std::cout << "Number of active cells: " << triangulation.n_active_cells()
             << std::endl
             << "Total number of cells: " << triangulation.n_cells()
             << std::endl
             << "Number of degrees of freedom: " << dof_handler.n_dofs()
             << " (" << n_u << '+' << n_p << ')' << std::endl;

   const std::vector<types::global_dof_index> block_sizes = {n_u, n_p};

   // Full coupling
   {
      BlockDynamicSparsityPattern dsp(block_sizes, block_sizes);

      DoFTools::make_sparsity_pattern(
         dof_handler, dsp, constraints, false);

      std::ofstream spfile("sparsity1.svg");
      dsp.print_svg(spfile);
   }

   // remove coupling in lower-right block
   {
      BlockDynamicSparsityPattern dsp(block_sizes, block_sizes);

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
      for(unsigned int c = 0; c < dim + 1; ++c)
         for(unsigned int d = 0; d < dim + 1; ++d)
            if(!((c == dim) && (d == dim)))
               coupling[c][d] = DoFTools::always;
            else
               coupling[c][d] = DoFTools::none;

      DoFTools::make_sparsity_pattern(
         dof_handler, coupling, dsp, constraints, false);

      std::ofstream spfile("sparsity2.svg");
      dsp.print_svg(spfile);
   }
}
