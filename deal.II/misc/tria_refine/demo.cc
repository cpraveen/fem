#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <fstream>
#include <iostream>

using namespace dealii;

int main()
{
   Triangulation<2> triangulation;

   {
      Triangulation<2> tmp_tria;
      GridGenerator::subdivided_hyper_cube (tmp_tria, 4, 0, 1);
      GridGenerator::convert_hypercube_to_simplex_mesh(tmp_tria, triangulation);
   }

   MappingFE<2> mapping(FE_SimplexP<2>(1));

   {
      std::ofstream outfile("grid1.svg");
      GridOut().write_svg(triangulation, outfile);
      std::cout << "   Number of active cells: "
                << triangulation.n_active_cells()
                << std::endl
                << "   Total number of cells: "
                << triangulation.n_cells()
                << std::endl;
   }

   triangulation.refine_global();

   {
      std::ofstream outfile("grid2.svg");
      GridOut().write_svg(triangulation, outfile);
      std::cout << "   Number of active cells: "
                << triangulation.n_active_cells()
                << std::endl
                << "   Total number of cells: "
                << triangulation.n_cells()
                << std::endl;
   }
}
