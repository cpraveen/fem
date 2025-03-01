#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <fstream>
#include <iostream>

using namespace dealii;

void test1()
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

void test2()
{
   Triangulation<2> triangulation;

   {
      Triangulation<2> tmp_tria;
      GridGenerator::subdivided_hyper_cube (tmp_tria, 4, 0, 1);
      GridGenerator::convert_hypercube_to_simplex_mesh(tmp_tria, triangulation);
   }

   MappingFE<2> mapping(FE_SimplexP<2>(1));

   const Point<2> center(0.5, 0.5);
   for (const auto &cell : triangulation.active_cell_iterators())
   {
      for (const auto v : cell->vertex_indices())
      {
         const double distance_from_center =
              center.distance(cell->vertex(v));
 
         if (distance_from_center < 0.25)
         {
            cell->set_refine_flag();
            break;
         }
      }
   }
 
   triangulation.execute_coarsening_and_refinement();

   {
      std::ofstream outfile("grid3.svg");
      GridOut().write_svg(triangulation, outfile);
      std::cout << "   Number of active cells: "
                << triangulation.n_active_cells()
                << std::endl
                << "   Total number of cells: "
                << triangulation.n_cells()
                << std::endl;
   }
}

int main()
{
   test1();
   test2();
}
