#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/fe/mapping_q.h>

#include <fstream>

using namespace dealii;

int main()
{
   Triangulation<2> triangulation;
   const double rin = 1.0;
   const double rout = 2.0;
   GridGenerator::hyper_shell(triangulation,
                              Point<2>(0.0, 0.0),
                              rin,
                              rout,
                              10,
                              true);
   GridOut grid_out;

   std::ofstream out1("grid_q1.gpl");
   grid_out.write_gnuplot(triangulation, out1); // uses MappingQ1

   std::ofstream out2("grid_q2.gpl");
   MappingQ<2> mapping(2);
   grid_out.write_gnuplot(triangulation, out2, &mapping);

   triangulation.refine_global(1);

   std::ofstream out3("grid_q2_ref.gpl");
   grid_out.write_gnuplot(triangulation, out3, &mapping);
}
