#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/function_parser.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>

using namespace dealii;

const int dim = 2;

int main()
{
  int nelem = 10;
  const int degree = 3;
  std::cout << nelem << " x " << nelem << " cells\n";
  std::cout << "Mapping degree = " << degree << std::endl;
  Triangulation<dim> triangulation;
  GridGenerator::subdivided_hyper_cube(triangulation, nelem);
  {
    // Save initial grid
    std::ofstream out("grid_init.eps");
    GridOut().write_eps(triangulation, out);
    std::cout << "Wrote grid_init.eps\n";
  }

  // Define function to transform grid
  std::map<std::string, double> constants;
  constants["pi"] = numbers::PI;
  std::string variables = "x,y";
  std::vector<std::string> expressions(2);
  expressions[0] = "x + 0.25*x*(1-x)*sin(2*pi*y)";
  expressions[1] = "y + 0.25*y*(1-y)*sin(2*pi*x)";
  // function parser with 2 variables and 2 components
  FunctionParser<dim> Transform(dim);
  // And populate it with the newly created objects.
  Transform.initialize(variables,
                       expressions,
                       constants);

  // Create MappingFEField using the above transformation
  const FE_Q<dim> feq(degree);
  const FESystem<dim> fesystem(feq, dim);
  DoFHandler<dim> dhq(triangulation);
  dhq.distribute_dofs(fesystem);
  Vector<double> eulerq(dhq.n_dofs());
  // Fills the euler vector with information from the Triangulation
  VectorTools::interpolate(dhq, Transform, eulerq);
  MappingFEField<dim> map(dhq, eulerq);
  {
    // This still shows Cartesian since we have not attached
    // any manifold to the grid.
    std::ofstream out("grid_mapped.eps");
    GridOut().write_eps(triangulation, out, &map);
    std::cout << "Wrote grid_mapped.eps\n";
  }
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dhq);
    std::vector<std::string> names{"dx","dy"};
    data_out.add_data_vector(eulerq, names);
    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);
    typename DataOut<dim>::CurvedCellRegion curved_region = DataOut<dim>::curved_inner_cells;
    data_out.build_patches(map, degree, curved_region);

    std::ofstream out("grid.vtu");
    data_out.write_vtu(out);
    std::cout << "Wrote grid.vtu  ==> open in paraview\n";
  }
}
