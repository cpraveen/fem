#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;

template <int dim, int spacedim>
void
print_statistics(const Triangulation<dim, spacedim> &tria, bool do_mg = false)
{
  deallog << "n_levels:                  " << tria.n_levels() << std::endl;
  deallog << "n_cells:                   " << tria.n_cells() << std::endl;
  deallog << "n_active_cells:            " << tria.n_active_cells()
          << std::endl;
  deallog << "has_hanging_nodes:         "
          << (tria.has_hanging_nodes() ? "true" : "false") << std::endl;

  if (do_mg)
    {
      for (auto level = 0u; level < tria.n_levels(); ++level)
        deallog << "n_cells on level=" << level << ":        "
                << tria.n_cells(level) << std::endl;

      for (auto level = 0u; level < tria.n_levels(); ++level)
        deallog << "n_active_cells on level=" << level << ": "
                << tria.n_active_cells(level) << std::endl;
    }

  deallog << std::endl;
}

template <int dim, int spacedim>
void
print_statistics(const DoFHandler<dim, spacedim> &dof_handler,
                 bool                             do_mg = false)
{
  deallog << "n_dofs:                             " << dof_handler.n_dofs()
          << std::endl;
  deallog << "n_locally_owned_dofs:               "
          << dof_handler.n_locally_owned_dofs() << std::endl;
  deallog << "has_hanging_nodes:                  "
          << (dof_handler.get_triangulation().has_hanging_nodes() ? "true" :
                                                                    "false")
          << std::endl;

  const auto n_levels = dof_handler.get_triangulation().n_levels();

  if (do_mg)
    {
      for (auto level = 0u; level < n_levels; ++level)
        deallog << "n_dofs on level=" << level << ":                  "
                << dof_handler.n_dofs(level) << std::endl;

      for (auto level = 0u; level < n_levels; ++level)
        deallog << "n_locally_owned_mg_dofs on level=" << level << ": "
                << dof_handler.locally_owned_mg_dofs(level).n_elements()
                << std::endl;
    }

  deallog << std::endl;
}

template <int dim>
void
test(MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria;

  GridIn<dim> grid_in;
  grid_in.attach_triangulation(basetria);
  std::ifstream gfile("grid.msh");
  AssertThrow(gfile.is_open(), ExcMessage("Grid file not found"));
  grid_in.read_msh(gfile);

  GridTools::partition_triangulation(
    Utilities::MPI::n_mpi_processes(comm), basetria);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  // test triangulation
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);

  // print statistics
  print_statistics(tria_pft);
  print_statistics(dof_handler);

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.build_patches (fe.degree);

  data_out.write_vtu_with_pvtu_record("./",
                                     "grid",
                                     0,
                                     comm,
                                     3);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  const MPI_Comm comm = MPI_COMM_WORLD;
  const auto rank = Utilities::MPI::this_mpi_process(comm);
  std::ofstream log("log_" + Utilities::int_to_string(rank,3) + ".txt");
  deallog.attach(log);

  {
    deallog.push("2d");
    test<2>(comm);
    deallog.pop();
  }

}
