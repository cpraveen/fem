// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a serial triangulation with periodic face in x-direction and copy it.
//
// Based on code taken from
// https://github.com/dealii/dealii/blob/master/tests/fullydistributed_grids/copy_serial_tria_04.cc

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

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
test(const int n_refinements, const int n_subdivisions, MPI_Comm comm)
{
  const double left  = 0;
  const double right = 1;

  auto add_periodicity = [&](dealii::Triangulation<dim> &tria) {
    std::vector<
      GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
         periodic_faces;
    auto cell = tria.begin();
    auto endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int face_number : GeometryInfo<dim>::face_indices())
        if (std::fabs(cell->face(face_number)->center()[0] - left) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(1);
        else if (std::fabs(cell->face(face_number)->center()[0] - right) <
                 1e-12)
          cell->face(face_number)->set_all_boundary_ids(2);
        else if (dim >= 2 &&
                 std::fabs(cell->face(face_number)->center()[1] - left) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(3);
        else if (dim >= 2 && std::fabs(cell->face(face_number)->center()[1] -
                                       right) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(4);
        else if (dim >= 3 &&
                 std::fabs(cell->face(face_number)->center()[2] - left) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(5);
        else if (dim >= 3 && std::fabs(cell->face(face_number)->center()[2] -
                                       right) < 1e-12)
          cell->face(face_number)->set_all_boundary_ids(6);

    if (dim >= 1)
      GridTools::collect_periodic_faces(tria, 1, 2, 0, periodic_faces);
    if (dim >= 2)
      GridTools::collect_periodic_faces(tria, 3, 4, 1, periodic_faces);
    if (dim >= 3)
      GridTools::collect_periodic_faces(tria, 5, 6, 2, periodic_faces);

    tria.add_periodicity(periodic_faces);
  };

  // create serial triangulation
  Triangulation<dim> basetria;
  GridGenerator::subdivided_hyper_cube(basetria, n_subdivisions);
  // new: add periodicity on serial mesh
  add_periodicity(basetria);
  basetria.refine_global(n_refinements);

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  // new: add periodicity on fullydistributed mesh (!!!)
  add_periodicity(tria_pft);

  // test triangulation
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);

  // print statistics
  print_statistics(tria_pft);
  print_statistics(dof_handler);

  if(dim == 2)
  {
     DataOut<dim> data_out;
     data_out.attach_dof_handler (dof_handler);
     data_out.build_patches (fe.degree);

     data_out.write_vtu_with_pvtu_record("./",
                                        "grid",
                                        0,
                                        comm,
                                        3);
  }
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
    deallog.push("1d");
    const int n_refinements  = 4;
    const int n_subdivisions = 8;
    test<1>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }

  {
    deallog.push("2d");
    const int n_refinements  = 3;
    const int n_subdivisions = 8;
    test<2>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }

  {
    deallog.push("3d");
    const int n_refinements  = 3;
    const int n_subdivisions = 4;
    test<3>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }
}
