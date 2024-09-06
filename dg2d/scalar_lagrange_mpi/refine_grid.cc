template <int dim>
void refine_grid(Triangulation<dim> &triangulation, double r0, double r1)
{
   const Point<dim> center(0.0, 0.0);
   for (auto &cell : triangulation.active_cell_iterators())
   {
      for (const auto v : cell->vertex_indices())
      {
         const double distance_from_center =
             center.distance(cell->vertex(v));

         if (distance_from_center <= r1 && distance_from_center >= r0)
         {
            cell->set_refine_flag();
            break;
         }
      }
   }
   triangulation.execute_coarsening_and_refinement();
}
