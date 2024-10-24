// Linear advection on circles

//------------------------------------------------------------------------------
namespace ProblemData
{
   const std::string name = "Advection on circles";
   const double xmin =-1.0;
   const double xmax = 1.0;
   const double ymin =-1.0;
   const double ymax = 1.0;
   const double final_time = 2.0 * M_PI;
   const bool periodic_x = false;
   const bool periodic_y = false;

   //---------------------------------------------------------------------------
   // Velocity field
   //---------------------------------------------------------------------------
   template <>
   void velocity(const Point<2>& p, Tensor<1,2>& v)
   {
      v[0] = -p[1];
      v[1] =  p[0];
   }
}

//------------------------------------------------------------------------------
template <int dim>
struct Problem : ProblemBase<dim>
{
   const double alpha = 50.0;
   const double x0 = 0.5;
   const double y0 = 0.0;

   void initial_value(const Point<dim>& p,
                      Vector<double>&   u) const override
   {
      const double x = p[0] - x0;
      const double y = p[1] - y0;
      u[0] = 1.0 + exp(-alpha * (x * x + y * y));
   }

   // Not needed if we have periodic bc
   void boundary_value(const int             /*boundary_id*/,
                       const Point<dim>&     /*p*/,
                       const double          /*t*/,
                       const Tensor<1, dim>& /*normal*/,
                       const Vector<double>& /*uin*/,
                       Vector<double>&       uout) const override
   {
      uout[0] = 1.0;
   }

   // See geo file or msh file for boundary ids
   // Here: 101 = inner circular boundary
   //       102 = outer circular boundary
   void set_manifolds(Triangulation<dim>& triangulation) const override
   {
      triangulation.set_all_manifold_ids(0);
      triangulation.set_all_manifold_ids_on_boundary(101, 1);
      triangulation.set_all_manifold_ids_on_boundary(102, 1);
      const Point<dim> center(0.0, 0.0);
      const SphericalManifold<dim> circle(center);
      triangulation.set_manifold(1, circle);
   }
};
