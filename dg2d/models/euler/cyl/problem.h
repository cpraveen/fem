//------------------------------------------------------------------------------
namespace ProblemData
{
   const std::string name = "CYLINDER";
   const double xmin = 0.0;
   const double xmax = 0.0;
   const double ymin = 0.0;
   const double ymax = 0.0;
   const double final_time = 100.0;
   const bool periodic_x = false;
   const bool periodic_y = false;

   const double gamma = 1.4;
}

//------------------------------------------------------------------------------
template <int dim>
struct Problem : ProblemBase<dim>
{
   const double gamma = ProblemData::gamma;
   const double mach = 0.3;
   const double rho_inf = 1.0;
   const double vel_inf = 1.0;
   const double pre_inf = 1.0/(gamma * mach * mach);

   //---------------------------------------------------------------------------
   void initial_value(const Point<dim>& /*p*/,
                      Vector<double>&   u) const override
   {
      const double rho = rho_inf;
      const double pre = pre_inf;
      Tensor<1,dim> vel;
      vel[0] = vel_inf;
      vel[1] = 0.0;

      PDE::prim2con(rho, vel, pre, u);
   }

   //---------------------------------------------------------------------------
   void boundary_value(const int             boundary_id,
                       const Point<dim>&     p,
                       const double          /*t*/,
                       const Tensor<1,dim>&  normal,
                       const Vector<double>& Uint,
                       Vector<double>&       Uout) const override
   {
      double rho, pre;
      Tensor<1,dim> vint;
      PDE::con2prim(Uint, rho, vint, pre);

      switch(boundary_id)
      {
         case 101: // cylinder
         {
            const double vn = vint * normal;
            const Tensor<1,dim> vout = vint - (2.0 * vn) * normal;
            Uout[0] = Uint[0];
            Uout[1] = Uint[0] * vout[0];
            Uout[2] = Uint[0] * vout[1];
            Uout[3] = Uint[3];
            break;
         }

         case 102: // farfield
         {
            initial_value(p, Uout);
            break;
         }

         default:
            DEAL_II_NOT_IMPLEMENTED();
      }
   }

   //---------------------------------------------------------------------------
   void set_manifolds(Triangulation<dim>& triangulation) const override
   {
      triangulation.set_all_manifold_ids(0);
      triangulation.set_all_manifold_ids_on_boundary(101, 1);
      const SphericalManifold<dim> cylinder(Point<dim>(0,0));
      triangulation.set_manifold(1, cylinder);
   }

};
