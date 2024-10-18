// Gaussian bump

//------------------------------------------------------------------------------
namespace ProblemData
{
   const std::string name = "GAUSSIAN BUMP";
   const double xmin = 0.0;
   const double xmax = 1.0;
   const double ymin = 0.0;
   const double ymax = 1.0;
   const double final_time = 100.0;
   const bool periodic_x = false;
   const bool periodic_y = false;

   const double gamma = 1.4;
}

//------------------------------------------------------------------------------
template <int dim>
struct Problem : ProblemBase<dim>
{
   const double gasGam = ProblemData::gamma;

   const double mach_inf = 0.5;
   const double rho_inf = 1.0;
   const double vel_inf = 1.0;
   const double pre_inf = 1.0/(gasGam * mach_inf * mach_inf);
   const double pre_out = pre_inf;

   //---------------------------------------------------------------------------
   void initial_value(const Point<dim>& /*p*/,
                      Vector<double>&   u) const override
   {
      u[0] = rho_inf;
      u[1] = rho_inf * vel_inf;
      u[2] = 0.0;
      u[3] = pre_inf / (gasGam - 1.0) + 0.5 * rho_inf * pow(vel_inf,2);
   }

   //---------------------------------------------------------------------------
   void boundary_value(const int             boundary_id,
                       const Point<dim>&     /*p*/,
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
         case 1: // bottom
         case 3: // top
         {
            const double vn = vint * normal;
            const Tensor<1,dim> vout = vint - (2.0 * vn) * normal;
            Uout[0] = Uint[0];
            Uout[1] = Uint[0] * vout[0];
            Uout[2] = Uint[0] * vout[1];
            Uout[3] = Uint[3];
            break;
         }

         case 2: // outflow
         {
            Uout = Uint;
            const double Eout = pre_out / (gasGam - 1.0) + 0.5 * rho * vint.norm_square();
            Uout[3] = Eout;
            break;
         }

         case 4: // inflow
         {
            Uout[0] = rho_inf;
            Uout[1] = rho_inf * vel_inf;
            Uout[2] = 0.0;
            Uout[3] = pre_inf / (gasGam - 1.0) + 0.5 * rho_inf * pow(vel_inf,2);
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
      triangulation.set_all_manifold_ids_on_boundary(1, 1);
      const FunctionManifold<2, 2, 1> bump("x; 0.0625*exp(-25.0*x*x)", "x");
      triangulation.set_manifold(1, bump);
   }
};
