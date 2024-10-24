//------------------------------------------------------------------------------
namespace ProblemData
{
   const std::string name = "FORWARD STEP";
   const double xmin = 0.0;
   const double xmax = 0.0;
   const double ymin = 0.0;
   const double ymax = 0.0;
   const double final_time = 3.5;
   const bool periodic_x = false;
   const bool periodic_y = false;

   const double gamma = 1.4;
}

//------------------------------------------------------------------------------
template <int dim>
struct Problem : ProblemBase<dim>
{
   const double gamma = ProblemData::gamma;
   const double mach = 3.0;
   const double rho_inf = gamma;
   const double vel_inf = mach;
   const double pre_inf = 1.0;

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
         case 1: // inlet
         {
            initial_value(p, Uout);
            break;
         }

         case 2: // bottom and top
         {
            const double vn = vint * normal;
            const Tensor<1,dim> vout = vint - (2.0 * vn) * normal;
            Uout[0] = Uint[0];
            Uout[1] = Uint[0] * vout[0];
            Uout[2] = Uint[0] * vout[1];
            Uout[3] = Uint[3];
            break;
         }

         case 3: // outlet
         {
            Uout = Uint;
            break;
         }

         default:
            DEAL_II_NOT_IMPLEMENTED();
      }
   }

};
