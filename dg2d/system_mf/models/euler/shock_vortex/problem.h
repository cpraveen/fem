//------------------------------------------------------------------------------
namespace ProblemData
{
   const std::string name = "SHOCK VORTEX";
   const double xmin = 0.0;
   const double xmax = 1.0;
   const double ymin = 0.0;
   const double ymax = 1.0;
   const double final_time = 0.35;
   const bool periodic_x = false;
   const bool periodic_y = false;

   const double gamma = 1.4;
}

//------------------------------------------------------------------------------
template <int dim>
struct Problem : ProblemBase<dim>
{
   const double gasGam = ProblemData::gamma;
   const double gasR = 1.0;

   // Vortex Perturbation Variables
   const double eps = 0.3;
   const double rc = 0.05;
   const double xc = 0.25;
   const double yc = 0.50;
   const double beta = 0.204;

   // Left State
   const double u_l   = sqrt(gasGam);
   const double v_l   = 0.0;
   const double rho_l = 1.0;
   const double p_l   = 1.0;
   const double T_l   = p_l/(gasR*rho_l);

   // Right State
   const double p_r   = 1.3;
   const double rho_r =   rho_l*(gasGam-1.0 + (gasGam + 1.0)*p_r) 
                        / (gasGam + 1.0 + (gasGam - 1.0)*p_r);
   const double u_r   = u_l + (1.0-p_r)*sqrt(2.0/(gasGam - 1.0 + p_r*(gasGam + 1.0)));

   // Location of jump
   const double xdia  = 0.5;

   //---------------------------------------------------------------------------
   void initial_value(const Point<dim>& p,
                      Vector<double>&   u) const override
   {
      const double x = p[0];
      const double y = p[1];
      double rho, pre;
      Tensor<1,dim> vel;

      if(x < xdia)
      {
         double r = (pow(x - xc,2) + pow(y - yc,2))/pow(rc,2);
         double du =  eps*((y - yc)/rc)*exp(beta*(1.0-r));
         double dv = -eps*((x - xc)/rc)*exp(beta*(1.0-r));
         double dt = -eps*eps*((gasGam-1.0)/(4.0*beta*gasGam))*exp(2.0*beta*(1.0-r));
         double dT_l = T_l + dt;
         rho = pow(dT_l, 1.0/(gasGam - 1.0));
         pre = pow(rho_l,gasGam);
         vel[0] = u_l + du;
         vel[1] = v_l + dv;
      }
      else
      {
         vel[0] = u_r;
         vel[1] = 0.0;
         rho = rho_r;
         pre = p_r;
      }

      PDE::prim2con(rho, vel, pre, u);
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
         case 2: // bottom
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

         case 0: // inflow
         case 1: // outflow
         {
            Uout = Uint;
            break;
         }

         default:
            DEAL_II_NOT_IMPLEMENTED();
      }
   }

};
