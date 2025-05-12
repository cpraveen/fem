// Performance of high-order Godunov-type methods in simulations of
// astrophysical low Mach number flows

//------------------------------------------------------------------------------
namespace ProblemData
{
   const std::string name = "KELVIN-HELMHOLTZ";
   const double xmin = 0.0;
   const double xmax = 2.0;
   const double ymin = -0.5;
   const double ymax =  0.5;
   const double M0 = 0.01;
   const double final_time = 0.8 / M0;
   const bool periodic_x = true;
   const bool periodic_y = true;

   const double gamma = 1.4;
}

//------------------------------------------------------------------------------
template <int dim>
struct Problem : ProblemBase<dim>
{
   const double gamma = ProblemData::gamma;

   const double M0 = ProblemData::M0;

   double eta(const double y) const
   {
      double X;

      if(y > -9.0/32.0 && y < -7.0/32.0)
      {
         X = 0.5 * (1.0 + sin(16 * M_PI * (y + 0.25)));
      }
      else if(y >= -7.0/32.0 && y <= 7.0/32.0)
      {
         X = 1.0;
      }
      else if(y > 7.0/32.0 && y < 9.0/32.0)
      {
         X = 0.5 * (1.0 - sin(16 * M_PI * (y - 0.25)));
      }
      else
      {
         X = 0.0;
      }

      return X;
   }

   //---------------------------------------------------------------------------
   void initial_value(const Point<dim>& p,
                      Vector<double>&   u) const override
   {
      const double x = p[0];
      const double y = p[1];

      const double rho = gamma + 1.0e-3 * (1.0 - 2.0 * eta(y));
      const double vex = M0 * (1.0 - 2.0 * eta(y));
      const double vey = (M0/10) * sin(2.0 * M_PI * x);
      const double pre = 1.0;

      u[0] = rho;
      u[1] = rho * vex;
      u[2] = rho * vey;
      u[3] = pre / (gamma - 1.0) + 0.5 * rho * (pow(vex,2) + pow(vey,2));
   }
};
