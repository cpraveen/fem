// Isentropic vortex

//------------------------------------------------------------------------------
namespace ProblemData
{
   const std::string name = "ISENTROPIC VORTEX";
   const double xmin = -10.0;
   const double xmax =  10.0;
   const double ymin = -10.0;
   const double ymax =  10.0;
   const double final_time = 40.0 * sqrt(2);
   const bool periodic_x = true;
   const bool periodic_y = true;

   const double gamma = 1.4;
}

//------------------------------------------------------------------------------
template <int dim>
struct Problem : ProblemBase<dim>
{
   const double gamma = ProblemData::gamma;

   const double mach_inf = 0.5;
   const double theta = 45.0 * (M_PI / 180.0);
   const double u0 = mach_inf * cos(theta);
   const double v0 = mach_inf * sin(theta);
   const double beta = 5.0;
   const double x0 = 0.0, y0 = 0.0;
   const double a1 = 0.5 * beta / M_PI;
   const double a2 = 0.5 * (gamma - 1.0) * pow(a1, 2) / gamma;

   //---------------------------------------------------------------------------
   void initial_value(const Point<dim>& p,
                      Vector<double>&   u) const override
   {
      double r2 = pow(p[0] - x0, 2) + pow(p[1] - y0, 2);
      double rho = pow(1.0 - a2 * exp(1.0 - r2), 1.0 / (gamma - 1.0));
      double vex = u0 - a1 * (p[1] - y0) * exp(0.5 * (1.0 - r2));
      double vey = v0 + a1 * (p[0] - x0) * exp(0.5 * (1.0 - r2));
      double pre = pow(rho, gamma);

      u[0] = rho;
      u[1] = rho * vex;
      u[2] = rho * vey;
      u[3] = pre / (gamma - 1.0) + 0.5 * rho * (pow(vex,2) + pow(vey,2));
   }
};
