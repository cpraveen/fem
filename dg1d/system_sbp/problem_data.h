using namespace dealii;

namespace Problem
{

   double xmin = 0.0;
   double xmax = 1.0;
   double final_time = 0.2;
   bool periodic = false;

// PDE data
   double gamma = 1.4;

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
   void
   initial_value(const Point<1>& p,
                 Vector<double>& values)
   {
      double x = p[0];
      double rho, vel, pre;
      if(x < 0.5)
      {
         rho = 1.0;
         vel = 0.0;
         pre = 1.0;
      }
      else
      {
         rho = 0.125;
         vel = 0.0;
         pre = 0.1;
      }

      values[0] = rho;
      values[1] = rho * vel;
      values[2] = pre / (gamma - 1.0) + 0.5 * rho * pow(vel, 2);
   }

//------------------------------------------------------------------------------
// Boundary condition
//------------------------------------------------------------------------------
   void
   boundary_value(const int             /*id*/,
                  const double          /*t*/,
                  const Vector<double>& ul,
                  Vector<double>&       ur)
   {
      // Neumann bc, we assume waves do not reach boundary
      ur = ul;
   }

}
