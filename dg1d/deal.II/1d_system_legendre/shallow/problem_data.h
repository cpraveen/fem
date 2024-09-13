using namespace dealii;

namespace Problem
{

   double xmin = -5.0;
   double xmax =  5.0;
   double final_time = 2.0;
   bool periodic = false;

// PDE data
   double g = 1.0;

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
   void
   initial_value(const Point<1>& p,
                 Vector<double>& values)
   {
      double x = p[0];
      if(x < 0.0)
      {
         values[0] = 3.0;
         values[1] = 0.0;
      }
      else
      {
         values[0] = 1.0;
         values[1] = 0.0;
      }
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
