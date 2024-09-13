using namespace dealii;

namespace Problem
{

double xmin = 0.0;
double xmax = 1.0;
bool periodic = true;

// PDE data
double rho = 1.0;
double bulk = 1.0;

// IC data
const double x0 = 0.5;
const double beta = 100.0;

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
void
initial_value(const Point<1>& p,
              Vector<double>&   values)
{
   double x = p[0];
   values[0] = exp(-beta*pow(x-x0,2));
   values[1] = 0.0;
}

}
